# ============================================================
# SPUDD: Stochastic Planning Using Decision Diagrams
# Full ADD representation — transitions are ADDs, not matrices
#
# State space : X in {0, 50, 100, 150}, Y in {0, 50, 100, 150}, Z in {0, 50, 100, 150}
#               Total: 4^3 = 64 states
# Actions     : a1, a2, a3
#
# Conditional independence structure:
#   P(X' | X, a)  — X' depends only on X (and action)
#   P(Y' | Y, a)  — Y' depends only on Y (and action)
#   P(Z' | Z, a)  — Z' depends only on Z (and action)
#
# Transition ADDs use "primed" variable names (X_p, Y_p, Z_p)
# to distinguish next-state variables from current-state variables.
# The regression step multiplies T_ADD(Xi, Xi_p) against V(Xi_p, ...)
# and sums out Xi_p — entirely symbolically on ADDs.
# ============================================================


# ------------------------------------------------------------
# 1. Multi-Valued ADD Representation
# ------------------------------------------------------------

N_STATES   <- 4L
STATE_VALS <- c(0, 50, 100, 150)
VARS       <- c("X",   "Y",   "Z")
VARS_P     <- c("X_p", "Y_p", "Z_p")  # primed (next-state) variable names

###
# functions for creating leaves and nodes
###

# Canonical leaf cache - ensures only one leaf per unique value
leaf_cache <- new.env(hash = TRUE, parent = emptyenv())
# node-level deduplication
node_cache <- new.env(hash = TRUE, parent = emptyenv())

leaf <- function(v) {
  key <- as.character(v)
  if (exists(key, envir = leaf_cache))
    return(get(key, envir = leaf_cache))
  l <- list(type = "leaf", value = v)
  assign(key, l, envir = leaf_cache)
  l
}

node <- function(var, children) {
  stopifnot(length(children) == N_STATES)
  
  # Reduction rule: all children identical -> return child directly
  if (length(unique(lapply(children, identity))) == 1L)
    return(children[[1]])
  
  # Deduplicate identical subtrees via cache
  key <- paste0(var, ":", paste(sapply(children, function(ch) {
    if (ch$type == "leaf") paste0("L", ch$value)
    else ch$.id
  }), collapse = ","))
  
  if (exists(key, envir = node_cache))
    return(get(key, envir = node_cache))
  
  n <- list(type = "node", var = var, children = children, .id = key)
  assign(key, n, envir = node_cache)
  n
}

build_add_from_fn <- function(f, var_order = VARS) {
  rm(list = ls(leaf_cache), envir = leaf_cache)
  rm(list = ls(node_cache), envir = node_cache)
  
  .build <- function(depth, partial) {
    if (depth > length(var_order)) return(leaf(f(partial)))
    var      <- var_order[[depth]]
    children <- setNames(
      lapply(STATE_VALS, function(v)
        .build(depth + 1L, c(partial, setNames(list(v), var)))),
      as.character(STATE_VALS))
    node(var, children)
  }
  .build(1L, list())
}


# ------------------------------------------------------------
# 3. Reward Function R(X, Y, Z, a)
# ------------------------------------------------------------

action_cost <- c(a1 = 10, a2 = 10, a3 = 00)

action_names <- c("a1", "a2", "a3")

reward_fn <- function(x, y, z, a) {
  
  R_S <- x + y + z
  R_a <- action_cost[[a]]
  
  return(R_S + R_a)
}

reward_adds <- setNames(
  lapply(action_names, function(a)
    build_add_from_fn(function(s) reward_fn(s$X, s$Y, s$Z, a))),
  action_names)



















node_vec <- function(var, child_list)
  node(var, setNames(child_list, as.character(STATE_VALS)))

# Evaluate an ADD at a variable assignment (may be partial)
eval_add <- function(add, assignment) {
  if (add$type == "leaf") return(add$value)
  key <- as.character(assignment[[add$var]])
  if (is.null(key) || is.na(key))
    stop(sprintf("Variable '%s' missing from assignment", add$var))
  eval_add(add$children[[key]], assignment)
}

# Pointwise binary operation on two ADDs.
# Variable ordering: current vars (X,Y,Z) precede primed vars (X_p,Y_p,Z_p).
# apply_op expands whichever variable it encounters first; the other ADD is
# passed through unchanged if it doesn't branch on that variable.
apply_op <- function(add1, add2, op) {
  if (add1$type == "leaf" && add2$type == "leaf")
    return(leaf(op(add1$value, add2$value)))
  
  var <- if (add1$type == "node") add1$var else add2$var
  
  get_child <- function(add, v) {
    if (add$type == "leaf")        return(add)
    if (add$var  != var)           return(add)
    add$children[[as.character(v)]]
  }
  
  node_vec(var, lapply(STATE_VALS, function(v)
    apply_op(get_child(add1, v), get_child(add2, v), op)))
}

scale_add <- function(add, scalar) {
  if (add$type == "leaf") return(leaf(add$value * scalar))
  node(add$var,
       setNames(lapply(add$children, scale_add, scalar = scalar),
                names(add$children)))
}

# Restrict: set variable `var` to integer `val`
restrict <- function(add, var, val) {
  if (add$type == "leaf") return(add)
  if (add$var  == var)    return(add$children[[as.character(val)]])
  node(add$var,
       setNames(lapply(add$children, restrict, var = var, val = val),
                names(add$children)))
}

# Sum out a variable: replace node with sum of all its children
# (used to marginalise over a primed variable after weighting by T)
sum_out <- function(add, var) {
  if (add$type == "leaf") return(add)   # variable absent — nothing to sum
  if (add$var  == var) {
    # Sum all N_STATES children
    return(Reduce(function(a, b) apply_op(a, b, `+`), add$children))
  }
  node(add$var,
       setNames(lapply(add$children, sum_out, var = var),
                names(add$children)))
}

# Build an ADD from a scalar function f(assignment)
build_add_from_fn <- function(f, var_order = VARS) {
  .build <- function(depth, partial) {
    if (depth > length(var_order)) return(leaf(f(partial)))
    var      <- var_order[[depth]]
    children <- setNames(
      lapply(STATE_VALS, function(v)
        .build(depth + 1L, c(partial, setNames(list(v), var)))),
      as.character(STATE_VALS))
    node(var, children)
  }
  .build(1L, list())
}


# ------------------------------------------------------------
# 2. Transition ADDs
#
#    For each (action, variable Xi), we build a two-level ADD:
#
#      T_add  branches first on Xi  (current state)
#             then  on Xi_p (next state, the "primed" variable)
#             leaf  = P(Xi_p | Xi, action)   [a scalar probability]
#
#    Structure:
#      node(Xi,
#        "0" -> node(Xi_p, "0"->leaf(p00), "1"->leaf(p01), ...),
#        "1" -> node(Xi_p, "0"->leaf(p10), ...),
#        ...
#      )
#
#    This is the natural ADD encoding of a stochastic matrix.
#    When parent sets are non-trivial (Xi_p depends on multiple
#    current variables), the Xi level would be replaced by a
#    sub-tree over the parent set — the regression code below
#    handles that automatically via apply_op / sum_out.
# ------------------------------------------------------------

# Build the underlying stochastic matrix (helper — internal use only)
.make_stoch_matrix <- function(shift = 0, spread = 0.15) {
  M <- matrix(0, N_STATES, N_STATES)
  for (xi in seq_along(STATE_VALS)) {
    target <- STATE_VALS[xi] + shift
    for (xp in seq_along(STATE_VALS)) {
      d <- abs(STATE_VALS[xp] - target)
      M[xi, xp] <- (1 - spread) * (d == 0) + (spread / 2) * (d == 1)
    }
    M[xi, ] <- M[xi, ] / sum(M[xi, ])   # renormalise at boundaries
  }
  M
}

# Convert a stochastic matrix to a two-level ADD T(Xi, Xi_p)
# var_cur  : current-state variable name  e.g. "X"
# var_next : primed variable name          e.g. "X_p"
make_trans_add <- function(shift = 0, spread = 0.15,
                           var_cur = "X", var_next = "X_p") {
  M <- .make_stoch_matrix(shift, spread)
  
  # For each current value xi, build the inner ADD over Xi_p
  outer_children <- setNames(
    lapply(seq_along(STATE_VALS), function(i) { # 'i' is the index (1, 2, 3...)
      xi <- STATE_VALS[i] 
      
      inner_children <- setNames(
        lapply(seq_along(STATE_VALS), function(j) { # 'j' is the index
          xp <- STATE_VALS[j]
          
          # Use the indices i and j to access the matrix
          leaf(M[i, j]) 
        }),
        as.character(STATE_VALS))
      
      node(var_next, inner_children)
    }),
    as.character(STATE_VALS))
  
  node(var_cur, outer_children)
}

# Build all transition ADDs
# transitions[[action]][[var]] = ADD over (var, var_p)
transitions <- list(
  a1 = list(
    X = make_trans_add(shift =  1, spread = 0.20, var_cur = "X", var_next = "X_p"),
    Y = make_trans_add(shift =  1, spread = 0.15, var_cur = "Y", var_next = "Y_p"),
    Z = make_trans_add(shift =  1, spread = 0.25, var_cur = "Z", var_next = "Z_p")
  ),
  a2 = list(
    X = make_trans_add(shift = -1, spread = 0.20, var_cur = "X", var_next = "X_p"),
    Y = make_trans_add(shift = -1, spread = 0.15, var_cur = "Y", var_next = "Y_p"),
    Z = make_trans_add(shift = -1, spread = 0.25, var_cur = "Z", var_next = "Z_p")
  ),
  a3 = list(
    X = make_trans_add(shift =  0, spread = 0.30, var_cur = "X", var_next = "X_p"),
    Y = make_trans_add(shift =  0, spread = 0.20, var_cur = "Y", var_next = "Y_p"),
    Z = make_trans_add(shift =  0, spread = 0.10, var_cur = "Z", var_next = "Z_p")
  )
)

# Sanity check: for each (xi, action, var), probabilities sum to 1
for (a in names(transitions)) {
  for (i in seq_along(VARS)) {
    T_add  <- transitions[[a]][[VARS[i]]]
    vp     <- VARS_P[i]
    for (xi in STATE_VALS) {
      total <- sum(sapply(STATE_VALS, function(xp)
        eval_add(T_add, setNames(list(xi, xp), c(VARS[i], vp)))))
      stopifnot(abs(total - 1) < 1e-10)
    }
  }
}
cat("Transition ADD sanity check passed.\n")


# ------------------------------------------------------------
# 3. Reward Function R(X, Y, Z, a) — unchanged
# ------------------------------------------------------------

action_cost <- c(a1 = 10, a2 = 10, a3 = 00)

reward_fn <- function(x, y, z, a) {
  
  R_S <- x + y + z
  R_a <- action_cost[[a]]
  
  return(R_S + R_a)
}

reward_adds <- setNames(
  lapply(names(transitions), function(a)
    build_add_from_fn(function(s) reward_fn(s$X, s$Y, s$Z, a))),
  names(transitions))


# ------------------------------------------------------------
# 4. ADD-Based Regression (Bellman Expectation)
#
#    Goal: given V(X, Y, Z) and T_add(Xi, Xi_p), compute
#          E_{Xi_p}[V] as a function of (X, Y, Z) with Xi
#          replaced by Xi_p's expectation.
#
#    Steps for variable Xi (e.g. X / X_p):
#
#    (a) Rename Xi -> Xi_p in V:
#          V_p = rename_var(V, "X", "X_p")
#        V_p now branches on X_p instead of X, so it can be
#        multiplied against T_add which also branches on X_p.
#
#    (b) Multiply:
#          weighted = apply_op(T_add, V_p, `*`)
#        This ADD branches on Xi (current) and Xi_p (next),
#        with leaves = T(xi, xi_p) * V(xi_p, y, z).
#
#    (c) Sum out Xi_p:
#          expectation = sum_out(weighted, "X_p")
#        Leaves are now sum_{xi_p} T(xi,xi_p)*V(xi_p,y,z)
#        = E_{X'}[V | X=xi].
#        The result branches on X (current) as desired.
#
#    This is purely symbolic — no explicit state enumeration.
# ------------------------------------------------------------

# Rename a variable throughout an ADD (renames var_old -> var_new at every node)
rename_var <- function(add, var_old, var_new) {
  if (add$type == "leaf") return(add)
  new_var      <- if (add$var == var_old) var_new else add$var
  new_children <- setNames(
    lapply(add$children, rename_var, var_old = var_old, var_new = var_new),
    names(add$children))
  node(new_var, new_children)
}

# Regress V through the transition ADD for one variable
# var_cur:  current variable name  (e.g. "X")
# var_next: primed variable name   (e.g. "X_p")
# T_add:    transition ADD over (var_cur, var_next)
regress_variable_add <- function(V, var_cur, var_next, T_add) {
  # Step (a): relabel Xi -> Xi_p in V so next-state var aligns with T_add
  V_p <- rename_var(V, var_cur, var_next)
  
  # Step (b): pointwise multiply T(Xi, Xi_p) * V(Xi_p, other_vars)
  weighted <- apply_op(T_add, V_p, `*`)
  
  # Step (c): marginalise out Xi_p -> result is E[V | Xi = xi]
  sum_out(weighted, var_next)
}

# Apply regression for all three variables (order matters for ADD structure
# but not for the numeric result under conditional independence)
regress_value <- function(V, trans_action) {
  current <- V
  for (i in seq_along(VARS))
    current <- regress_variable_add(current,
                                    var_cur  = VARS[i],
                                    var_next = VARS_P[i],
                                    T_add    = trans_action[[VARS[i]]])
  current
}


# ------------------------------------------------------------
# 5. Bellman Backup and SPUDD Value Iteration
# ------------------------------------------------------------

bellman_backup <- function(V, reward_adds, transitions, gamma) {
  Q_adds <- mapply(function(trans_a, R_a) {
    regressed <- regress_value(V, trans_a)
    apply_op(R_a, scale_add(regressed, gamma), `+`)
  }, transitions, reward_adds, SIMPLIFY = FALSE)
  
  Reduce(function(a, b) apply_op(a, b, max), Q_adds)
}

spudd <- function(reward_adds, transitions, gamma = 0.9,
                  max_iter = 100, tol = 1e-4) {
  
  all_states <- expand.grid(X = STATE_VALS, Y = STATE_VALS, Z = STATE_VALS)
  
  V     <- leaf(0)
  diffs <- numeric(max_iter)
  
  cat("\nStarting SPUDD Value Iteration (ADD transitions)\n")
  cat(sprintf("  Variables : %s (each with %d states)\n",
              paste(VARS, collapse = ", "), N_STATES))
  cat(sprintf("  State space: %d^%d = %d states\n",
              N_STATES, length(VARS), N_STATES^length(VARS)))
  cat(sprintf("  Actions   : %s\n", paste(names(transitions), collapse = ", ")))
  cat(sprintf("  gamma=%.2f  max_iter=%d  tol=%.1e\n\n",
              gamma, max_iter, tol))
  
  for (iter in seq_len(max_iter)) {
    V_new <- bellman_backup(V, reward_adds, transitions, gamma)
    
    max_diff <- max(sapply(seq_len(nrow(all_states)), function(i) {
      s <- as.list(all_states[i, ])
      abs(eval_add(V_new, s) - eval_add(V, s))
    }))
    
    diffs[iter] <- max_diff
    cat(sprintf("  Iter %3d | max delta = %.6f\n", iter, max_diff))
    
    V <- V_new
    if (max_diff < tol) {
      cat(sprintf("\nConverged after %d iterations.\n", iter))
      break
    }
  }
  
  list(V = V, diffs = diffs[1:iter])
}


# ------------------------------------------------------------
# 6. Policy Extraction
# ------------------------------------------------------------

extract_policy <- function(V, reward_adds, transitions, gamma) {
  all_states   <- expand.grid(X = STATE_VALS, Y = STATE_VALS, Z = STATE_VALS)
  action_names <- names(transitions)
  
  Q_adds <- setNames(
    mapply(function(trans_a, R_a)
      apply_op(R_a, scale_add(regress_value(V, trans_a), gamma), `+`),
      transitions, reward_adds, SIMPLIFY = FALSE),
    action_names)
  
  policy <- apply(all_states, 1, function(s_vec) {
    s    <- as.list(s_vec)
    vals <- sapply(action_names, function(a) eval_add(Q_adds[[a]], s))
    action_names[which.max(vals)]
  })
  
  cbind(all_states, policy = policy)
}


# ------------------------------------------------------------
# 7. Inspect a Transition ADD (diagnostic helper)
#
#    Prints P(Xi_p | Xi = xi_val) for every next-state value —
#    confirms the ADD encodes the same distribution as the
#    equivalent stochastic matrix.
# ------------------------------------------------------------

inspect_trans_add <- function(T_add, var_cur, var_next, xi_val) {
  cat(sprintf("  P(%s' | %s=%d): ", var_cur, var_cur, xi_val))
  probs <- sapply(STATE_VALS, function(xp)
    eval_add(T_add, setNames(list(xi_val, xp), c(var_cur, var_next))))
  cat(round(probs, 4), "\n")
}

cat("\nSample transition ADD rows for action a1:\n")
for (xi in STATE_VALS)
  inspect_trans_add(transitions$a1$X, "X", "X_p", xi)


# ------------------------------------------------------------
# 8. Run, Report, and Plot
# ------------------------------------------------------------

result <- spudd(reward_adds, transitions, gamma = 0.9,
                max_iter = 100, tol = 1e-4)
V_star <- result$V

cat("\n--- Optimal Value Function V*(X, Y, Z) [first 20 rows] ---\n")
all_states        <- expand.grid(X = STATE_VALS, Y = STATE_VALS, Z = STATE_VALS)
all_states$V_star <- sapply(seq_len(nrow(all_states)), function(i)
  eval_add(V_star, as.list(all_states[i, ])))
print(head(all_states, 20))

cat("\n--- Optimal Policy (first 20 rows) ---\n")
policy_table <- extract_policy(V_star, reward_adds, transitions, gamma = 0.9)
print(head(policy_table, 20))

cat("\n--- Policy Action Frequencies ---\n")
print(table(policy_table$policy))

cat("\n--- V*(X, Y, Z=2) slice ---\n")
slice_z2 <- subset(all_states, Z == 2, select = c(X, Y, V_star))
V_matrix <- matrix(slice_z2$V_star,
                   nrow = N_STATES, ncol = N_STATES,
                   dimnames = list(paste0("X=", STATE_VALS),
                                   paste0("Y=", STATE_VALS)))
print(round(V_matrix, 3))

cat("\n--- Convergence residuals ---\n")
print(round(result$diffs, 6))

if (requireNamespace("graphics", quietly = TRUE)) {
  plot(result$diffs, type = "b", pch = 19, col = "steelblue",
       xlab = "Iteration", ylab = "Max |V_new - V_old|",
       main = "SPUDD Convergence — ADD Transitions",
       log  = "y")
  grid()
}

