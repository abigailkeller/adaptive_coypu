import operator
from functools import reduce
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from utils_ADD import *
import numpy as np

########################################
# helper functions for SPUDD algorithm #
########################################

def restrict(mgr, root, subsystem, val):
    """
    Fix subsystem to `val` in ADD, returning reduced ADD.
    Example: restrict(mgr, root, "X", 50) returns an ADD where P(X_p | X=50, a) for all values of X_p.
    """
    cache = {}
    def _restrict(n):
        if id(n) in cache:
            return cache[id(n)]
        if isinstance(n, Leaf):
            return n
        if n.var == subsystem:
            idx = mgr.state_vals.index(val)
            result = _restrict(n.children[idx])
        else:
            result = mgr.node(n.var, tuple(_restrict(c) for c in n.children))
        cache[id(n)] = result
        return result
    return _restrict(root)

def sum_out(mgr, root, subsystem):
    """
    Sum ADD over all values of a subsystem, marginalising it out.
    Example: if subsystem = "X_p", this computes sum_{X'} P(X_p | X, a) * V(X_p, Y_p, Z_P), returning an ADD that no longer contains X_p    
    """
    result = mgr.leaf(0.0)
    for val in mgr.state_vals:
        branch = restrict(mgr, root, subsystem, val)
        result = mgr.apply(operator.add, result, branch)
    return result

def add_max(mgr, a, b):
    """
    Return an ADD where each leaf is the maximum value of the leaves in input ADD; used to compute V(s) = max_a Q(s, a) across all actions.
    Example: if ADD a represents Q(s, a1) for state s and action a1, and ADD b represents Q(s, a2) for state s and action a2, 
             this function returns max(Q(s, a1), Q(s, a2))
    """
    return mgr.apply(lambda x, y: max(x, y), a, b)

def scalar_mul(mgr, root, scalar):
    """
    Multiply every leaf of ADD by a constant scalar; used to apply the discount factor to the expected future value.
    Example: scalar_mul(mgr, root, 0.95) will return 0.95 * E[V(s') | s, a]
    """
    cache = {}
    def _mul(n):
        if id(n) in cache:
            return cache[id(n)]
        if isinstance(n, Leaf):
            result = mgr.leaf(n.value * scalar)
        else:
            result = mgr.node(n.var, tuple(_mul(c) for c in n.children))
        cache[id(n)] = result
        return result
    return _mul(root)

def max_diff(mgr, a, b):
    """
    Computes the maximum absolute difference between two ADDs, used to check whether value iteration has converged.
    """
    diff = mgr.apply(lambda x, y: abs(x - y), a, b)
    vals = []
    def _collect(n):
        if isinstance(n, Leaf):
            vals.append(n.value)
        else:
            for c in n.children:
                _collect(c)
    _collect(diff)
    return max(vals)

def regress(mgr, V, subsystems, transition_adds_a, discount):
    """
    Compute expected discounted future value for action a.
    Compute Q(s, a) = discount * sum_{s'} P(s'|s,a) * V(s')
    
    For independent subsystems, joint transition factors:
    P(s'|s,a) = P(X'|X,a) * P(Y'|Y,a) * P(Z'|Z,a)
    
    We compute sum_{X',Y',Z'} P(X'|X,a)*P(Y'|Y,a)*P(Z'|Z,a)*V(X',Y',Z')
    by summing out primed variables one at a time.
    """
    # Start with V(X', Y', Z') and work through each subsystem
    result = V

    # loop through one subsystem at a time
    for subsystem in subsystems:

        var = subsystem
        var_p = subsystem + "_p"

        # get transition ADD for this subsystem
        T_var = transition_adds_a[var]  # P(var_p | var, a)
        
        # Rename V's var to var_p by building a new ADD
        # Then multiply P(var_p | var) * V(..., var_p, ...)
        # and sum out var_p
        
        # Multiply transition by current result
        # P(X_p | X, a) * V(X', Y', Z')
        product = mgr.apply(operator.mul, T_var, result)
        
        # Sum out the primed variable
        result = sum_out(mgr, product, var_p)
    
    return scalar_mul(mgr, result, discount)

def rename_to_primed(mgr, V, subsystems):
    """Return copy of V with subsystem variables replaced by primed versions."""
    cache = {}
    def _rename(n):
        if id(n) in cache:
            return cache[id(n)]
        if isinstance(n, Leaf):
            return n
        new_var = n.var + "_p" if n.var in subsystems else n.var
        result = mgr.node(new_var, tuple(_rename(c) for c in n.children))
        cache[id(n)] = result
        return result
    return _rename(V)

    
###################################
# SPUDD value iteration algorithm #
###################################

def spudd(mgr, subsystems, reward_adds, transition_adds, discount, tolerance, max_iter=1000):
    """
    SPUDD value iteration.
    
    Args:
        mgr:              ADD manager
        reward_adds:      dict {action_name -> reward ADD over states}
        transition_adds:  dict {action_idx -> {var -> transition ADD}}
        discount:         float, discount factor
        tolerance:        float, convergence threshold
        max_iter:         int, maximum iterations, default = 1000
    
    Returns:
        V:       value function ADD over states
        policy:  policy ADD over states (action index at each state)
    """
    # Map action names to indices
    action_names = list(reward_adds.keys())  # ['a1', 'a2', 'a3']
    
    # Initialise V = 0
    V = mgr.leaf(0.0)
    
    for iteration in range(max_iter):
        V_primed = rename_to_primed(mgr, V, subsystems)

        # empty list to store one Q-value ADD per action
        Q_adds = []
        
        for a_idx, a_name in enumerate(action_names):
            
            # Regress V through transition
            # Compute the expected discounted future value
            # discount * sum_{s'} P(s'|s,a) * V(s')
            future = regress(mgr, V_primed, subsystems, transition_adds[a_idx], discount)
            
            # Add immediate reward to future value: Q(s,a) = R(s,a) + discount * E[V(s')]
            Q_a = mgr.apply(operator.add, reward_adds[a_name], future)

            # Store the Q-value ADD for this action
            # Q_adds contains one ADD per action, where each ADD maps every state s to the value of taking that action from s
            Q_adds.append(Q_a)
        
        # V_new = max_a Q(s, a)
        # takes the element-wise maximum across all the Q-value ADDs for every state s
        V_new = reduce(lambda a, b: add_max(mgr, a, b), Q_adds)
        
        # Check convergence
        delta = max_diff(mgr, V, V_new)
        if (iteration + 1) % 10 == 0:
            print(f"Iteration {iteration+1}: delta = {delta:.6f}, "
                  f"nodes = {mgr.count_nodes(V_new)}")
        
        V = V_new
        if delta < tolerance:
            print(f"Converged after {iteration+1} iterations.")
            break

    else:
        raise RuntimeError(
            f"SPUDD failed to converge after {max_iter} iterations. "
            f"Final delta = {delta:.6f}. "
            f"Consider increasing max_iter or tolerance."
        )
    
    # Extract policy: for each state, which action gives max Q?
    policy_adds = []
    V_primed = rename_to_primed(mgr, V, subsystems)
    for a_idx, a_name in enumerate(action_names):
        future = regress(mgr, V_primed, subsystems, transition_adds[a_idx], discount)
        Q_a = mgr.apply(operator.add, reward_adds[a_name], future)
        policy_adds.append(Q_a)
    
    return V, policy_adds

    
#######################################
# Utils for interpreting SPUDD output #
#######################################

def get_optimal_action(mgr, reward_adds, policy_adds, state):
    """
    Extract optimal action for a given state
    state: dict e.g. {"X": 50, "Y": 100, "Z": 0}
    Returns index of optimal action.
    """
    
    # get action names
    action_names = list(reward_adds.keys())

    # evaluate each Q-value ADD at a given state
    q_vals = [mgr.evaluate(Q, state) for Q in policy_adds]

    # get index of the highest Q-value
    best = int(np.argmax(q_vals))
    
    print(f"State {state}: Q-values = {dict(zip(action_names, q_vals))}, "
          f"optimal action = {action_names[best]}")
    return best

def get_value_df(mgr, V, states):
    """
    Get optimal value function for each state, return a dataframe.
    """
    rows = []
    
    for x, y, z in itertools.product(states, repeat=3):
        state = {"X": x, "Y": y, "Z": z}
        val = mgr.evaluate(V, state)
        rows.append({"X": x, "Y": y, "Z": z, "V": val})

    df = pd.DataFrame(rows)

    return(df)

def plot_qvalues(mgr, reward_adds, policy_adds, states):
    """
    Plot optimal values for each state.
    """
    action_names = list(reward_adds.keys())

    fig, axes = plt.subplots(len(action_names), len(states), 
                             figsize=(14, 3 * len(action_names)))

    for a_idx, a_name in enumerate(action_names):
        Q = policy_adds[a_idx]
        for col, z in enumerate(states):
            ax = axes[a_idx, col]
            grid = np.array([[mgr.evaluate(Q, {"X": x, "Y": y, "Z": z})
                              for y in states]
                             for x in states])
            im = ax.imshow(grid, origin="lower")
            ax.set_title(f"{a_name}, Z={z}")
            ax.set_xticks(range(len(states)))
            ax.set_xticklabels(states)
            ax.set_yticks(range(len(states)))
            ax.set_yticklabels(states)
            ax.set_xlabel("Y")
            ax.set_ylabel("X")
            plt.colorbar(im, ax=ax)

    plt.suptitle("Q-values per action")
    plt.tight_layout()

    return(plt)