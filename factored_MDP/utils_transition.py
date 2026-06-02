######################################################
# build transition ADD for one action, one subsystem #
######################################################

def build_transition_add_single(mgr, get_transition, params, states, a_binary,
                                subsystem):
    """
    Build ADD for one action based on transition matrix, P(var_next | var_current, a).
    subsystem: name of subsystem e.g. "X", "Y", "Z"
    """

    var_current = subsystem
    var_next = subsystem + "_p"

    T = get_transition(params, states, a_binary)

    # T[i, j] = P(s'=states[i] | s=states[j])
    # ADD maps (var_current, var_next) -> probability
    def prob_fn(assignment):
        # get column index
        j = states.index(assignment[var_current])
        # get row index
        i = states.index(assignment[var_next])
        return T[i, j]

    return mgr.build(prob_fn, [var_current, var_next])

##################################################
# build transition for all actions and subsystem #
##################################################

def build_transition_adds(mgr, get_transition, params, states, action_set, n_actions,
                          subsystems, verbose = True):
    """
    Build ADD for all actions based on transition matrices, P(sprime | s, a_i).
    subsystems: list of names of subsystem e.g. ["X", "Y", "Z"]
    """

    transition_adds = {}   # (action_idx, var) -> ADD
    for a_idx in range(n_actions):
        transition_adds[a_idx] = {}
        for sub_idx, subsystem in enumerate(subsystems):
            a_binary = int(action_set[sub_idx, a_idx])
            add = build_transition_add_single(
                mgr, get_transition, params, states, a_binary,
                subsystem
            )
            transition_adds[a_idx][subsystem] = add
            if verbose:
                print(f"Action {a_idx}, {subsystem}: "
                      f"{mgr.count_nodes(add)} nodes, "
                      f"{mgr.count_leaves(add)} leaves")
    return transition_adds