from dataclasses import dataclass

#####################################
# Create classes for ADD components #
#####################################

# An ADD is a tree made of two kinds of nodes:
## 1. Leaf — a terminal node holding a number (the function's output value)
## 2. Node — an internal node representing a variable, with one child per possible value of that variable
### note that frozen=True makes them immmutable and hashable

@dataclass(frozen=True)
class Leaf:
    value: float

    # print 
    def __repr__(self):
        return f"Leaf({self.value})"

@dataclass(frozen=True)
class Node:
    var:      str # node name
    children: tuple  # node children; one per state value, ordered by state_vals

    def __repr__(self):
        return f"Node({self.var})"

#######################################################
# Create ADD manager to build and operate on diagrams #
#######################################################

class ADD:
    def __init__(self, state_vals: list):
        self.state_vals  = state_vals

        # hash-consing: if you ask for the same leaf or node twice, the same identity in memory is returned
        # this is what makes ADDs compact: the structure is shared, not duplicated
        self._leaf_cache = {}   # value -> Leaf
        self._node_cache = {}   # (var, children_tuple) -> Node

    def leaf(self, value: float) -> Leaf:
        """Add value to leaf cache."""
        if value not in self._leaf_cache:
            self._leaf_cache[value] = Leaf(value)
        return self._leaf_cache[value]

    def node(self, var: str, children: tuple) -> Leaf | Node:
        """Add state variable and children to node cache."""
        # Reduction rule: all children identical -> return child
        if len(set(id(c) for c in children)) == 1:
            return children[0]
        key = (var, children)
        if key not in self._node_cache:
            self._node_cache[key] = Node(var, children)
        return self._node_cache[key]

    def build(self, fn, var_order: list) -> Leaf | Node:
        """Build an ADD from a function f(assignment) -> float."""
        def _build(depth, assignment):
            """
            depth: tracks variables
            assigment: dictionary of values (eg. {'x': 0, 'y': 50, 'z': 100})
            """
            if depth == len(var_order):
                # add leaf based on assignment and fn()
                return self.leaf(fn(assignment))
            var = var_order[depth]
            children = tuple(
                _build(depth + 1, {**assignment, var: v})
                for v in self.state_vals
            )
            return self.node(var, children)
        # begin at depth = 0 and an empty assignment, {}
        return _build(0, {})

    def evaluate(self, root: Leaf | Node, assignment: dict) -> float:
        """Evaluate ADD at a given assignment; involves traversing the tree based on the assignment and returning the float value at the end."""
        # start at the top of the tree
        node = root
        while isinstance(node, Node):
            val = assignment[node.var]
            idx = self.state_vals.index(val)
            node = node.children[idx]
        return node.value

    def apply(self, op, a: Leaf | Node, b: Leaf | Node) -> Leaf | Node:
        """Apply a binary operation to two ADDs."""
        cache = {}
        def _apply(x, y):
            
            key = (id(x), id(y))

            # if this pair has already been combined, return the cache
            if key in cache:
                return cache[key]
            if isinstance(x, Leaf) and isinstance(y, Leaf):
                result = self.leaf(op(x.value, y.value))
            elif isinstance(x, Node) and isinstance(y, Node) and x.var == y.var:
                result = self.node(x.var, tuple(
                    _apply(cx, cy)
                    for cx, cy in zip(x.children, y.children)
                ))
            elif isinstance(x, Node):
                result = self.node(x.var, tuple(
                    _apply(cx, y) for cx in x.children
                ))
            else:
                result = self.node(y.var, tuple(
                    _apply(x, cy) for cy in y.children
                ))
            cache[key] = result
            return result
        return _apply(a, b)

    def count_nodes(self, root: Leaf | Node) -> int:
        """Count the number of nodes in an ADD."""
        seen = set()
        def _walk(n):
            if id(n) in seen:
                return 0
            seen.add(id(n))
            if isinstance(n, Leaf):
                return 1
            return 1 + sum(_walk(c) for c in n.children)
        return _walk(root)

    def count_leaves(self, root: Leaf | Node) -> int:
        """Count the number of leaves in an ADD."""
        seen = set()
        def _walk(n):
            if id(n) in seen:
                return 0
            seen.add(id(n))
            if isinstance(n, Leaf):
                return 1
            return sum(_walk(c) for c in n.children)
        return _walk(root)

    def to_dot(self, root: Leaf | Node, title: str = "ADD") -> str:
        """Visualize ADD"""
        lines = [
            f'digraph "{title}" {{',
            '  node [fontname="Helvetica"];',
            '  edge [fontname="Helvetica"];',
            ]
        seen = set()

        def node_id(n):
            return f"n{id(n)}"

        def _walk(n):
            if id(n) in seen:
                return
            seen.add(id(n))
            nid = node_id(n)
            if isinstance(n, Leaf):
                lines.append(
                    f'  {nid} [shape=square, style=filled, '
                    f'fillcolor=lightblue, label="{n.value}"];'
                )
            else:
                lines.append(f'  {nid} [shape=ellipse, label="{n.var}"];')
                for val, child in zip(self.state_vals, n.children):
                    lines.append(
                        f'  {nid} -> {node_id(child)} [label="{val}"];'
                    )
                    _walk(child)

        _walk(root)
        lines.append("}")
        return "\n".join(lines)