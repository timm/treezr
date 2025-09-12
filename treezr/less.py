#!/usr/bin/env python3
from types import SimpleNamespace as o
import math

# Config
the_leaf = 3
big = 1e32

# Incremental stats helpers
def sd(n, m2): 
    return 0 if n < 2 else (m2 / (n - 1))**0.5

def adds(lst,stats=None,inc=1):
   n,mu,m2=stats or (0,0,0)
   for y in lst:
      n += inc
      d = y - mu
      mu += inc * (d / n)
      stats = (n.mu, m2 + inc * (d * (y - mu)))
   return stats

def subs(lst,stats): return adds(lst,stats,-1)
 
def entropy(lst):
    if not lst: return 0
    counts = {}
    for x in lst: counts[x] = counts.get(x, 0) + 1
    n = len(lst)
    return -sum((c/n) * math.log(c/n, 2) for c in counts.values())

# Tree operations
def selects(row, op, at, y):
  x = row[at]
  if x == "?": return True
  if op == "<=": return x <= y
  if op == "==": return x == y
  if op == ">": return x > y

def ranges(rows, names):
    "Pre-compute lo/hi for all y-columns"
    y_cols = [i for i, name in enumerate(names) if name[-1] in "+-"]
    rngs = {}
    for at in y_cols:
        vals = [r[at] for r in rows if r[at] != "?"]
        if vals: rngs[at] = (min(vals), max(vals))
    return rngs

def disty(row, rngs, names):
    "Distance from row to best y-values (using pre-computed ranges)"
    y_cols = [i for i, name in enumerate(names) if name[-1] in "+-"]
    if not y_cols: return 0
    
    dists = []
    for at in y_cols:
        if at not in rngs: continue
        lo, hi = rngs[at]
        norm_val = (row[at] - lo) / (hi - lo + 1e-32) if row[at] != "?" else 0.5
        ideal = 1.0 if names[at][-1] == "+" else 0.0
        dists.append(abs(norm_val - ideal))
    
    return sum(d*d for d in dists) ** 0.5 if dists else 0

# Core functions (given)
def coerce(s):
    for fn in [int, float]:
        try: return fn(s)
        except: pass
    s = s.strip()
    return {'True': True, 'False': False}.get(s, s)

def csv(file):
    with open(file, encoding="utf-8") as f:
        for line in f:
            if (line := line.split("%")[0]):
                yield [coerce(s) for s in line.split(",")]

def treeCuts(at, is_sym, rows, Y):
    "Return best cut for column at position 'at'"
    valid = [(r[at], Y(r)) for r in rows if r[at] != "?"]
    if not valid: return big, []
    
    if is_sym:
        groups = {}
        for x, y in valid:
            if x not in groups: groups[x] = []
            groups[x].append(y)
        diversity = sum(len(ys)/len(valid) * entropy(ys) 
                       for ys in groups.values())
        return diversity, [("==", at, x) for x in groups]
    
    # Numeric: fast incremental stats with tuples
    valid.sort()
    best_div, best_cuts = big, []
    
    # Initialize RHS with all values
    lhs,rhs = (0,0,0), (0, 0, 0)
    for _, y in valid: rhs = adds([y],rhs)
    
    
    for i, (x, y) in enumerate(valid[:-1]):
        # Move y from RHS to LHS
        lhs, rhs = adds([y],lhs), subs([y],rhs)
        if (x != valid[i+1][0] and 
            the_leaf <= lhs[0] <= len(valid) - the_leaf):
            div = (lhs[0] * sd(lhs[0], lhs[2]) + rhs[0] * sd(rhs[0], rhs[2])) / len(valid)
            if div < best_div:
                best_div = div
                best_cuts = [("<=", at, x), (">", at, x)]
    
    return best_div, best_cuts
   

def Tree(rows, names, Y=None, how=None):
    "Create tree from list of lists"
    Y = Y or (lambda row: 0)
    tree = o(rows=rows, how=how, kids=[], 
             mu=adds(Y(row) for row in rows])[1])
u)
    
    if len(rows) >= the_leaf:
        x_cols = [i for i, name in enumerate(names) 
                 if name[-1] not in "X-+"]
        
        best_div, best_cuts = min(
            treeCuts(at, not names[at][0].isupper(), rows, Y)
            for at in x_cols)
        
        if best_div < big:
            for cut in best_cuts:
                subset = [r for r in rows if selects(r, *cut)]
                if the_leaf <= len(subset) < len(rows):
                    tree.kids.append(Tree(subset, names, Y, cut))
    
    return tree

def treeLeaf(tree, row):
    "Find which leaf a row belongs to"
    for kid in tree.kids:
        if selects(row, *kid.how): 
            return treeLeaf(kid, row)
    return tree

def treeNodes(tree, lvl=0):
    "Iterate over all tree nodes"
    yield lvl, tree
    for kid in tree.kids:
        yield from treeNodes(kid, lvl + 1)

def treeShow(tree, names):
    "Display tree structure with Y means"
    for lvl, node in treeNodes(tree):
        if lvl == 0: continue
        op, at, y = node.how
        indent = '|  ' * (lvl - 1)
        rule = f"if {names[at]} {op} {y}"
        leaf = ";" if not node.kids else ""
        print(f"{len(node.rows):4} {node.mu:6.2f} {indent}{rule}{leaf}")

# Example usage:
if __name__ == "__main__":
    # Read data
    data = list(csv("auto93.csv"))
    names = data[0]  # header row: ['Clndrs','Volume','Hp+','Lbs-',...]
    rows = data[1:]  # data rows
    
    # Pre-compute ranges once for efficiency  
    rngs = ranges(rows, names)
    
    # Build tree optimizing single objective (Mpg+)
    tree1 = Tree(rows, names, Y=lambda row: row[7])
    
    # Build tree optimizing multiple objectives using disty
    tree2 = Tree(rows, names, Y=lambda row: disty(row, rngs, names))
    
    # Show both trees
    print("Single objective (Mpg+):")
    print(f"{'#rows':>4} {'Y_mu':>6}")
    treeShow(tree1, names)
    
    print("\nMulti-objective (all +/- columns):")  
    print(f"{'#rows':>4} {'Y_mu':>6}")
    treeShow(tree2, names)
    
    # Find leaf for a specific row
    leaf = treeLeaf(tree2, rows[0])
    print(f"\nMulti-objective prediction: {leaf.mu:.3f}")
