#!/usr/bin/env python3
from types import SimpleNamespace as o
import math

# Config
the_leaf = 3
big = 1e32

from typing import NamedTuple, get_type_hints

def of(fn, doc=""):
  cls = next(iter(get_type_hints(fn).values()), None)  # first type hint
  if cls: setattr(cls, fn.__name__, fn); fn.__doc__ = doc 
  return fn

# Incremental stats helpers
def Stats(n=0, mu=0, m2=0, lo=big, hi=-big): return (n,mu,m2,lo,hi)

def sd(n,_mu,m2,*_): return 0 if n < 2 else (m2 / (n - 1))**0.5

def add(src, stats=None, inc=1):
   n,mu,m2,lo,hi = stats or Stats()
   for x in hasattr(src,"__iter__") or [src]:
      n  += inc
      d   = x - mu
      mu += inc * (d / n)
      m2 += inc * (d * (x - mu))
      lo  = min(lo,x)
      hi  = max(hi,x)
   return (n, mu, m2, lo, hi)

def sub(x,stats): return add(x,stats,-1)
 
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

def entropy(d):
  N = sum(d.values())
  return -sum(p*math.log(p,2) for n in d.values() if (p:=n/N) > 0)

def treeCuts(at, is_sym, rows, Y):
  "Return best cut for column at position 'at'"
  div, cuts = big, []
  xys = [(r[at], Y(r)) for r in rows if r[at] != "?"]
  if is_sym:
    d = {}
    for x, y in xys:
      d[x]    = d.get(x) or {}
      d[x][y] = d[x].get(y,0) + 1
    here = sum(sum(ys.values())/len(xys) * entropy(ys) for ys in d.values())
    div, cuts = here, [("==", at, x) for x in d]
  else:
    xys.sort()
    l,r = Stats(), Stats()
    for _, y in xys: r = add(y,r)
    for i, (x, y) in enumerate(xys[:-1]):
      l = add(y,l)
      r = sub(y,r)
      if x != xys[i+1][0]:
        if the_leaf <= l[0] <= len(xys) - the_leaf:
          here = (i*sd(*l) + (len(xys) - i)*sd(*r)) / len(xys)
          if here < div:
            div, cuts = here, [("<=", at, x), (">", at, x)]
  return div, cuts

def Tree(rows, names, Y=None, how=None):
  "Create tree from list of lists"
  tree = o(rows=rows, how=how, kids=[], mu=add(Y(r) for r in rows)[1])
  if len(rows) >= the_leaf:
    div, cuts = min(treeCuts(at, s[0].islower(), rows, Y)
                    for at,s in enumerate(names) if s[-1] not in "X-+" )
    if div < big:
      for cut in cuts:
        subset = [r for r in rows if selects(r, *cut)]
        if the_leaf <= len(subset) < len(rows):
          tree.kids += [Tree(subset, names, Y, cut)]
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
