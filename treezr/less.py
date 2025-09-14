#!/usr/bin/env python3
from typing import Iterator,Iterable
from math import log,sqrt
import random,sys

class o(dict): __getattr__, __setattr__ = dict.get, dict.__setitem__

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
big  = 1e32  # some very large number
the  = o(seed     = 1234567891,
         decimals = 2, 
         file     = "../../moot/optimize/config/SSA-csv")

#--------------------------------------------------------------------
def Num(n=0, mu=0, m2=0, lo=big, hi=-big): return (n,mu,m2,lo,hi)
def Op(op=">=",at=0, value=0)            : return (op,at,value)
def Data(rows, dy, dx, win, names)       : return (rows,dy,dx,win,names)

#--------------------------------------------------------------------
def selects(row, op, at, y):
  x = row[at]
  if x  == "?" : return True
  if op == "<=": return x <= y
  if op == "==": return x == y
  if op == ">" : return x > y
 
def stdev(n,_mu,m2,*_)    : return 0 if n < 2 else (m2 / (n - 1))**0.5

def adds(src, it=None):
  it = it or Num()
  for x in src: it = add(x, *it)
  return it

def sub(x, *num): return add(x, *num, inc=-1)

def add(x, n, mu, m2, lo, hi, inc=1):
  if inc < 0 and n < 2: return Num()
  n  += inc
  d   = x - mu
  mu += inc * (d / n)
  m2 += inc * (d * (x - mu))
  lo  = min(lo,x)
  hi  = max(hi,x)
  return (n, mu, m2, lo, hi)

def eg__stats():
  lst = [random.gauss(10,1) for _ in range(1000)]
  stats = adds(lst)
  for i,x in enumerate(lst): 
    stats = add(x,*stats)
    if i==500: sd = stdev(*stats)
  for i,x in enumerate(lst[::-1]):
    stats = sub(x,*stats)
    if i==500: print(round(sd,4), round(stdev(*stats),4))

def coerce(s):
  try: return int(s)
  except: 
    try: return float(s)
    except: return s.strip()

def Data(src):
  def add(row):
    rows += [row]
    for c in lo:
      if (v := row[c]) != "?": 
        lo[c],hi[c] = min(v,lo[c]), max(v,hi[c])
  
  def _dist(src):
    d,n = 0,0
    for x in src: n += 1; d += x ** the.p
    return (d/n) ** (1/the.p)

  def norm(c, x): 
    return x if x == "?" else (x - lo[c]) / (hi[c] - lo[c] + 1e-32)
  
  def fx1(c, a, b):
    if a == b == "?": return 1
    if c not in w: return a != b
    a, b = norm(c, a), norm(c, b)
    a = a if a != "?" else (0 if b > 0.5 else 1)
    b = b if b != "?" else (0 if a > 0.5 else 1)
    return abs(a - b)
  
  def winner():
    _,__,___,lo,hi = adds(map(data.dy, data.rows))
    return lambda z: int(100*(x - lo) / (hi - lo + 1e-32))

  names, *rows = list(src)
  lo, hi, w = {}, {},{}
  x,y = [],[]
  for c, s in enumerate(names):
    if s[-1] == "X": continue
    (y if s[-1] in "-+" else x).append(c)
    if s.isupper():
      lo[c], hi[c], w[c] = big, -big, 0 if s[-1] == "-" else 1
  [add(row) for row in rows]
  data= o(rows=rows, names=names, add=add,
    cols=o(lo=lo, hi=hi, w=w, x=x, y=y), winner=winnner,
    fx=lambda a,b: _dist(fx1(c, a[c], b[c]) for c in x),
    fy=lambda row: _dist(abs(norm(c, row[c]) - w[c]) for c in y))
  return data

if __name__ == "__main__":
  for n,s in enumerate(sys.argv):
    if (fn := globals().get(f"eg{s.replace('-', '_')}")):
      random.seed(the.seed); fn()
    else:
      for k in the:
        if s=="-"+k[0]: the[k] = coerce(sys.argv[n+1])

exit()

# Tree operations
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

def entropy(d):
  N = sum(d.values())
  return -sum(p*log(p,2) for n in d.values() if (p:=n/N) > 0)

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

def csv(file: str) -> Iterator[Row]:
  with open(file, encoding="utf-8") as f:
    for line in f:
      if (line := line.split("%")[0]): # skip comments
        yield [coerce(s) for s in line.split(",")]

if __name__ == "__main__":
  for n,s in enumerate(sys.argv):
    if (fn := globals().get(f"eg{s.replace('-', '_')}")):
      random.seed(the.seed); fn()
    else:
      for k in vars(the):
        if s=="-"+k[0]: the[k] = coerce(sys.argv[n+1])
      print(the)

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
