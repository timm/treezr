#!/usr/bin/env python3
from typing import Iterator,Iterable
from math import log,sqrt
import random,sys

class o(dict): __getattr__, __setattr__ = dict.get, dict.__setitem__

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
the  = o(seed     = 1234567891,
         decimals = 2, 
         bins     = 9,
         file     = "../../moot/optimize/config/SSA-csv")

big  = 1e32  # some very large number

#--------------------------------------------------------------------
def is_sym(col): return len(col)==2

def Num(n=0, mu=0, m2=0): return (n, mu, m2)
def Sym(n=0, has=None): return (n, has or {})

def adds(src, col=None):
  col = col or Num()
  for x in src: col = add(x, col)
  return col

def add(x, col):
  def _num(n, mu, m2):
    N = n + 1
    d = x - mu
    MU = mu + d / N
    return (N, MU, m2 + d * (x - MU))
  def _sym(n, has):
    N = n + 1
    has[x] = 1 + has.get(x, 0)
    return (N, has)
  return (_sym if is_sym(col) else _num)(*col)

def cols(rxs):
  if is_sym(rxs[0]):
    N, HAS = 0, {}
    for n, has in rxs:
      N += n
      for s, c in has.items(): 
        HAS[s] = HAS.get(s, 0) + c
    return (N, HAS)
  else:
    N = sum(n for n, mu, m2 in rxs)
    MU = sum(n * mu for n, mu, m2 in rxs) / N
    return (N, MU, sum(m2 + n * (mu - MU)**2 for n, mu, m2 in rxs))

def cdf(x, n,mu,m2):
  fn = lambda z: 1 - 0.5 * exp(-(0.717*z + 0.416*z*z))
  z  = (x - mu) / (m2 / (n - 1))**0.5
  return fn(z) if z > 0 else 1 - fn(-z)

def bin(x, col):
  return x if is_sym(col) else max(bins-1, int(the.bins * cdf(x,*col)))

def div(col):
  def _num(n, _, m2): return 0 if n < 2 else (m2 / (n - 1))**0.5
  def _sym(n, has):
    N = sum(has.values())
    return -sum(p*log(p) for n in has.values() if (p:=n/N) > 0)
  return (_sym if is_sym(col) else _num)(*col)

def bins(col,rows,Y):
  xys = [(x,Y(row)) for row in rows if (x:=row[col.at]) != "?"]
  if is_sym(col):
    d = {}
    for x,y in xys: d[x] = add(y, d.get(x) or Sym())
    return (sum(sym[0]/N * div(sym) for sym in d.values()),
            [("==",col.at,x) for x in d])
  else:
    xys.sort()
    b4, lhs, out = None, Klass(), (big, [])
    for x,y in xys[:-1]:
      if x != b4 and the.leaf <= lhs.n <= len(xys) - the.leaf:
        now = (lhs.n * div(lhs) + rhs.n * div(rhs)) / len(xys)
        if not out or now < out[0]:
          out = (now, [("<=",col.at,b4), (">",col.at,b4)])
      add(lhs, sub(rhs, y))
    b4 = x
  return out
       
a sum
#--------------------------------------------------------------------
def Op(op=">=",at=0, value=0): return (op,at,value)

def selects(row, op, at, y):
  if (x:=row[at]) == "?" : return True
  if op == "<="          : return x <= y
  if op == "=="          : return x == y
  if op == ">"           : return x > y
 def Num(n=0, mu=0, m2=0, lo=big, hi=-big): return (n,mu,m2,lo,hi)

#--------------------------------------------------------------------
def Cols(names):
  lo, hi, w, x, y = {}, {}, {}, [], []
  for c, s in enumerate(names):
    if s[-1] == "X": continue
    (y if s[-1] in "-+" else x).append(c)
    if s.isupper():
      lo[c], hi[c], w[c] = big, -big, 0 if s[-1] == "-" else 1
  cols = o(it=Cols, names=names, lo=lo, hi=hi, w=w, x=x, y=y)
  return cols

def colsClone(i): return Cols(i.names)

def colsAdd(i,row):
  for c in i.lo:
    if (v := row[c]) != "?":
      i.lo[c] = min(v, i.lo[c])
      i.hi[c] = max(v, i.hi[c])

def colsNorm(i, c, x):
  return x if x == "?" else (x - i.lo[c]) / (i.hi[c] - i.lo[c] + 1e-32)

def colsWin(i,rows):
  _,__,___,lo,hi = adds(colsy(i,row) for row in rows)
  return lambda v: int(100 * (v-lo) / (hi - lo + 1e-32))

def colsx(i, r1, r2):
  def fn(c, a, b):
    if a == b == "?": return 1
    if c not in i.w: return a != b
    a, b = colsNorm(i, c, a), colsNorm(i, c, b)
    a = a if a != "?" else (0 if b > 0.5 else 1)
    b = b if b != "?" else (0 if a > 0.5 else 1)
    return abs(a - b)

  return dist(fn(c, r1[c], r2[c]) for c in i.x)

def colsy(i, row):
  return dist(abs(colsNorm(i,c,row[c]) - i.w[c]) for c in i.y)

#--------------------------------------------------------------------
def csv(file: str) -> Iterator[Row]:
  with open(file, encoding="utf-8") as f:
    for line in f:
      if (line := line.split("%")[0]): # skip comments
        yield [coerce(s) for s in line.split(",")]

def entropy(d):
  N = sum(d.values())
  return -sum(p*log(p,2) for n in d.values() if (p:=n/N) > 0)

def coerce(s):
  try: return int(s)
  except: 
    try: return float(s)
    except: return s.strip()

def dist(src, p=2):
  d, n = 0, 0
  for v in src: n, d = n + 1, d + v**p
  return (d/n) ** (1/p)

#--------------------------------------------------------------------
def eg__stats():
  lst = [random.gauss(10,1) for _ in range(1000)]
  stats = adds(lst)
  for i,x in enumerate(lst): 
    stats = add(x,*stats)
    if i==500: sd = stdev(*stats)
  for i,x in enumerate(lst[::-1]):
    stats = sub(x,*stats)
    if i==500: print(round(sd,4), round(stdev(*stats),4))

#--------------------------------------------------------------------
if __name__ == "__main__":
  for n,s in enumerate(sys.argv):
    if (fn := globals().get(f"eg{s.replace('-', '_')}")):
      random.seed(the.seed); fn()
    else:
      for k in the:
        if s=="-"+k[0]: the[k] = coerce(sys.argv[n+1])

exit()

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
