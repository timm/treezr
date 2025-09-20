#!/usr/bin/env python3 
from types import SimpleNamespace as o
from typing import Any, Iterator, Iterable
import traceback, random, time, sys, re
from math import log

sys.dont_write_bytecode = True

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
Rows = list[Row]
Num  = list
Sym  = dict
Col  = Num | Sym
Data = o

the = o(p=2, few=64, seed=1234567891, Few=256, leaf=4,
        file="../../moot/optimize/misc/auto93.csv")

big = 1e32

#--------------------------------------------------------------------
def Num(at=0,s=" "): 
  "Create a numeric column summarizer"
  return o(it=Num, at=at, txt=s, n=0, mu=0, m2=0, sd=0, 
           hi=-big, lo=big, more = 0 if s[-1] == "-" else 1)

def Sym(at=0,s=" "): 
  "Create a symbolic column summarizer"
  return o(it=Sym, at=at, txt=s, n=0, has={})

def Cols(names : list[str]) -> o:
  "Create column summaries from column names"
  all=[(Num if s[0].isupper() else Sym)(c,s) 
        for c,s in enumerate(names)]
  return o(it=Cols, names = names, all = all, 
           x = [col for col in all if col.txt[-1] not in "X-+"],
           y = [col for col in all if col.txt[-1] in "-+"])

def Data(src) -> o:
  "Create data structure from source rows"
  src = iter(src)
  return adds(src, o(it=Data, n=0, mid=None, rows=[], kids=[], 
                     ys=None, cols=Cols(next(src)))) 

def clone(data:Data, rows=None) -> o:
  "Create new Data with same columns but different rows"
  return adds(rows or [], Data([data.cols.names]))

#--------------------------------------------------------------------
def adds(src, it=None) -> o:
  "Add multiple items to a summarizer"
  it = it or Num()
  [add(it,x) for x in src]
  return it

def sub(x:o, v:Any, zap=False) -> Any: 
  "Remove value from summarizer"
  return add(x,v,-1,zap)

def add(x: o, v:Any, inc=1, zap=False) -> Any:
  "incrementally update Syms,Nums or Datas"
  if v == "?": return v
  x.n += inc
  if   x.it is Sym: x.has[v] = inc + x.has.get(v,0)
  elif x.it is Num:
    x.lo, x.hi = min(v, x.lo), max(v, x.hi)
    if inc < 0 and x.n < 2:
      x.sd = x.m2 = x.mu = x.n = 0
    else:
      d     = v - x.mu
      x.mu += inc * (d / x.n)
      x.m2 += inc * (d * (v - x.mu))
      x.sd  = 0 if x.n < 2 else (max(0,x.m2)/(x.n-1))**.5
  elif x.it is Data:
    x.mid = None
    if inc > 0: x.rows += [v]
    elif zap: x.rows.remove(v) # slow for long rows
    [add(col, v[col.at], inc) for col in x.cols.all]
  else: raise TypeError(f"cannot add to {type(x)}")
  return v

#--------------------------------------------------------------------
def norm(num:Num, v:float) -> float:  
  "Normalize a value to 0..1 range"
  return v if v=="?" else (v - num.lo) / (num.hi - num.lo + 1E-32)

def mid(col: o) -> Atom:
  "Get central tendency of one column"
  return max(col.has, key=col.has.get) if col.it is Sym else col.mu

def div(col:o) -> float:
  "Return the central tendnacy for one column."
  if col.it is Num: return col.sd
  vs = col.has.values()
  N  = sum(vs)
  return -sum(p*log(p,2) for n in vs if (p:=n/N) > 0)

#--------------------------------------------------------------------
def dist(src) -> float:
  "Calculate Minkowski distance"
  d,n = 0,0
  for v in src: n,d = n+1, d + v**the.p;
  return (d/n) ** (1/the.p)

def disty(data:Data, row:Row) -> float:
  "Distance from row to best y-values"
  return dist(abs(norm(c, row[c.at]) - c.more) for c in data.cols.y)

def distx(data:Data, row1:Row, row2:Row) -> float:
  "Distance between two rows using x-values"
  return dist(_aha(c, row1[c.at], row2[c.at])  for c in data.cols.x)

def _aha(col, a,b):
  "David Aha's distance function."
  if a==b=="?": return 1
  if col.it is Sym: return a != b
  a,b = norm(col,a), norm(col,b)
  a = a if a != "?" else (0 if b>0.5 else 1)
  b = b if b != "?" else (0 if a>0.5 else 1)
  return abs(a - b)

#--------------------------------------------------------------------
def project(data,rows):
  zero, *few = random.choices(rows, k=the.Few)
  D  = lambda r1,r2:distx(data,r1,r2)
  lo = max(few, key= lambda r: D(zero,r))
  hi = max(few, key= lambda r: D(lo,r))
  c  = D(lo,hi)
  x  = lambda row: (D(row,lo)**2 +c*c - D(row,hi)**2)/(2*c + 1e-32)
  return sorted(rows, key=x)

def cluster(data, rows=None, stop=None):
  def go(rows, cid,stop):
    if len(rows) >= stop:
      rows = project(data,rows)
      n = len(rows)//2
      return go(rows[n:], 1 + go(rows[:n], cid,stop), stop)
    for row in rows: ids[id(row)] = cid
    return cid
  ids = {}
  rows = rows or data.rows
  go(shuffle(rows),1,stop or len(rows)**.5)
  return ids

#-------------------------------------------------------------
def treeSelects(row:Row, op:str, at:int, y:Atom) -> bool: 
  "Have we selected this row?"
  if (x:=row[at]) == "?" : return True
  if op == "<="          : return x <= y
  if op == "=="          : return x == y
  if op == ">"           : return x > y

def Tree(data, rows=None, Y=None, Klass=Num, how=None):
  "Create tree from list of lists"
  rows = rows or data.rows
  Y    = Y or (lambda row: disty(data,row))
  tree = o(rows=rows, how=how, kids=[], 
           mu=mid(adds(Y(r) for r in rows)))
  if len(rows) >= the.leaf:
    spread, cuts = min(treeCuts(c,rows,Y,Klass) for c in data.cols.x)
    if spread < big:
      for cut in cuts:
        subset = [r for r in rows if treeSelects(r, *cut)]
        if the.leaf <= len(subset) < len(rows):
          tree.kids += [Tree(data, subset, Y, Klass, cut)]
  return tree

def treeCuts(col, rows, Y:callable, Klass:callable):
  "Return best cut for column at position 'at'"
  xys = sorted([(r[col.at], Y(r)) for r in rows if r[col.at] != "?"])
  return (_symCuts if col.it is Sym else _numCuts)(col.at,xys,Y,Klass)

def _symCuts(at,xys,Y,Klass) -> (float, list):
  "Cuts for symbolic column."
  d = {}
  for x, y in xys:
    d[x] = d.get(x) or Klass()
    add(d[x], y)
  here = sum(ys.n/len(xys) * div(ys) for ys in d.values())
  return here, [("==", at, x) for x in d]

def _numCuts(at,xys,Y,Klass) -> (float, list):
  "Cuts for numeric columns."
  spread, cuts, left, right = big, [], Klass(), Klass()
  [add(right,y) for _, y in xys]
  for i, (x, y) in enumerate(xys[:-1]):
    add(left, sub(right, y))
    if x != xys[i+1][0]:
      if the.leaf <= i < len(xys) - the.leaf:
        now = (left.n*div(left) + right.n*div(right)) / (left.n+right.n)
        if now < spread:
          spread = now
          cuts = [("<=", at, x), (">", at, x)]
  return spread, cuts

# ## Tree Processing -------------------------------------------------
def treeLeaf(tree, row):
  "Find which leaf a row belongs to"
  for kid in tree.kids:
    if treeSelects(row, *kid.how): return treeLeaf(kid, row)
  return tree

def treeNodes(tree, lvl=0):
  "Iterate over all tree nodes"
  yield lvl, tree
  for kid in sorted(tree.kids, key=lambda kid: kid.mu):
    yield from treeNodes(kid, lvl + 1)

def treeShow(data,tree,win=None):
  "Display tree structure with Y means"
  win = win or (lambda v:int(100*v))
  n   = {s:0 for s in data.cols.names}
  print(" ")
  for lvl, node in treeNodes(tree):
    if lvl == 0: continue
    op, at, y = node.how
    indent = '|  ' * (lvl - 1)
    rule = f"if {data.cols.names[at]} {op} {y}"
    n[data.cols.names[at]] += 1
    leaf = ";" if not node.kids else ""
    print(f"n:{len(node.rows):4}   win:{win(node.mu):5}     ",end="")
    print(f"{indent}{rule}{leaf}")
  print("\nUsed: ",*sorted([k for k in n.keys() if n[k]>0],
                           key=lambda k: -n[k]))

#--------------------------------------------------------------------
def coerce(s:str) -> Atom:
  for fn in [int,float]:
    try: return fn(s)
    except Exception as _: pass
  s = s.strip()
  return {'True':True,'False':False}.get(s,s)

def csv(file: str ) -> Iterator[Row]:
  with open(file,encoding="utf-8") as f:
    for line in f:
      if (line := line.split("%")[0]):
        yield [coerce(s) for s in line.split(",")]

def shuffle(lst:list) -> list:
  random.shuffle(lst); return lst

def dataDwin(file=None):
  data = Data(csv(file or the.file))
  D    = lambda row: disty(data,row)
  b4   = adds(D(row) for row in data.rows)
  return data, D, lambda v: 100*(1 - (v - b4.lo)/(b4.mu - b4.lo))
 
def eg__data():
  print(Data(csv(the.file)).cols.all)
 
def eg__cluster():
  print(cluster(Data(csv(the.file))).keys())
 
def eg__tree():
  data  = Data(csv(the.file))
  data1 = clone(data, random.choices(data.rows, k=32))
  tree  = Tree(data1)
  treeShow(data1, tree)

def eg__demo():
  random.seed(the.seed)
  data,D,win = dataDwin(the.file)
  for b in range(10,200,20):
    alls=[]
    somes=[]
    for _ in range(50):
      somes += [int(win(D(min(random.choices(data.rows,k=b//2),key=D))))]
      rows = random.choices(data.rows,k=b)
      ids  = distClusters(data,rows)
      tree = Tree(clone(data, rows),
                Y=lambda row: ids[id(row)],
                Klass=Sym)
      mid=lambda a: a[len(a)//2]
      alls += [max([int(win(D(mid(x.rows)))) for n,x in treeNodes(tree)])]
    print(b//2, mid(sorted(somes)), mid(sorted(alls)))

#---------------------------------------------------------------------
def main(settings : o, funs: dict[str,callable]) -> o:
  for n,s in enumerate(sys.argv):
    if (fn := funs.get(f"eg{s.replace('-', '_')}")):
      random.seed(settings.seed); fn()
    else:
      for key in vars(settings):
        if s=="-"+key[0]: 
          settings.__dict__[key] = coerce(sys.argv[n+1])

if __name__ == "__main__": 
  main(the, globals())
