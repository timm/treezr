#!/usr/bin/env python3 
from types import SimpleNamespace as o
from typing import Iterator, Iterable
import traceback, random, time, sys, re
from math import log

sys.dont_write_bytecode = True

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
Rows = list[Row]
Num  = tuple
Sym  = dict
Col  = Num | Sym
Data = o

the = o(p=2, few=64, seed=1234567891, Few=256, leaf=4,
        file="../../moot/optimize/misc/auto93.csv")

big = 1e32

#--------------------------------------------------------------------
def Data(src: Iterable):
  rows = iter(src)
  data = o(cols=Cols(next(rows)), rows=[])
  [adds(data,row) for row in rows]
  return data

def adds(data,row):
  data.rows += [row]
  for c in data.cols.all: add(data.cols.all, c, row[c]) 
  return row

def add(cols,c,v):
  if v != "?":
    if type(cols[c]) is Sym: 
      cols[c][v] = 1 + cols[c].get(v,0)
    else: 
      cols[c] = (min(v,cols[c][0]), max(v,cols[c][1]))
  return v

def clone(data,rows=[]):
  return Data([data.cols.names] + rows)

def Cols(txt : list[str]) -> o:
  tmp = {c for c,s in enumerate(txt) if s[-1] != "X"}
  all = {c:(big,-big) if txt[c][0].isupper() else Sym() for c in tmp} 
  y   = {c:txt[c][-1]!="-" for c in tmp if txt[c][-1] in "-+"}
  x   = {c for c in tmp if c not in y}
  return o(names=txt, all=all, y=y, x=x)

#--------------------------------------------------------------------
def norm(x,lo,hi): return (x - lo) / (hi - lo + 1e-32) 

def dist(src):
  n,d = 0,0
  for x in src: n,d = n + 1, d + x ** the.p
  return (d/n) ** (1/the.p)

def disty(data:Data, row:Row) -> float:
  return dist(abs(norm(row[c], *data.cols.all[c]) - w)
              for c,w in data.cols.y.items())

def distx(data:Data, row1:Row, row2:Row) -> float:
  def _aha(col, a,b):
    if a==b=="?": return 1
    if type(col) is Sym: return a != b
    a,b = norm(a, *col), norm(b, *col)
    a = a if a != "?" else (0 if b>0.5 else 1)
    b = b if b != "?" else (0 if a>0.5 else 1)
    return abs(a - b)
  return dist(_aha(data.cols.all[x], row1[x], row2[x])  
              for x in data.cols.x)

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
# klass 
def treeSelects(row:Row, op:str, at:int, y:Atom) -> bool: 
  "Have we selected this row?"
  if (x:=row[at]) == "?" : return True
  if op == "<="          : return x <= y
  if op == "=="          : return x == y
  if op == ">"           : return x > y

def Tree(data, rows=None, Y=None,  how=None):
  "Create tree from list of lists"
  rows = rows or data.rows
  Y    = Y or (lambda row: disty(data,row))
  tree = o(rows=rows, how=how, kids=[], mu=mean([Y(r) for r in rows]))
  if len(rows) >= the.leaf:
    spread, cuts = min(treeCuts(x, data.cols.all[x], rows,Y) 
                       for x in data.cols.x)
    if spread < big:
      for cut in cuts:
        subset = [r for r in rows if treeSelects(r, *cut)]
        if the.leaf <= len(subset) < len(rows):
          tree.kids += [Tree(data, subset, Y, cut)]
  return tree

def treeCuts(at, col, rows, Y:callable):
  "Return best cut for column at position 'at'"
  xys = sorted([(r[at], Y(r)) for r in rows if r[at] != "?"])
  return (_symCuts if type(col) is Sym else _numCuts)(at,xys)

def _symCuts(at,xys) -> (float, list):
  "Cuts for symbolic column."
  d = {}
  for x, y in xys:
    d[x]    = d.get(x) or {}
    d[x][y] = 1 + d[x].get(y,0) 
  here = sum(sum(ys.values())/len(xys) * ent(ys) for ys in d.values())
  return here, [("==", at, x) for x in d]

def _numCuts(at,xys) -> (float, list):
  "Cuts for numeric columns."
  spread, cuts, left, right = big, [], {}, {}
  for _, y in xys: 
    right[y] = 1 + right.get(y,0)
  for i, (x, y) in enumerate(xys[:-1]):
    right[y] -= 1
    left[y]   = 1 + left.get(y,0)
    if x != xys[i+1][0]:
      if the.leaf <= i < len(xys) - the.leaf:
        now = (i*ent(left) + (len(xys) - i)*ent(right)) / len(xys)
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
def mean(lst): return sum(lst) / len(lst)

def ent(d):
  N = sum(d.values())
  return -sum(p*log(p,2) for n in d.values() if (p:=n/N) if n>0)

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
  print(len([leaf for _,leaf in treeNodes(tree) if not leaf.kids]))
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
