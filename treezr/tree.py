#!/usr/bin/env python3 
from types import SimpleNamespace as o
from typing import Any,List,Iterator
import traceback, random, time, math, sys, re
sys.dont_write_bytecode = True

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
Num  = list
Sym  = dict
Col  = Num | Sym

the = o(p=2, few=64, file="../../moot/optimize/misc/auto93.csv")

#--------------------------------------------------------------------
def Data(src):
  data = o(rows=[], cols={})
  for n,row in enumerate(src):
    if n==0: data.cols = Cols(row)
    else:
      [add(col,row[c]) for c,col in data.cols.all()]
      data.rows += [row]        
  for col in data.cols.all: ok(col)
  return data

def clone(data,rows=[]):
   return Data([data.cols.rows] + rows)

def Cols(names : List[str]) -> o:
  what = lambda s: Num() if s[0].isupper() else Sym()
  all  = {c:what(s) for c,s in enumerate(names)}
  return o(names = names, all = all, 
    x = {c:all[n] for c,s in enumerate(names) if s[-1] not in "X-+"]},
    y = {c:all[n] for c,s in enumerate(names) if s[-1] in "-+"]})

def adds(src, col=None):
  for x in src:
    col = col or (Sym if type(x) is str else Num)()
    add(col,x)
  return col

def add(col,x):
  if x != "?": 
    if type(col) is Num: col.append(x)
    else: col[x] = 1 + col.get(x,0)

def ok(col,x):
  return sorted(col) if type(col) is Num else col

def div(col):
  if type(col) is Num:
    n=len(col)//10; return (col[9*n] - col[n])/2.56
  else:
    n=len(col)
    return -sum(p*math.log(p) for m in col.values() if (p : =m/n) >0)

def norm(col,x):
  lo,hi = col[0],col[-1]
  return x if x=="?" else (x - lo)/(hi - lo + 1e-32)

def dist(data:Data, row1:Row, row2:Row) -> float:
  n,d = 0,0
  for c,col in data.cols.x.items():
    n += 1
    d += _dist(col,row1[c], row2[c])**the.p
  return (d/n) ** (1/the.p)

def _dist(col, a,b):
  if a==b=="?": return 1
  if type(col) is Sym: return a != b
  a,b = norm(col,a), norm(col,b)
  a = a if a != "?" else (0 if b>0.5 else 1)
  b = b if b != "?" else (0 if a>0.5 else 1)
  return abs(a - b)

def projectToLine(data,rows):
  zero, *few = random.choices(rows, k=the.Few)
  D  = lambda r1,r2:dist(data,r1,r2)
  lo = max(few, key= lambda r: D(zero,r))
  hi = max(few, key= lambda r: D(lo,r))
  c  = D(lo,hi)
  x  = lambda row: (D(row,lo)**2 +c*c - D(row,hi)**2)/(2*c + 1e-32)
  return sorted(rows, key=x)

def cut1(rows): return len(rows)//2

def cut2(ids):
  def fn(rows):
    N         = len(rows)
    out,e,eps = None, 1e32, N**.5
    Y         = lambda row: ids[id(row)]
    l,r       = Sym(), adds([Y(row) for row in rows], Sym())
    for n,row in enumerate(rows):
      r[Y(row)] -= 1
      l[Y(row)] += 1
      if n >= eps and N-n >= eps:
        if (now := (n*div(l) + (N-n)*div(r))/N) < e:
          out,e = n,now
    return out
  return fn

def cluster(data, rows=None, stop=4, cut=cut1):
  def go(rows, cid):
    if len(rows) >= stop:
      if n := cut(projectToLine(data,rows))
        return go(rows[n:]], 1 + go(rows[:n], cid))
    for row in rows: ids[id(row)] = cid
    return cid
  ids  = {}
  go(shuffle(rows or data.rows),1)
  return ids

def clusters(data, rows=None, stop=4):
  return cluster(data, rows, stop, 4, 
                 cut2( cluster(data,rows,stop,4,cut1)))

#--------------------------------------------------------------------
treeOps = {'<=' : lambda x,y: x <= y, 
           '==' : lambda x,y:x == y, 
           '>'  : lambda x,y:x > y}

def treeSelects(row:Row, op:str, at:int, y:Atom) -> bool: 
  return (x := row[at]) == "?" or treeOps[op](x, y)

def Tree(data:Data, Klass=Num, Y=None, how=None) -> Data:
  Y = Y or (lambda row: disty(data, row))
  data.kids, data.how = [], how
  data.ys = adds(Y(row) for row in data.rows)
  if len(data.rows) >= the.leaf:
    hows = [how for col in data.cols.x 
            if (how := treeCuts(col,data.rows,Y,Klass))]
    if hows:
      for how1 in min(hows, key=lambda c: c.xpect).hows:
        rows1 = [r for r in data.rows if treeSelects(r, *how1)]
        if the.leaf <= len(rows1) < len(data.rows):
          data.kids += [Tree(clone(data,rows1), Klass, Y, how1)]
  return data

def treeCuts(col:o, rows:list[Row], Y:callable, Klass:callable) -> o:
  def _sym(sym):
    d, n = {}, 0
    for row in rows:
      if (x := row[col.at]) != "?":
        n += 1
        d[x] = d.get(x) or Klass()
        add(d[x], Y(row))
    return o(xpect = sum(c.n/n * div(c) for c in d.values()),
             hows = [("==",col.at,x) for x in d])

  def _num(num):
    out, b4, lhs, rhs = None, None, Klass(), Klass()
    xys = [(row[col.at], add(rhs, Y(row))) # add returns the "y" value
           for row in rows if row[col.at] != "?"]
    for x, y in sorted(xys, key=lambda z: z[0]):
      if x != b4 and the.leaf <= lhs.n <= len(xys) - the.leaf:
        now = (lhs.n * div(lhs) + rhs.n * div(rhs)) / len(xys)
        if not out or now < out.xpect:
          out = o(xpect=now, hows=[("<=",col.at,b4), (">",col.at,b4)])
      add(lhs, sub(rhs, y))
      b4 = x
    return out

  return (_sym if col.it is Sym else _num)(col)

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

def shuffle(lst:List) -> List:
  random.shuffle(lst); return lst

def value(x):
  return -1e32 if x=="?" else x


