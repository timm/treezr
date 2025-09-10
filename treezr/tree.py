#!/usr/bin/env python3 
from types import SimpleNamespace as o
from typing import Any,Iterator
import traceback, random, time, math, sys, re
sys.dont_write_bytecode = True

Qty  = int | float
Atom = Qty | str | bool
Row  = list[Atom]
Num  = list
Sym  = dict
Col  = Num | Sym

the = o(p=2, few=64, file="../../moot/optimize/misc/auto93.csv",
        seed=1234567891)

#--------------------------------------------------------------------
def Data(src):
  data = o(rows=[], cols={})
  for n,row in enumerate(src):
    if n==0: data.cols = Cols(row)
    else:
      [add(col,row[c]) for c,col in data.cols.all.items()]
      data.rows += [row]        
  for col in data.cols.all: ok(col)
  return data

def clone(data,rows=[]):
   return Data([data.cols.rows] + rows)

def Cols(names : list[str]) -> o:
  what = lambda s: Num() if s[0].isupper() else Sym()
  all  = {c:what(s) for c,s in enumerate(names)}
  y    = {c:all[c]  for c,s in enumerate(names) if s[-1] in "-+"},
  w    = {c:(0 if s[-1] =="-" else 1) for c,s in enumerate(names) if c in y}
  return o(names = names, all = all, w=w, y=y,
    x = {c:all[c] for c,s in enumerate(names) if s[-1] not in "X-+"})

#--------------------------------------------------------------------
def adds(src, col=None):
  for x in src:
    col = col or (Sym if type(x) is str else Num)()
    add(col,x)
  return col

def add(col,x):
  if x != "?": 
    if type(col) is Num: col.append(x)
    else: col[x] = 1 + col.get(x,0)

def ok(col):
  return sorted(col) if type(col) is Num else col

def div(col):
  if type(col) is Num:
    n=len(col)//10; return (col[9*n] - col[n])/2.56
  else:
    N = sum(col.values())
    return -sum(p*math.log(p) for n in col.values() if (p:=n/N) >0)

def norm(col,x):
  lo,hi = col[0],col[-1]
  return x if x=="?" else (x - lo)/(hi - lo + 1e-32)

#--------------------------------------------------------------------
def dist(src) -> float:
  d,n = 0,0
  for v in src: n,d = n+1, d + v**the.p;
  return (d/n) ** (1/the.p)

def disty(data:Data, row:Row) -> float:
  return dist(abs(norm(col, row[c]) - w[c]) 
              for c,col in data.cols.y.items())

def besty(data:Data,rows=None) -> list[Row]:
  return min(rows or data.rows, key=lambda r: disty(data,r))

def distx(data:Data, row1:Row, row2:Row) -> float:
  def _aha(col, a,b):
    if a==b=="?": return 1
    if type(col) is Sym: return a != b
    a,b = norm(col,a), norm(col,b)
    a = a if a != "?" else (0 if b>0.5 else 1)
    b = b if b != "?" else (0 if a>0.5 else 1)
    return abs(a - b)
  return dist(_aha(col, row1[c], row2[c])  
              for c,col in data.cols.x.items())

#--------------------------------------------------------------------
def projectToLine(data,rows):
  zero, *few = random.choices(rows, k=the.Few)
  D  = lambda r1,r2:distx(data,r1,r2)
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
      if n := cut(projectToLine(data,rows)):
        return go(rows[n:], 1 + go(rows[:n], cid))
    for row in rows: ids[id(row)] = cid
    return cid
  ids = {}
  go(shuffle(rows or data.rows),1)
  return ids

def clusters(data, rows=None, stop=4):
  ids= cluster(data,rows,stop,4,cut1)
  return ids, cluster(data, rows, stop, 4, cut2(ids))

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

def value(x):
  return -1e32 if x=="?" else x

def dataDwin(file=None):
  data = Data(csv(file or the.file))
  D    = lambda row: disty(data,row)
  b4   = adds(D(row) for row in data.rows)
  return data, D, lambda v: 100*(1 - (v - b4.lo)/(b4.mu - b4.lo))
  
def eg__demo():
  random.seed(the.seed)
  data,D,win = dataDwin(the.file)
  for b in range(10,200,20):
    print(".")
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
