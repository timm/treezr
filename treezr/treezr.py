#!/usr/bin/env python3 
"""
treezr.py: unsupervised tree learning for multi-objective optimization   
(c) 2025, Tim Menzies <timm@ieee.org>, MIT license.    
code: http://github.com/timm/treezr   
data: http://github.com/timm/moot  
Options:

      -A  Any=4             on init, how many initial guesses?   
      -B  Budget=30         when growing theory, how many labels?      
      -C  Check=5           budget for checking learned model   
      -F  Few=64            sample size of data random sampling     
      -l  leaf=3            min items in tree leaves   
      -p  p=2               distance coeffecient   
      -s  seed=1234567891   random number seed      
      -f  file=../moot/optimize/misc/auto93.csv    data file   
      -h                     show help   

"""
from types import SimpleNamespace as o
from typing import Any,Iterator
import traceback, random, time, math, sys, re

sys.dont_write_bytecode = True

Qty     = int | float
Atom    = Qty | str | bool
Row     = list[Atom]

big     = 1e32

#--------------------------------------------------------------------
def label(row:Row) -> Row: 
  "Stub. Ensure a row is labelled."
  return row

#--------------------------------------------------------------------
def Num(at=0,s=" ") -> o: 
  "Create a numeric column summarizer"
  return o(it=Num, at=at, txt=s, n=0, mu=0, m2=0, sd=0, 
           hi=-big, lo=big, more = 0 if s[-1] == "-" else 1)

def Sym(at=0,s=" ") -> o: 
  "Create a symbolic column summarizer"
  return o(it=Sym, at=at, txt=s, n=0, has={})

def Cols(names : list[str]) -> o:
  "Create column summaries from column names"
  all=[(Num if s[0].isupper() else Sym)(c,s) for c,s in enumerate(names)]
  klass=None
  for col in all: 
    if col.txt[-1]=="!": klass=col
  return o(it=Cols, names = names, all = all, klass = klass,
           x = [col for col in all if col.txt[-1] not in "X-+"],
           y = [col for col in all if col.txt[-1] in "-+"])

def Data(src) -> o:
  "Create data structure from source rows"
  src = iter(src)
  return adds(src, o(it=Data, n=0, mid=None, rows=[], kids=[], ys=None,
                     cols=Cols(next(src)))) # kids=[], how=[], ys=None))

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
  return  v if v=="?" else (v - num.lo) / (num.hi - num.lo + 1E-32)

def mids(data: Data) -> Row:
  "Get central tendencies of all columns"
  data.mid = data.mid or [mid(col) for col in data.cols.all]
  return data.mid

def mid(col: o) -> Atom:
  "Get central tendency of one column"
  return max(col.has, key=col.has.get) if col.it is Sym else col.mu

def divs(data:Data) -> float:
  "Return the central tendency for each column."
  return [div(col) for col in data.cols.all]

def div(col:o) -> float:
  "Return the central tendnacy for one column."
  if col.it is Num: return col.sd
  return -sum(p*math.log(p,2) for n in col.has.values() if (p:=n/col.n) > 0)

#--------------------------------------------------------------------
def dist(src) -> float:
  "Calculate Minkowski distance"
  d,n = 0,0
  for v in src: n,d = n+1, d + v**the.p;
  return (d/n) ** (1/the.p)

def disty(data:Data, row:Row) -> float:
  "Distance from row to best y-values"
  return dist(abs(norm(c, row[c.at]) - c.more) for c in data.cols.y)

def distysort(data:Data,rows=None) -> list[Row]:
  "Sort rows by distance to best y-values"
  return sorted(rows or data.rows, key=lambda r: disty(data,r))

def distx(data:Data, row1:Row, row2:Row) -> float:
  "Distance between two rows using x-values"
  def _aha(col, a,b):
    if a==b=="?": return 1
    if col.it is Sym: return a != b
    a,b = norm(col,a), norm(col,b)
    a = a if a != "?" else (0 if b>0.5 else 1)
    b = b if b != "?" else (0 if a>0.5 else 1)
    return abs(a - b)
  return dist(_aha(col, row1[col.at], row2[col.at])  
              for col in data.cols.x)

def distFastmap(data,rows):
  "Sort rows along a line between 2 distant points."
  zero, *few = random.choices(rows, k=the.Few)
  D  = lambda r1,r2:distx(data,r1,r2)
  lo = max(few, key= lambda r: D(zero,r))
  hi = max(few, key= lambda r: D(lo,r))
  c  = D(lo,hi)
  x  = lambda row: (D(row,lo)**2 +c*c - D(row,hi)**2)/(2*c + 1e-32)
  return sorted(rows, key=x)

def distClusters(data, rows=None):
  "Assign rows to clusters using recursive random projections."
  ids  = {}
  rows = shuffle(rows or data.rows)
  stop = 2
  def worker(rows, cid):
    n = len(rows) // 2
    if n >= stop:
      rows = distFastmap(data,rows)
      return worker(rows[:n], 1 + worker(rows[n:], cid))
    else:
      for row in rows: ids[id(row)] = cid
      return cid
  worker(rows,1)
  return ids

#--------------------------------------------------------------------
treeOps = {'<=' : lambda x,y: x <= y, 
           '==' : lambda x,y:x == y, 
           '>'  : lambda x,y:x > y}

def treeSelects(row:Row, op:str, at:int, y:Atom) -> bool: 
  "Have we selected this row?"
  return (x := row[at]) == "?" or treeOps[op](x, y)

def Tree(data:Data, Klass=Num, Y=None, how=None) -> Data:
  "Create regression or decision tree."
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
  "Divide a col into ranges."
  rhs = None, Klass()
  xys = [(r[col.at], add(rhs,Y(r))) for r in rows if r[col.at] != "?"]
  if col.it is Sym:
    d = {}
    [add((d[x] := d.get(x) or Klass()), y) for x,y in xys]
    return o(xpect = sum(c.n/len(xys) * div(c) for c in d.values()),
             hows  = [("==",col.at,x) for x in d])
  out, b4, lhs = None, None, Klass()
  for x,y in sorted(xys):
    if x != b4 and the.leaf <= lhs.n <= len(xys) - the.leaf:
      now = (lhs.n * div(lhs) + rhs.n * div(rhs)) / len(xys)
      if not out or now < out.xpect:
        out = o(xpect=now, hows=[("<=",col.at,b4), (">",col.at,b4)])
    add(lhs, sub(rhs, y))
    b4 = x
  return out

#--------------------------------------------------------------
def treeNodes(data:Data, lvl=0, key=None) -> Data:
  "iterate over all treeNodes"
  yield lvl, data
  for j in sorted(data.kids, key=key) if key else data.kids:
    yield from treeNodes(j,lvl + 1, key)

def treeLeaf(data:Data, row:Row, lvl=0) -> Data:
  "Select a matching leaf"
  for j in data.kids:
    if treeSelects(row, *j.how): return treeLeaf(j,row,lvl+1)
  return data

def treeShow(data:Data, key=lambda d: d.ys.mu) -> None:
  "Display tree with rows and win columns"
  ats = {}
  print(f"{'#rows':>6} {'win':>4}")
  for lvl, d in treeNodes(data, key=key): #\n{100}#
    if lvl == 0: continue
    op, at, y = d.how
    name = data.cols.names[at]
    indent = '|  ' * (lvl - 1)
    expl = f"if {name} {op} {y}"
    score = int(100 * (1 - (d.ys.mu - data.ys.lo) /
             (data.ys.mu - data.ys.lo + 1e-32)))
    leaf = ";" if not d.kids else ""
    print(f"{d.ys.n:6} {score:4}    {indent}{expl}{leaf}")
    ats[at] = 1
  used = [data.cols.names[at] for at in sorted(ats)]
  print(len(data.cols.x), len(used), ', '.join(used))

#--------------------------------------------------------------------
def cdf(x,mu,sd):
  def cdf1(z): return 1 - 0.5*2.718**(-0.717*z - 0.416*z*z)
  z = (x - mu) / sd
  return cdf1(z) if z >= 0 else 1 - cdf1(-z)

#--------------------------------------------------------------------
def fyi(s, end=""):
  "write the standard error (defaults to no new line)"
  print(s, file=sys.stderr, flush=True, end=end)

def coerce(s:str) -> Atom:
  "coerce a string to int, float, bool, or trimmed string"
  for fn in [int,float]:
    try: return fn(s)
    except Exception as _: pass
  s = s.strip()
  return {'True':True,'False':False}.get(s,s)

def csv(file: str ) -> Iterator[Row]:
  "Returns rows of a csv file."
  with open(file,encoding="utf-8") as f:
    for line in f:
      if (line := line.split("%")[0]):
        yield [coerce(s) for s in line.split(",")]

def shuffle(lst:list) -> list:
  "shuffle a list, in place"
  random.shuffle(lst); return lst

def _main(settings : o, funs: dict[str,callable]) -> o:
  "from command line, update config find functions to call"
  for n,s in enumerate(sys.argv):
    if (fn := funs.get(f"eg{s.replace('-', '_')}")):
     try: random.seed(settings.seed); fn()
     except Exception as e:
       print("Error:", e)
       traceback.print_exc()
    else:
      for key in vars(settings):
        if s=="-"+key[0]: 
          settings.__dict__[key] = coerce(sys.argv[n+1])

def dataDwin(file=None):
  data = Data(csv(file or the.file))
  D    = lambda row: disty(data,row)
  b4   = adds(D(row) for row in data.rows)
  return data, D, lambda v: 100*(1 - (v - b4.lo)/(b4.mu - b4.lo))
  
def demo():
  "The usual run"
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

def main():
  "top-level call"
  _main(the,globals()); demo()

#---------------------------------------------------------------------
the = o(**{k:coerce(v) for k,v in re.findall(r"(\w+)=(\S+)",__doc__)})
if __name__ == "__main__": main();
