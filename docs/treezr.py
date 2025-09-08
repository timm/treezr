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
from typing import Any,List,Iterator
import traceback, random, time, math, sys, re

sys.dont_write_bytecode = True

Number = int|float
Atom   = Number|str|bool
Row    = List[Atom]

big    = 1e32

#--------------------------------------------------------------------
# Stub. Ensure a row is labelled..
def label(row): 
  return row

#--------------------------------------------------------------------
# Create a numeric column summarizer.
def Num(at=0,s=" ") -> o: 
  return o(it=Num, at=at, txt=s, n=0, mu=0, m2=0, sd=0, 
           hi=-big, lo=big, more = 0 if s[-1] == "-" else 1)

# Create a symbolic column summarizer.
def Sym(at=0,s=" ") -> o: 
  return o(it=Sym, at=at, txt=s, n=0, has={})

# Create column summaries from column names.
def Cols(names : List[str]) -> o:
  all=[(Num if s[0].isupper() else Sym)(c,s) for c,s in enumerate(names)]
  klass=None
  for col in all: 
    if col.txt[-1]=="!": klass=col
  return o(it=Cols, names = names, all = all, klass = klass,
           x = [col for col in all if col.txt[-1] not in "X-+"],
           y = [col for col in all if col.txt[-1] in "-+"])

# Create data structure from source rows.
def Data(src) -> o:
  src = iter(src)
  return adds(src, o(it=Data, n=0, mid=None, rows=[], kids=[], ys=None,
                     cols=Cols(next(src)))) # kids=[], how=[], ys=None))

# Create new Data with same columns but different rows.
def clone(data:Data, rows=None) -> o:
  return adds(rows or [], Data([data.cols.names]))

#--------------------------------------------------------------------
# Add multiple items to a summarizer.
def adds(src, it=None) -> o:
  it = it or Num()
  [add(it,x) for x in src]
  return it

# Remove value from summarizer.
def sub(x:o, v:Any, zap=False) -> Any: 
  return add(x,v,-1,zap)

# incrementally update Syms,Nums or Datas.
def add(x: o, v:Any, inc=1, zap=False) -> Any:
  if v == "?": return v
  x.n += inc
  if   x.it is Sym: tmp = x.has[v] = inc + x.has.get(v,0)
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
# Normalize a value to 0..1 range.
def norm(num:Num, v:float) -> float:  
  return  v if v=="?" else (v - num.lo) / (num.hi - num.lo + 1E-32)

# Get central tendencies of all columns.
def mids(data: Data) -> Row:
  data.mid = data.mid or [mid(col) for col in data.cols.all]
  return data.mid

# Get central tendency of one column.
def mid(col: o) -> Atom:
  return max(col.has, key=col.has.get) if col.it is Sym else col.mu

# Return the central tendency for each column..
def divs(data:Data) -> float:
  return [div(col) for col in data.cols.all]

# Return the central tendnacy for one column..
def div(col:o) -> float:
  if col.it is Num: return col.sd
  return -sum(p*math.log(p,2) for n in col.has.values() if (p:=n/col.n) > 0)

#--------------------------------------------------------------------
# Calculate Minkowski distance.
def dist(src) -> float:
  d,n = 0,0
  for v in src: n,d = n+1, d + v**the.p;
  return (d/n) ** (1/the.p)

# Distance from row to best y-values.
def disty(data:Data, row:Row) -> float:
  return dist(abs(norm(c, row[c.at]) - c.more) for c in data.cols.y)

# Sort rows by distance to best y-values.
def distysort(data:Data,rows=None) -> List[Row]:
  return sorted(rows or data.rows, key=lambda r: disty(data,r))

# Distance between two rows using x-values.
def distx(data:Data, row1:Row, row2:Row) -> float:
  def _aha(col, a,b):
    if a==b=="?": return 1
    if col.it is Sym: return a != b
    a,b = norm(col,a), norm(col,b)
    a = a if a != "?" else (0 if b>0.5 else 1)
    b = b if b != "?" else (0 if a>0.5 else 1)
    return abs(a - b)
  return dist(_aha(col, row1[col.at], row2[col.at])  
              for col in data.cols.x)

# Sort rows along a line between 2 distant points..
def distFastmap(data,rows):
  zero, *few = random.choices(rows, k=the.Few)
  D  = lambda r1,r2:distx(data,r1,r2)
  lo = max(few, key= lambda r: D(zero,r))
  hi = max(few, key= lambda r: D(lo,r))
  c  = D(lo,hi)
  x  = lambda row: (D(row,lo)**2 +c*c - D(row,hi)**2)/(2*c + 1e-32)
  return sorted(rows, key=x)

# Assign rows to clusters using recursive random projections..
def distClusters(data, rows=None):
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

# Have we selected this row?.
def treeSelects(row:Row, op:str, at:int, y:Atom) -> bool: 
  return (x := row[at]) == "?" or treeOps[op](x, y)

# Create regression tree..
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

# Divide a col into ranges..
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

#--------------------------------------------------------------
# iterate over all treeNodes.
def treeNodes(data:Data, lvl=0, key=None) -> Data:
  yield lvl, data
  for j in sorted(data.kids, key=key) if key else data.kids:
    yield from treeNodes(j,lvl + 1, key)

# Select a matching leaf.
def treeLeaf(data:Data, row:Row, lvl=0) -> Data:
  for j in data.kids:
    if treeSelects(row, *j.how): return treeLeaf(j,row,lvl+1)
  return data

# Display tree with rows and win columns.
def treeShow(data:Data, key=lambda d: d.ys.mu) -> None:
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
# write the standard error (defaults to no new line).
def fyi(s, end=""):
  print(s, file=sys.stderr, flush=True, end=end)

# coerce a string to int, float, bool, or trimmed string.
def coerce(s:str) -> Atom:
  for fn in [int,float]:
    try: return fn(s)
    except Exception as _: pass
  s = s.strip()
  return {'True':True,'False':False}.get(s,s)

# Returns rows of a csv file..
def csv(file: str ) -> Iterator[Row]:
  with open(file,encoding="utf-8") as f:
    for line in f:
      if (line := line.split("%")[0]):
        yield [coerce(s) for s in line.split(",")]

# shuffle a list, in place.
def shuffle(lst:List) -> List:
  random.shuffle(lst); return lst

# from command line, update config find functions to call.
def _main(settings : o, funs: dict[str,callable]) -> o:
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
  
# The usual run.
def demo():
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

# top-level call.
def main():
  _main(the,globals()); demo()

#---------------------------------------------------------------------
the = o(**{k:coerce(v) for k,v in re.findall(r"(\w+)=(\S+)",__doc__)})
if __name__ == "__main__": main();

