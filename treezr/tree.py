#!/usr/bin/env python3 
from types import SimpleNamespace as o
from typing import Any,List,Iterator
import traceback, random, time, math, sys, re
sys.dont_write_bytecode = True

Number = int|float
Atom   = Number|str|bool
Row    = List[Atom]

the = o(keep=64, file="../../moot/optimize/misc/auto93.csv")

#--------------------------------------------------------------------
def Cols(names : List[str]) -> o:
  all = [[] for _ in enumerates(names)]
  return o(it=Cols, names = names, all = all, 
           x = [col for s,col in zip(names,all) if s[-1] not in "X-+"],
           y = [col for s,col in zip(names,all) if s[-1]     in "-+" ])

def Data(src):
  data = o(it=Data, rows=[], cols={})
  for m,row in enumerate(src):
    if n==0: 
      data.cols=Cols(row)
    else: 
      [col.append(x) for col,x in zip(data.cols.all, row)]
      data.rows += [rpw]
  for col in data.cols.all: col.sort(key=value)
  return data

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

def value(x):
  return -1e32 if x=="?" else x


