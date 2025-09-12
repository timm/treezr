def is_sym(x): return len(x)==2

def add(i,v):
  if v != "?":
    match i
      case (n, mu, m2): 
        n+=1; d=v-mu; mu+=d/n; return (n,mu,m2+d*(v-mu))
      case (n,has):
        n+=1; has[v] = 1 + has.get(v,0); return (n,has,most,most)
  return v

def div(it)
  match it:
    case (n,_,m2):  return 0 if n < 2 else (max(0,m2)/(n-1))**.5
    case (n,has):   return -sum(p*math.log(p) for m in has.values() if (p:=m/n)>0)


