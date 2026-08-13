#!/usr/bin/env python3
"""Cumulative-implied hazard offsets (measured-negative experiment; see NOTE.md).
Usage: build_offset.py <april_surface_dump.csv> <april_design.csv> <june_design.csv> <out.csv>"""
import csv, math, bisect, sys

surface, april, june, out = sys.argv[1:5]
G={}; xs=set(); ys=set()
for r in csv.DictReader(open(surface)):
    G[(float(r['wala']),float(r['wac']))]=float(r['g_hat'])
    xs.add(float(r['wala'])); ys.add(float(r['wac']))
xs=sorted(xs); ys=sorted(ys)
tk=tn=0.0
with open(april) as f:
    next(f)
    for line in f:
        p=line.split(','); tn+=float(p[2]); tk+=float(p[3])
off_apr=math.log((tk/tn)/(1-tk/tn))

def pava(v):
    out=[]
    for x in v:
        out.append([x,1.0])
        while len(out)>1 and out[-2][0]>out[-1][0]:
            a=out.pop(); b=out.pop()
            out.append([(a[0]*a[1]+b[0]*b[1])/(a[1]+b[1]), a[1]+b[1]])
    res=[]
    for val,wt in out: res.extend([val]*int(wt))
    return res

Fg={}
for c in ys:
    col=[1/(1+math.exp(-(off_apr+G[(w,c)]/100.0))) for w in xs]
    for w,Fv in zip(xs,pava(col)): Fg[(w,c)]=Fv

def interpF(w,c):
    w=min(max(w,xs[0]),xs[-1]); c=min(max(c,ys[0]),ys[-1])
    i=min(max(bisect.bisect_right(xs,w)-1,0),len(xs)-2)
    j=min(max(bisect.bisect_right(ys,c)-1,0),len(ys)-2)
    x0,x1=xs[i],xs[i+1]; y0,y1=ys[j],ys[j+1]
    tx=(w-x0)/(x1-x0); ty=(c-y0)/(y1-y0)
    return (Fg[(x0,y0)]*(1-tx)*(1-ty)+Fg[(x1,y0)]*tx*(1-ty)
            +Fg[(x0,y1)]*(1-tx)*ty+Fg[(x1,y1)]*tx*ty)

def h_implied(w,c,delta=6.0):
    hi=min(max(w,0.0)+delta/2, xs[-1]); lo=max(hi-delta, 0.0)
    Slo=1-interpF(lo,c); Shi=1-interpF(hi,c)
    ratio=max(min(Shi/max(Slo,1e-9),1.0),1e-9)
    return min(max(1-ratio**(1.0/(hi-lo if hi>lo else 1.0)),2e-4),0.25)

rows=[]; tk=tn=0.0
with open(june) as f:
    next(f)
    for line in f:
        p=line.strip().split(',')
        rows.append(tuple(float(v) for v in p[:4]))
        tn+=float(p[2]); tk+=float(p[3])
off_jun=math.log((tk/tn)/(1-tk/tn))
raw=[math.log(h_implied(w,c)/(1-h_implied(w,c)))-off_jun for w,c,_,_ in rows]
mean=sum(o*r[2] for o,r in zip(raw,rows))/tn
with open(out,'w') as f:
    f.write("wala,wac,n0,k_term,off5\n")
    for (w,c,n,k),o in zip(rows,raw):
        f.write(f"{w},{c},{n:.0f},{k:.0f},{o-mean:.6f}\n")
print("wrote", out)
