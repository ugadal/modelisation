#!/usr/bin/env python3
from hmmlibdp import *
fn="obs.txt"
ttobs=0
LO=[]
for line in open(fn):
	a,b=line.split()
	a=int(a)
	LO.append(SO(b,a))
	ttobs+=a
# ~ po=LO[0]
# ~ print(po.so,po.pobs,ttobs)
for obs in LO:
	obs.pobs/=ttobs
# ~ print(po.so,po.pobs,ttobs)
um=mkrndmodel(2,"NSP")
for obs in LO:
	obs.gammab(um)
	print(obs.so,obs.logcpb)
