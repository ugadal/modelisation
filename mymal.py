#!/usr/bin/env python3
#
#  mymal.py
#  
from hmm4malib  import *
random.seed(.31415)
LS=[Seq("ABC"),Seq("BCD")]
print(LS)
model=mkrndmodel(LS)
print (",".join(list(model.PI.keys())))
print (",".join(map(str,list(model.PI.values()))))
for s in (I,M):
	print("Em,"+",".join(model.LM[0].M.E.keys()))
	for mod in model.LM:
		print(f"{s}{mod.col},"+",".join(map(str,list(mod.st(s).E.values()))))
for s in (M,D):
	print("tr,"+",".join(model.LM[0].M.T.keys()))
	for mod in model.LM[:-1]:
		print(f"{s}{mod.col},"+",".join(map(str,list(mod.st(s).T.values()))))
for s in (I):
	print("tr,"+",".join(model.LM[0].M.T.keys()))
	for mod in model.LM:
		print(f"{s}{mod.col},"+",".join(map(str,list(mod.st(s).T.values()))))
fs=LS[0]
print(fs)
fs.alphab(model)
