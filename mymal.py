#!/usr/bin/env python3
#
#  mymal.py
#  
from hmm4malib  import *
random.seed(.31415)
LS=[Seq("ABC"),Seq("BCD")]
# ~ print(LS)
model=mkrndmodel(LS)
print (";".join(list(model.PI.keys())))#Pi R1
print (";".join(map(str,list(model.PI.values()))))
# ~ for s in (I,M,D):
	# ~ print(";".join())
for s in (I,M): #Em ix row 4+0.1.2.3... , mx 12=4+N+2 
	print("Em;"+";".join(model.LM[0].M.E.keys()))
	for mod in model.LM:
		print(f"{s}{mod.col},"+";".join(map(str,list(mod.st(s).E.values()))))
for s in (M,D): #tr Mx row 20=4+2*N+2+2, Dx = 27
	print("tr;"+";".join(model.LM[0].M.T.keys()))
	for mod in model.LM[:-1]:
		print(f"{s}{mod.col},"+";".join(map(str,list(mod.st(s).T.values()))))
for s in (I):
	print("tr;"+",".join(model.LM[0].M.T.keys()))
	for mod in model.LM:
		print(f"{s}{mod.col},"+";".join(map(str,list(mod.st(s).T.values()))))
exit()
fs=LS[0]
print(fs)
fs.alphab(model)
