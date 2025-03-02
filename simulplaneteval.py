from hmmlib import *
LO=[SO(seq,1) for defl,seq in getseq(open("allp.fa"))]
nbseq=len(LO)
for obs in LO:obs.pobs/=nbseq
allseq="".join("".join(obs.so) for obs in LO)
alphabet=set(list(allseq))
um=mkrndmodel(5,alphabet)
pause=True
while True:
	for so in LO:so.gammab(um)
	print(sum(so.logcpb for so in LO))
	um.rep()
	if pause:
		input()
		pause=False
	um.updatetr(LO)
	um.updatepi(LO)
	um.updateem(LO)
