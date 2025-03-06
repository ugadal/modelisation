from hmmlib import *
LO=[SO(seq,1) for defl,seq in getseq(open("allp.fa"))]
nbseq=len(LO)
for obs in LO:obs.pobs/=nbseq
allseq="".join("".join(obs.so) for obs in LO)
alphabet=set(list(allseq))
# ~ um=mkrndmodel(5,alphabet)
pause=True
log=open("perm.log","a")
while True:
	um=mkrndmodel(5,alphabet)
	pp=float("-Inf")
	worse=-pp
	seed=random.random()
	random.seed(seed)
	cyc=0
	worsened=0
	while True:
		cyc+=1
		for so in LO:so.gammab(um)
		tp=sum(so.logcpb for so in LO)
		print(tp)
		if tp<worse:worse=tp
		if tp-pp<0:
			print("******************* WORSENED ****************************")
			worsened+=1
		if abs(tp-pp)<1e-8:
			log.write(f"{seed},{worse},{tp},{cyc},{worsened}\n")
			break
		pp=tp
		um.rep()
		# ~ if pause:
			# ~ input()
			# ~ pause=False
		um.updatetr(LO)
		um.updatepi(LO)
		um.updateem(LO)
	
