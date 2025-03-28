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
	pp=float("-Inf")
	ppb=float("-Inf")
	worse=-pp
	seed=random.random()
	random.seed(seed)
	um=mkrndmodel(5,alphabet)
	cyc=0
	worsened=0
	while True:
		cyc+=1
		for so in LO:so.gammab(um)
		tp=sum(so.logcpb for so in LO)
		# ~ um.temp=1/abs(tp-pp)
		if tp-2*ppb+pp>0:um.nerv-=0.0001
		else:um.nerv+=.0001
		um.temp=exp(um.nerv)
		print(tp,"tampering" if tp-2*ppb+pp>0 else "boosting","delta",tp-ppb,"temp",um.temp,"nerv",um.nerv)
		if tp<worse:worse=tp
		if tp-pp<0:
			print("******************* WORSENED ****************************")
			worsened+=1
		if abs(tp-ppb)<1e-8:
			log=open("perm.log","a")
			log.write(f"{seed},{worse},{tp},{cyc},{worsened}\n")
			log.close()
			break
		pp,ppb=ppb,tp
		um.rep()
		# ~ if pause:
			# ~ input()
			# ~ pause=False
		um.updatetr(LO)
		um.updatepi(LO)
		um.updateem(LO)
	
