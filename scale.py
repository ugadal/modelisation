#!/usr/bin/env python3
#
from hmmlib import *
LO=[SO(seq,1) for defl,seq in getseq(open("allp.fa"))]
nbseq=len(LO)
for obs in LO:obs.pobs/=nbseq
allseq="".join("".join(obs.so) for obs in LO)
alphabet=set(list(allseq))
um=mkrndmodel(5,alphabet)
