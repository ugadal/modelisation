#!/usr/bin/env python3
from hmmlibdp import *
pa=state("A")
pb=state("B")
pa.E={"N":.5,"S":.3,"P":.2}
pb.E={"N":.3,"S":.5,"P":.2}
pa.T={pa:.8,pb:.2}
pb.T={pa:.4,pb:.6}
Pi={pa:.6,pb:.4}
thismodel=model([pa,pb],Pi,list("NSP"))
sss=SO("SSS",1)
sss.gammab(thismodel)
sss.rep(sss.Ab,thismodel)
sss.rep(sss.Bb,thismodel)
sss.rep(sss.Gb,thismodel)
print(sss.cpb)
print(sss.logcpb)
