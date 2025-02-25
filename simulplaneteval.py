from hmmlib import *
pa=state("A")
pb=state("B")
pa.E={"N":.5,"S":.3,"P":.2}
pb.E={"N":.3,"S":.5,"P":.2}
pa.T={pa:.8,pb:.2}
pb.T={pa:.4,pb:.6}
Pi={pa:.6,pb:.4}
thismodel=model([pa,pb],Pi,list("NSP"))
fn="obs.txt"
LO=[]
tobs=0
alphabet=[]
for l in open(fn):
	n,so=l.split()
	LO.append(SO(so,int(n)))
	tobs+=int(n)
	alphabet.extend(so)
for so in LO:so.pobs/=tobs
alphabet=set(alphabet)
rmodel=mkrndmodel(8,alphabet)

um=rmodel
# ~ um=thismodel
# ~ for so in LO:so.gammab(um)
# ~ for so in LO[7:8]:
	# ~ print(so.so,so.pobs,so.cpb)
	# ~ so.rep(so.Ab,um)
	# ~ so.rep(so.Bb,um)
	# ~ so.rep(so.Gb,um)
	# ~ um.updatetr(LO)
	# ~ printtrd(so.esttr(um))
	# ~ printemd(so.estem(um))
	# ~ print(so.estpi(um))
# ~ print("done")
# ~ exit()

while True:
	# ~ for so in LO:so.gamma(um)
	for so in LO:so.gammab(um)
	print(sum(so.decal for so in LO))
	um.rep()
	# ~ input()
	# ~ nnn=LO[11]
	# ~ d=nnn.esttr(um)
	# ~ nnn.estem(um)
	# ~ exit()
	# ~ nnn.rep(nnn.A,um)
	# ~ nnn.rep(nnn.B,um)
	# ~ nnn.rep(nnn.G,um)
	# ~ nnn.rep(nnn.Ab,um)
	# ~ nnn.rep(nnn.Bb,um)
	# ~ nnn.rep(nnn.Gb,um)
	# ~ print(nnn.BF.values())
	um.updatetr(LO)
	um.updatepi(LO)
	um.updateem(LO)

	# ~ exit()

