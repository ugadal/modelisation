import random
from math import log
from functools import reduce
from functools import cache
I="I"
M="M"
D="D"
random.seed(.31415)
def getseq(df):
	line="#"
	while not line.startswith(">"):line=df.readline().strip()
	defline=line[1:]
	line=df.readline().strip()
	seq=[]
	while line:
		if line.startswith(">"):
			yield defline,"".join(seq)
			seq=[]
			defline=line
		else:seq.append(line)
		line=df.readline().strip()
	yield defline,"".join(seq)
"""
time_obs/module		0		1		2		3		4		5		EIC
___________________________________________________________________________
		0			0		1		2		3		4		5		6
		1			1		2		3		4		5		6		7
		2			2		3		4		5		6		7		8
		EDR			3		4		5		6		7		8		

nbmod=6  ; 0->5
nbobs=3  ; 0->2
for layer in range (nbmod+nbos):
	if layer <nbobs: 
		edr=eic=None
	else:	
		edr= module:layer-nbobs, time:nbobs
	if layer <nbmod:
		eic=None
	else:
		eic=  module:nbmod, time: layer-nbmod
	inner=[(layer-t,t) for t in range(nbobs) if layer-t<nbmod]
	return [eic,*inner,edr]
"""
@cache
def getlayers(layer,nbmod,nbobs):
	print(f"layer : {layer} ,{nbmod}, {nbobs}")
	if layer<nbobs:edr=eic=None
	else:edr=(D,layer-nbobs,nbobs)
	if layer<nbmod:eic=None
	else:eic=(I,nbmod,layer-nbmod)
	inner=[(st,layer-t,t) for t in range(nbobs) if 0<=layer-t<nbmod for st in (I,M,D)]
	return[eic,*inner,edr]

def normd(D):
	trsh=1e-5
	if any(v==1.0 for v in D.values()):
		new={k:0 if v!=1.0 else 1.0 for k,v in D.items()}
		return new
	t=sum(D.values())
	new={k:v/t for k,v in D.items()}
	# ~ if any(v==1.0 for v in new.values()):
		# ~ new={k:0 if v!=1.0 else 1.0 for k,v in new.items()}
	while True:
		V=list(new.values())
		V.sort()
		mi=next(v for v in V if v>0)
		mx=max(new.values())
		if mi/mx<trsh:
			new={k:0 if v/mx<trsh  else v for k,v in new.items()}
			# ~ new=normd(new)
			continue
		break
	return new
	# ~ return {k:v/t for k,v in D.items()}
def one(*args):return 1
class state():
	def __init__(self,n="pas de nom"):
		# ~ print("initialisation d'une instance d'etat",n)
		self.name=n
		self.E={}
		self.T={}
	def emit(self,symb):return self.E.get(symb)
	def transit(self,target):return self.T.get(target)
	def repemit(self):
		print(self.name,1/(1-self.T[self]),"".join([s*int(100*v) for s,v in self.E.items()]))
class module():
	def __init__(self,c,alphabet):
		self.col=c
		self.I=state()
		self.I.E=makernddic(alphabet)
		self.M=state()
		self.M.E=makernddic(alphabet)
		self.D=state()
		self.D.emit=one
		self.to_mod=None
		self.from_mod=None
	def st(self,s):
		match s:
			case "I":return self.I
			case "D":return self.D
			case "M":return self.M
		
class model():
	def __init__(self,LM,PI,LSEQ):
		self.LM=LM
		self.PI=PI
		self.LSEQ=LSEQ
		self.nbmod=len(LM)-1
		print("Pi:",PI,len(LM),self.nbmod)
		# ~ self.alphabet=alphabet
	# ~ def updatepi(self,LO):
		# ~ newpi={st:0 for st in self.LE}
		# ~ for so in LO:
			# ~ tt=so.estpi(self)
			# ~ for k,v in tt.items():
				# ~ newpi[k]+=v*so.pobs
		# ~ self.PI=normd(newpi)
	# ~ def updateem(self,LO):
		# ~ newe={st:{observable:0 for observable in self.alphabet} for st in self.LE}
		# ~ SG={st:0 for st in self.LE}
		# ~ for so in LO:
			# ~ tt,tsg=so.estem(self)
			# ~ for st in self.LE:
				# ~ for obs in self.alphabet:
					# ~ newe[st][obs]+=tt[st][obs]*so.pobs
		# ~ for st in self.LE:
			# ~ st.E=normd(newe[st])

	# ~ def updatetr(self,LO):
		# ~ newtr={source:{k:0 for k in self.LE} for source in self.LE}
		# ~ SG={st:0 for st in self.LE}
		# ~ for so in LO:
			# ~ nt,tsg=so.esttr(self)
			# ~ for source in self.LE:
				# ~ for target in self.LE:
					# ~ newtr[source][target]+=nt[source][target]*so.pobs
		# ~ for source in self.LE:
			# ~ source.T=normd(newtr[source])
	def rep(self):
		print([(k.name,v) for k,v in self.PI.items()])
		for s in self.LE:
			s.repemit()
			# ~ print(s.name,[(k,v) for k,v in s.E.items()])
		for s in self.LE:
			print(s.name,[(k.name,v) for k,v in s.T.items()])
				
def makernddic(L):
	P=[random.random() for _ in list(L)]
	t=sum(P)
	P=[p/t for p in P]
	return {o:p for o,p in zip(L,P)}
def mkrndmodel(LSEQ):
	ttl=sum(seq.L for seq in LSEQ)
	alphabet=set()
	for seq in LSEQ:alphabet=alphabet or set(seq.so)
	alphabet=list(alphabet)
	alphabet.sort()
	LM=[module(x,alphabet) for x in range(ttl)]
	for source in LM:
		source.M.T=makernddic(list("IMD"))
		source.D.T=makernddic(list("IMD"))
		source.I.T=makernddic(list("IMD"))
	lastmod=LM[-1]
	lastmod.D.T={I:1}
	lastmod.M.T={I:1}
	final=module(ttl+1,alphabet)
	final.I.T=makernddic(list("IMD"))
	print(final.I.T)
	final.I.T["I"]=1.0
	print(final.I.T)
	final.I.T=normd(final.I.T)
	print(final.I.T)
	LM.append(final)
	for source,dest in zip(LM,LM[1:]):
		source.to_mod=dest
		dest.from_mod=source
	pi=makernddic((list("IMD")))
	return model(LM,pi,LSEQ)
def printtrd(D):
	for k,d in D.items()		:
		print(k.name,[(t.name,v) for t,v in d.items()])
def printemd(D):
	for k,d in D.items()		:
		print(k.name,[(t,v) for t,v in d.items()])
class Seq():
	def __init__(self,so,pobs=1):
		self.so=list(so)
		self.L=len(so)
		self.pobs=pobs
		self.Ab={}
		self.Bb={}
		self.Gb={}
		self.BF={}
	def alphab(self,model):
		# ~ layer 0
		current=model.LM[0]
		self.Ab[I,0,0]=model.PI[I] * current.I.emit(self.so[0])
		self.Ab[M,0,0]=model.PI[M] * current.M.emit(self.so[0])
		self.Ab[D,0,0]=model.PI[D] 
		# ~ print(model.PI[I] , current.I.emit(self.so[0]),self.so[0])
		# ~ print("IOO",self.Ab[I,0,0])
		# ~ print("MOO",self.Ab[M,0,0])
		# ~ print("DOO",self.Ab[D,0,0])
		bf=self.Ab[I,0,0]+self.Ab[M,0,0]+self.Ab[D,0,0]
		self.Ab[I,0,0]/=bf
		self.Ab[M,0,0]/=bf
		self.Ab[D,0,0]/=bf
		self.BF[0]=bf
		# ~ print(bf)
		# ~ layer 1
		current=current.to_mod
		previous=current.from_mod
		bf=0
		v = self.Ab[D,0,0] * previous.D.transit(I) * current.I.emit(self.so[0])
		# ~ print(self.Ab[D,0,0] , previous.D.transit(I) , current.I.emit(self.so[0]))
		bf+=v
		self.Ab[I,1,0] = v
		v = self.Ab[D,0,0] * previous.D.transit(M) * current.M.emit(self.so[0])
		# ~ print(self.Ab[D,0,0] , previous.D.transit(M) , current.M.emit(self.so[0]))
		bf+=v
		self.Ab[M,1,0] = v
		v = self.Ab[D,0,0] * previous.D.transit(D) 
		# ~ print(self.Ab[D,0,0] , previous.D.transit(D))
		bf+=v
		self.Ab[D,1,0] = v
		
		v = self.Ab[I,0,0] * previous.I.transit(I) * previous.I.emit(self.so[1])
		bf+=v
		self.Ab[I,0,1] = v
		v = self.Ab[I,0,0] * previous.I.transit(M) * previous.M.emit(self.so[1])
		bf+=v
		self.Ab[M,0,1] = v
		v = self.Ab[I,0,0] * previous.I.transit(D) 
		bf+=v
		self.Ab[D,0,1] = v
		self.Ab[I,1,0]/=bf
		self.Ab[M,1,0]/=bf
		self.Ab[D,1,0]/=bf
		self.Ab[I,0,1]/=bf
		self.Ab[M,0,1]/=bf
		self.Ab[D,0,1]/=bf
		self.BF[1]=bf
		for k,v in self.Ab.items():print(k,v)
		for layer in range(2,self.L+model.nbmod):
			cells=getlayers(layer,model.nbmod,self.L)
			print(cells)
			bf=0
			if cells[0]: # I run
				state,modnum,t=cells[0]
				current=model.LM[modnum]
				previous=current.from_mod
				v = self.Ab[D,modnum-1,t] * previous.D.transit(I) * current.I.emit(self.so[t])
				if t>0:
					v+=self.Ab[I,modnum,t-1] * current.I.transit(I) * current.I.emit(self.so[t])
					v+=self.Ab[M,modnum-1,t-1] * previous.M.transit(I) * current.I.emit(self.so[t])
				self.Ab[state,modnum,t]=v
				bf+=v
			if cells[-1]: # D run
				state,modnum,t=cells[-1]
				current=model.LM[modnum]
				v = self.Ab[I,modnum,t-1] * current.I.transit(D)
				if modnum>0:
					previous=current.from_mod
					v+=self.Ab[D,modnum-1,t] * previous.D.transit(D)
					v+=self.Ab[M,modnum-1,t-1] * previous.M.transit(D) 
				self.Ab[state,modnum,t]=v
				bf+=v
			for state,modnum,t in cells[1:-1]:
				print(f"fiddling on {state}, {modnum},{t}")
				current=model.LM[modnum]
				if t==0:
					v=self.Ab[D,modnum-1,t] * current.from_mod.D.transit(state) * current.st(state).emit(self.so[t])
				elif modnum==0:
					v=self.Ab[I,modnum,t-1] * current.I.transit(state) * current.st(state).emit(self.so[t])
				else:
					v=self.Ab[D,modnum-1,t] * current.from_mod.D.transit(state) * current.st(state).emit(self.so[t])
					v+=self.Ab[I,modnum,t-1] * current.I.transit(state) * current.st(state).emit(self.so[t])
					v+=self.Ab[M,modnum-1,t-1] * current.from_mod.M.transit(state) * current.st(state).emit(self.so[t])*self.BF[layer-1]
				bf+=v
				self.Ab[state,modnum,t]=v
			for state,modnum,t in [c for c in cells if c]:
				self.Ab[state,modnum,t]/=bf
				print(state,modnum,t,self.Ab[state,modnum,t])
			self.BF[layer]=bf
			print(self.BF)
		# ~ for layer in range(2:self.L+model.nbmod):
			# ~ cells=getlayers(layer,model.nbmod,self.L)
			# ~ if cells[0]: #EIC True
				# ~ mod,t=cells[0]
				# ~ module=model.LM[module]
				# ~ previous=module.from_mod
				# ~ self.Ab[I,mod,t]=self.Ab[D,mod-1,t] * previous.D.transits(I) *
		# ~ layer 0
		# ~ model.
		# ~ T=[model.PI[s] for s in model.LE]
		# ~ for t,symb in enumerate(self.so):
			# ~ T=[t*s.emit(symb) for t,s in zip(T,model.LE)] # ~ emission
			# ~ bf=sum(T)
			# ~ T= [t/bf for t in T]
			# ~ self.BF[t]=bf
			# ~ for a,s in zip(T,model.LE):self.Ab[s,t]=a
			# ~ T=[sum(self.Ab[source,t] * source.transit(target) for source in model.LE) for target in model.LE] # ~ transitions
		# ~ self.cpb=reduce(lambda x,y:x*y,self.BF.values())
		# ~ self.logcpb=sum(map(log,self.BF.values()))
		# ~ self.decal=abs(self.cpb-self.pobs)
	# ~ def betab(self,model):
		# ~ t=len(self.so)-1
		# ~ for s in model.LE:self.Bb[s,t]=1
		# ~ z=list(enumerate(self.so))
		# ~ z.reverse()
		# ~ for t,symb in z[:-1]:
			# ~ for source in model.LE:
				# ~ self.Bb[source,t-1]=sum(source.transit(target)*target.emit(symb)*self.Bb[target,t] for target in model.LE)/self.BF[t]
	# ~ def betabsw(self,model):
		# ~ t=len(self.so)-1
		# ~ for s in model.LE:self.Bb[s,t]=1/self.BF[t]
		# ~ z=list(enumerate(self.so))
		# ~ for t,symb in z[-1:0:-1]:
			# ~ for source in model.LE:
				# ~ self.Bb[source,t-1]=sum(source.transit(target)*target.emit(symb)*self.Bb[target,t] for target in model.LE)/self.BF[t-1]
	def gammab(self,model):
		self.alphab(model)
		self.betabsw(model)
		self.Gb={k:self.Ab[k]*self.Bb[k] * self.BF[k[1]] for k in self.Ab.keys()}
	def norm(self):
		for k in self.G.keys():self.G[k]/=self.cp
	def esttr(self,model):
		ntt={source:{target:0 for target in model.LE} for source in model.LE}
		sg={st:0 for st in model.LE}
		for t,o in enumerate(self.so[:-1]):
			for source in model.LE:
				sg[source]+=self.Gb[source,t]
		for source in model.LE:
			for target in model.LE:
				for t,observable in enumerate(self.so[1:]):
					ntt[source][target] += self.Ab[source,t] * source.transit(target) * target.emit(observable) * self.Bb[target,t+1]
		# ~ for source in model.LE:
			# ~ for target in model.LE:
				# ~ ntt[source][target]/=sumgamma[source]
		return ntt,sg
	def estpi(self,model):
		return {st:self.Gb[st,0] for st in model.LE}
	def estem(self,model):
		newte={st:{symb:0 for symb in model.alphabet} for st in model.LE}
		sg={st:0 for st in model.LE}
		for t,symb in enumerate (self.so):
			for source in model.LE:
				newte[source][symb]+=self.Gb[source,t]
				sg[source]+=self.Gb[source,t]
		return newte,sg
	def rep(self,D,model):
		print("\t".join(s.name for s in model.LE))
		for t in range(len(self.so)):
			print("\t".join(map(str,(D[s,t] for s in model.LE))))
		print("============================================")
if __name__ == '__main__':
	LS=[Seq("ABC"),Seq("BCD")]
	# ~ LS=[Seq("ABCDEF"),Seq("BCDEFG")]
	print(LS)
	model=mkrndmodel(LS)
	# ~ print(model.LM[-1].I.T)
	print("D0T",model.LM[0].D.T)
	print("M0T",model.LM[0].M.T)
	print("M0E",model.LM[0].M.E)
	print("I0T",model.LM[0].I.T)
	print("I0E",model.LM[0].I.E)
	print("D1T",model.LM[1].D.T)
	print("M1T",model.LM[1].M.T)
	print("M1E",model.LM[1].M.E)
	print("I1T",model.LM[1].I.T)
	print("I1E",model.LM[1].I.E)
	print (model.PI)
	# ~ nbobs=3
	# ~ nbmod=6
	# ~ for layer in range(nbmod+nbobs):
		# ~ for t in getlayers(layer,nbmod,nbobs):print(t)
	LS[0].alphab(model)
