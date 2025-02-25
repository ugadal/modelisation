import random
from math import log
from functools import reduce
def normd(D):
	t=sum(D.values())
	return {k:v/t for k,v in D.items()}
class state():
	def __init__(self,n="pas de nom"):
		# ~ print("initialisation d'une instance d'etat",n)
		self.name=n
		self.E={}
		self.T={}
	def emit(self,symb):return self.E.get(symb)
	def transit(self,target):return self.T.get(target)
class model():
	def __init__(self,LE,PI,alphabet):
		self.LE=LE
		self.PI=PI
		self.alphabet=alphabet
	def rep(self):
		print([(k.name,v) for k,v in self.PI.items()])
		for s in self.LE:
			print(s.name,[(k,v) for k,v in s.E.items()])
		for s in self.LE:
			print(s.name,[(k.name,v) for k,v in s.T.items()])
def makernddic(L):
	P=[random.random() for _ in L]
	t=sum(P)
	P=[p/t for p in P]
	return {o:p for o,p in zip(L,P)}
def mkrndmodel(nbe,lobservable):
	SL=[state(f"st-{x}") for x in range(nbe)]
	for s in SL:
		s.E=makernddic(lobservable)
		s.T=makernddic(SL)
	pi=makernddic(SL)
	return model(SL,pi,lobservable)
def printtrd(D):
	for k,d in D.items()		:
		print(k.name,[(t.name,v) for t,v in d.items()])
def printemd(D):
	for k,d in D.items()		:
		print(k.name,[(t,v) for t,v in d.items()])
class SO():
	def __init__(self,so,pobs):
		self.so=list(so)
		self.pobs=pobs
		self.Ab={}
		self.Bb={}
		self.Gb={}
		self.BF={}
	def alphab(self,model):
		T=[model.PI[s] for s in model.LE]
		for t,symb in enumerate(self.so):
			T=[t*s.emit(symb) for t,s in zip(T,model.LE)] # ~ emission
			bf=sum(T)
			T= [t/bf for t in T]
			self.BF[t]=bf
			for a,s in zip(T,model.LE):self.Ab[s,t]=a
			T=[sum(self.Ab[source,t] * source.transit(target) for source in model.LE) for target in model.LE] # ~ transitions
		self.cpb=reduce(lambda x,y:x*y,self.BF.values())
		self.logcpb=sum(map(log,self.BF.values()))
		self.decal=abs(self.cpb-self.pobs)
	def betab(self,model):
		t=len(self.so)-1
		for s in model.LE:self.Bb[s,t]=1
		z=list(enumerate(self.so))
		z.reverse()
		for t,symb in z[:-1]:
			for source in model.LE:
				self.Bb[source,t-1]=sum(source.transit(target)*target.emit(symb)*self.Bb[target,t] for target in model.LE)/self.BF[t]
	def gammab(self,model):
		self.alphab(model)
		self.betab(model)
		self.Gb={k:self.Ab[k]*self.Bb[k] for k in self.Ab.keys()}
	def norm(self):
		for k in self.G.keys():self.G[k]/=self.cp
	def rep(self,D,model):
		print("\t".join(s.name for s in model.LE))
		for t in range(len(self.so)):
			print("\t".join(map(str,(D[s,t] for s in model.LE))))
		print("============================================")
