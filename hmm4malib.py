import random
from math import log
from functools import reduce

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
def getlayers(l,m,t):
	for mod in range(l):
		if mod>m:break
		if l-mod>t:break
		yield mod,l-mod
def normd(D):
	trsh=1e-5
	t=sum(D.values())
	new={k:v/t for k,v in D.items()}
	if any(v==1.0 for v in new.values()):
		new={k:0 if v!=1.0 else 1.0 for k,v in new.items()}
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
		self.to_mod=None
		self.from_mod=None
		
class model():
	def __init__(self,LM,PI,LSEQ):
		self.LM=LM
		self.PI=PI
		self.LSEQ=LSEQ
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
	P=[random.random() for _ in L]
	t=sum(P)
	P=[p/t for p in P]
	return {o:p for o,p in zip(L,P)}
def mkrndmodel(LSEQ):
	ttl=sum(seq.L for seq in LSEQ)
	alphabet=set()
	for seq in LSEQ:alphabet=alphabet or set(seq.so)
	LM=[module(x,alphabet) for x in range(ttl+1)]
	for source,dest in zip(LM,LM[1:]):
		source.M.T=makernddic(list("IMD"))
		source.D.T=makernddic(list("IMD"))
		source.I.T=makernddic(list("IMD"))
	final=module(ttl+1,alphabet)
	final.I.T=makernddic(list("IMD"))
	final.I.T["I"]=1.0
	final.I.T=normd(final.I.T)
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
	# ~ def alphab(self,model):
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
	LS=[Seq("ABC"),Seq("BCD"),Seq("BC")]
	print(LS)
	mkrndmodel(LS)
	
