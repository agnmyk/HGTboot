#!/usr/bin/python

from collections import defaultdict
import getopt
import sys
import itertools
import time

LOSS = 1
DUP = 1#2
HGT = 1#3

deltatreev={}
deltav = {}
deltaupv = {}

deltares={}
deltatreeres={}
deltaupres={}
deltauprespath={}

transfers=defaultdict(lambda : defaultdict(int))

recon={}

INFTY=1e100
reconstruct = True




def losscount(s,x):		 
	if (not s.hgtStart and not s.hgtEnd) or hgtcount(s,x): return len(s.kids)-1
	return 0

def hgtcount(s,x): return int(s.hgtStart and x==s.hgtLink)


def locate_min(a):
    smallest = min(a)
    return smallest, [index for index, element in enumerate(a) 
                      if smallest == element]

def delta(g,s):
	if (g,s) in deltav: 
		return deltav[g, s]

	if g.leaf and s.leaf: 
	
		if g.map.label==s.label: res=0
		else: res=INFTY
		
		if reconstruct:
			deltares[g, s]=[]

	else:
		alphav=betav=gammav=INFTY*3
		if not s.hgtStart and not s.hgtEnd:
			alphav = alpha(g,s)
			betav = beta(g,s)
			res = min(alphav, betav)
		elif s.hgtStart:
			betav = beta(g,s)
			gammav = gamma(g,s)
			res = min(betav, gammav)

		else: res = INFTY
		deltav[g, s]=res

		if reconstruct:
			start = time.time()
			deltares[g, s]=[]
			if not s.hgtStart and not s.hgtEnd:
				if alphav==res:
					for pp in [p for p in [zip(g.sons,x) for x in itertools.permutations(s.sons,2)]]:

						if res == (len(s.sons)-len(g.sons))*LOSS + sum([deltaup(*ppp) for ppp in pp]):
							for x,y in itertools.product(deltaupres[pp[0][0], pp[0][1]],deltaupres[pp[1][0],pp[1][1]]):
								deltares[g, s].append((x,y,pp[0][1],pp[1][1],"SPEC"))
			if not s.hgtEnd:
				if betav == res and not g.leaf:
					if (g,s) not in recon: recon[g, s] = []
					pairs = [ zip(g.sons, [s,s]) ] if s.leaf else [p for p in [zip(g.sons, x) for x in set(itertools.permutations(s.kids+[s,s],2)) if s in x]]
					for pair in pairs:
						ssum=0
						for p in pair:
							ssum+= delta(*p) if p[1] == s else deltaup(*p) + LOSS*losscount(s,p[1])+HGT*hgtcount(s,p[1])
						if res == DUP + ssum:
							if pair[0][1] == s and pair[1][1] != s:
								for x,y in itertools.product([s],deltaupres[pair[1][0],pair[1][1]]):
									deltares[g, s].append((x,y,pair[0][1],pair[1][1],"DUP"))
							elif pair[1][1] == s and pair[0][1] != s:
								for x,y in itertools.product(deltaupres[pair[0][0], pair[0][1]],[s]):
									deltares[g, s].append((x,y,pair[0][1],pair[1][1],"DUP"))
							elif pair[0][1] == s and pair[1][1] == s:
								for x,y in itertools.product([s],[s]):
									deltares[g, s].append((x,y,pair[0][1],pair[1][1],"DUP"))
							
			if s.hgtStart:
				
				if gammav == res and not g.leaf:
					pairs = [p for p in [zip(g.sons,x) for x in itertools.permutations(s.kids,2)]]
					for pair in pairs:
						if res == HGT + deltaup(*pair[0]) + deltaup(*pair[1]):
							for x,y in itertools.product(deltaupres[pair[0][0],pair[0][1]],deltaupres[pair[1][0],pair[1][1]]):
								deltares[g, s].append((x,y,pair[0][1],pair[1][1],"HGT"))


	end = time.time()
	return res


def alpha(g,s):
	if s.leaf or g.leaf: 
		return INFTY
	pairs = [p for p in [zip(g.sons,x) for x in itertools.permutations(s.sons,2)]]
	return (len(s.sons)-len(g.sons))*LOSS + min([sum([deltaup(*ppp) for ppp in pp]) for pp in pairs])




def beta(g,s):
	if g.leaf: return INFTY
 
	pairs = [ zip(g.sons, [s,s]) ] if s.leaf else [p for p in [zip(g.sons, x) for x in set(itertools.permutations(s.kids+[s,s],2)) if s in x]]		
	sums = []
	for pair in pairs:
		ssum = 0
		for p in pair:
			ssum+= delta(*p) if p[1] == s else deltaup(*p) + LOSS*losscount(s,p[1])+HGT*hgtcount(s,p[1])
		sums.append(ssum)
	return  DUP+min(sums)

def gamma(g,s):	
	if g.leaf: return INFTY
	r= HGT + min([deltaup(*pair[0])+deltaup(*pair[1]) for pair in [p for p in [zip(x,s.sons+[s.hgtLink]) for x in itertools.permutations(g.sons,2)]]])
	return r




def deltaup(g,s):
	start = time.time()
	
	if (g,s) not in deltaupv:
		
		if s.leaf:
			deltaupv[g, s] = res = delta(g,s)
		else:
			deltaupv[g, s] = res = min(INFTY,delta(g,s), 
				min( LOSS*losscount(s,x) + HGT*hgtcount(s,x) + deltaup(g,x) for x in s.kids))
		if reconstruct: 
			deltaupres[g, s]=[]
			if res!=INFTY: 
				if res==delta(g,s): 
					deltaupres[g, s].append(s)
					
				if not s.leaf:
						for x in s.kids:						
							if res==LOSS*losscount(s,x) + HGT*hgtcount(s,x) + deltaup(g,x):
								deltaupres[g, s].extend( deltaupres[g, x] )

	end = time.time()
	return deltaupv[g, s]

def deltatree(g,s):
		if (g,s) not in deltatreev: 
			deltatreev[g, s]=res1=min(delta(g,c) for c in s.nodes())		

		return deltatreev[g, s]


def pt(n): return n.printTreeNewick([])

deltaoptpath={}
optpathres={}

def optpath(s,s1): # opt. sciezki od s do s1
	if (s,s1) in deltaoptpath: return deltaoptpath[s,s1]
	if s1==s: 
		deltaoptpath[s,s1]=0
	elif s.leaf: 
		deltaoptpath[s,s1]=INFTY
	else:
		deltaoptpath[s,s1]=min(INFTY,  min( LOSS*losscount(s,x) + HGT*hgtcount(s,x) + optpath(x,s1) for x in s.kids))
	
	optpathres[s,s1]=[]
	if deltaoptpath[s,s1]!=INFTY:
		if s==s1: optpathres[s,s1].append(None)
		else:
			optpathres[s,s1]=[x for x in s.kids if deltaoptpath[s,s1]==LOSS*losscount(s,x) + HGT*hgtcount(s,x) + optpath(x,s1) ]
	return deltaoptpath[s,s1]
		


def _gensolution(g,s):
		def addlosses(s1,s,t,dupcase=0):
			while (not dupcase and s1 not in s.sons) or (dupcase and s1!=s):
				t="(%s,%s-)~"%(t,s1.strsiblingcluster())
				s1=s1.fath
			return t
		def _genpath(s1,s,t):
			optpath(s,s1)
			if s1==s: yield t
			else:
				for x in optpathres[s,s1]:
					for i in _genpath(s1,x,t):
						if hgtcount(s,x):
							yield "(%s-,%s^%s)"%(s.strcluster(),i,s.hgtlabel)
						elif s.hgtEnd or s.hgtStart:
							yield i
						else:
							yield "(%s-,%s)"%(x.strsiblingcluster(),i)

					

		if g.leaf: 
			yield "%s"%g
		else:
			g1,g2=g.sons
			for x,y,x1,y1,tp in deltares[g, s]:
				if tp=='SPEC': 
					for t1,t2 in itertools.product(_gensolution(g1,x),_gensolution(g2,y)):
						for p1,p2 in itertools.product(_genpath(x,x1,t1),_genpath(y,y1,t2)):

							yield "(%s,%s)"%(p1,p2)					

				elif tp=="DUP": 

					for t1,t2 in itertools.product(_gensolution(g1,x),_gensolution(g2,y)):						
						for p1,p2 in itertools.product(_genpath(x,s,t1),_genpath(y,s,t2)):
							
							yield "(%s,%s)"%(p1,p2)
				elif tp=="HGT": 

					for t1,t2 in itertools.product(_gensolution(g1,x),_gensolution(g2,y)):						
						for p1,p2 in itertools.product(_genpath(x,x1,t1),_genpath(y,y1,t2)):
							if hgtcount(s,x1):
								yield "(%s,%s^%s)"%(p2,p1,s.hgtlabel)
							else:
								yield "(%s,%s^%s)"%(p1,p2,s.hgtlabel)


def clear_dict(d):
	for k in d.keys(): del d[k]

def find_opt_sol(G,S):
	clear_dict(deltatreev)
	clear_dict(deltav)
	clear_dict(deltaupv)

	clear_dict(deltares)
	clear_dict(deltatreeres)
	clear_dict(deltaupres)
	clear_dict(deltauprespath)

	for k in transfers: transfers[k] = defaultdict(int)
	G.lca(S)
	deltatree(G.root(), S.root())

	return min(deltav[G.root(),s] for s in S.nodes())



def opt_reconstruct(G,S,d=2,l=1,t=3):
	import sys
	global DUP, HGT, LOSS
	DUP, HGT, LOSS = d, l, t
	
	optsol = find_opt_sol(G,S)
	solutions = []
	for s in S.nodes():
		if optsol==deltav[G.root(),s]:
			for t in _gensolution(G.root(), s):
				solutions.append(t)
	return solutions, optsol

def transfers_usage():
	return transfers



