

def read(cluster, tree, i):
	while i < len(cluster):
		x = cluster[i]
		
		if x == "(":
			newtree, i = read(cluster, [], i+1)
			tree.append(newtree)
			i+=1
		elif x == ")":
			return tree, i
		elif x == "," or x == ";" or x == " ":
			i+=1
		else:				
			s = ""
			while x != "," and x != "]" and x != "[" and x != ")" and x != "(" and x != ";" and x != " ":
				s+=x
				i+=1
				x = cluster[i]	
			ss = s.split(":")[0] #removing branch lenghts
			if ss != "": tree.append(ss)
	return tree, i



class Node:
	import copy
	import sys
	w=0
	def __init__(self, n, f=None, d=0):
		self.nr = Node.w
		self.src = n
		self.str = ""
		Node.w+=1
		self.fath = f
		self.map = None
		self.deep = d
		self.time = 0
		self.hgtStart = False
		self.hgtEnd = False
		self.leaf = False
		self.hgtLink=None
		self.sons = []

		if type(n)==str:
			
			if "hgt" in n:
				self.hgtStart = True
				self.label = n
				self.sons = []

			elif n.isupper():
				self.hgtEnd = True
				self.label = n
				self.sons=[]
			
			else:
				self.leaf = True
				self.label = n
			
		else:
			if len([i for i in n if (type(i) == str and "hgt" not in i and not i.isupper()) or type(i) == list])   > 1 and [i for i in n if type(i) == str and ( "hgt" in i or i.isupper())]:
				self.sons = []
				self.label = "n%d" %self.nr
				nodes = []
				i=0
				
				while i < len(n):
					nnodes = [n[i]]
					j=i+1
					while j < len(n) and type(n[j]) == str and (n[j].isupper() or "hgt" in n[j]):
						nnodes.append(n[j])
						j+=1
					i=j
					nodes.append(nnodes) if len(nnodes) > 1 else nodes.append(nnodes[0])
				self.sons = [Node(s, self, d+1)  for s in nodes]
				
			elif len([i for i in n if (type(i) == str and "hgt" not in i and not i.isupper()) or type(i) == list])   == 1:
				if True:	
					nodes = n[1:] + [n[0]]
				
					if "hgt" in nodes[0]:
						self.hgtStart = True
						self.label = nodes[0]
						self.hgtlabel = self.label.split("=")[1]
						
						sons = [nodes[-1]]+nodes[1:len(nodes)-1] if len(nodes) > 2 else nodes[1]
						self.sons=[Node(sons, self, d+1)]

					elif type(nodes[0]) == str and nodes[0].isupper():
						self.label = nodes[0]
						self.hgtEnd = True
						sons = [nodes[-1]]+nodes[1:len(nodes)-1] if len(nodes) > 2 else nodes[1]
						self.sons=[Node(sons, self, d+1)]
			

						
			else:
				self.label = "n%d" %self.nr	
				self.sons = [Node(s, self, d+1)  for s in n]
		self.kids=self.sons[:]
		

		
	@classmethod
	def readTree(cls,src):
		q=read(src.rstrip(), [], 0)[0][0]
		cls.w=0
		S=cls(q)
		S.connectHSTs()
		S.check_cycles()
		for n in S.nodes():	n.str = n.printTreeNewick([])
		return S
	
	def __str__(self):
		return self.label	

	def __repr__(self):
		return self.label	

	def siblingcluster(self):
		return list(set(self.fath.leaves()).difference(self.leaves()))	

	def strsiblingcluster(self):
		return " ".join("%s"%x for x in self.siblingcluster())

	def strcluster(self):
		return " ".join("%s"%x for x in set(self.leaves()))

		

	def connectHSTs(self):
		for n in self.hgtStarts():
			name = n.label.replace("hgt=", "")
			n.hgtLink = self.nodeByName(name)
			n.hgtLink.hgtLink = n
			n.kids.append(n.hgtLink)  
	
			
	def root(self):
		return self.nodeById(0) 																									
		
	def leaves(self):
		if self.leaf:	return [self]
		else: 
			l = []
			for s in self.sons: l+=s.leaves()
			return l

	def inner(self):
		if not self.leaf: 
			l = [self]
			for s in self.sons: 
				if not s.leaf: l+=s.inner()
			return l

	def nodes(self):
		if self.leaf: return [self]
		return self.leaves() + self.inner()

	def hgtStarts(self):
		l = [self] if self.hgtStart else []
		if not self.leaf: 
			for s in self.sons: l+=s.hgtStarts()
		return l

	def hgtEnds(self):
		l = [self] if self.hgtEnd else []
		if not self.leaf: 
			for s in self.sons: l+=s.hgtEnds()
		return l

	def HGTs_num(self):
		return len(self.hgtStarts())

	def all_edges(self):
		edges = {}
		for n in self.inner():
			for s in n.sons:
				edges[(n.label, s.label, False)] = (n, s)
		for n in self.hgtStarts():
			edges[(n.label, n.hgtLink.label, True)] = (n, n.hgtLink)
		return edges 

	def all_tree_edges(self):
		edges = {}
		for n in self.inner():
			for s in n.sons:
				if not (n.hgtStart and s.hgtEnd):
					edges[(n.label, s.label, False)] = (n, s)
		return edges

	def add_transfer(self, e1, e2, hgtName):

		nstart, nend = [self.nodeByName(o.label) for o in e1]
		newHgtStart = Node("hgt=%s" %hgtName, nstart, 0)
		newHgtStart.hgtlabel = newHgtStart.label.split("=")[1]
		newHgtStart.sons = [nend]

		
		nstart.sons.remove(nend)
		nstart.sons.append(newHgtStart)
		newHgtStart.fath = nstart
		nstart.kids = nstart.sons + [nstart.hgtLink] if nstart.hgtStart else nstart.sons
		nend.fath = newHgtStart


		nstart, nend = [self.nodeByName(o.label) for o in e2]
		newHgtEnd = Node("%s" %hgtName, nstart, 0)
		newHgtEnd.sons = [nend]
		

		nstart.sons.remove(nend)
		nstart.sons.append(newHgtEnd)
		newHgtEnd.fath = nstart
		nstart.kids =  nstart.sons + [nstart.hgtLink] if nstart.hgtStart else nstart.sons
		nend.fath = newHgtEnd

		newHgtStart.hgtLink = newHgtEnd
		newHgtEnd.hgtLink = newHgtStart
		newHgtStart.deep = newHgtStart.fath.deep
		newHgtEnd.deep = newHgtEnd.fath.deep

		newHgtStart.kids = newHgtStart.sons + [newHgtEnd]
		newHgtEnd.kids = newHgtEnd.sons

		for n in newHgtStart.nodes() + newHgtEnd.nodes(): n.deep += 1

	def generate_HGTs(self, name):
		import copy
		import sys
		trees = []
		c=0
		edges1 = self.all_tree_edges().values()
		for e in edges1:
			for e2 in edges1:
				if e != e2 :
					estart, e2start = self.nodeByName(e[0].label), self.nodeByName(e2[0].label)
					eend, e2end = self.nodeByName(e[1].label), self.nodeByName(e2[1].label)
					if (estart.hgtStart and e2start.hgtEnd and estart.label.replace("hgt=", "") == e2start.label) or (eend.hgtStart and e2end.hgtEnd and eend.label.replace("hgt=", "") == e2end.label): 
						continue
					newTree = copy.deepcopy(self)
					newTree.add_transfer(e, e2, name)
					try:
						newTree.check_cycles()
						for n in newTree.nodes(): n.str = n.printTreeNewick([])
						trees.append(newTree)
					except ValueError as err: 
						continue
					c+=1			
		return trees
			




	def lca(self, node):
		if self.leaf:
			for i in node.leaves():		
				if i.label==self.label or "?" in self.label:
					self.map=i
					break
		else:
			for s in self.sons: 
				if s.map==None: s.lca(node)
			r = [s.map for s in self.sons]
			
			min_r = min([(rr, rr.deep) for rr in r], key=lambda x:x[1])[0]
			for rr in r:
				while rr.deep > min_r.deep:
					rr = rr.fath
				while rr != min_r:
					rr = rr.fath
					min_r= min_r.fath

			self.map=min_r
	


	
		

	def printTree(self):
		if self.hgtStart:
			return "%s %s" %(self.sons[0].printTree(), self.label)
		elif self.hgtEnd:
			return "%s %s" %(self.sons[0].printTree(), self.label)
		elif self.leaf:
			return self.label + "#" + str(self.nr)
		else:
			return '('+' '.join([s.printTree() for s in self.sons])+')'

	def printTreeD(self):
		for n in self.nodes():
			print ">", n.label, n.deep

	def printTreeNewick(self, src):
		if self.hgtStart or self.hgtEnd:
			
			if self.sons[0].hgtStart or self.sons[0].hgtEnd: self.sons[0].printTreeNewick(src)
			else:	src.append(self.sons[0].printTreeNewick(src))
			
			src.insert(1, self.label)			
			return ' '.join(src)
		
	
		elif self.leaf:
			return self.label
		else:
			sons = sorted(self.sons, key=lambda s: s.hgtStart or s.hgtEnd)
			return '('+', '.join([s.printTreeNewick([]) for s in self.sons])+')'




	def printMap(self):
		print str(self.label)+'->'+str(self.map.label)
		if not self.leaf:
			for s in self.sons: s.printMap()


	def printMap1(self):
		if type(self.map) == list: print str(self.nr)+'-> !!!!',[i.nr for i in self.map]
		else: print str(self.nr)+'->'+str(self.map.nr)
		if not self.leaf:
			self.sons[0].printMap1()
			self.sons[1].printMap1()

	def printTime(self):
		if not self.leaf: 
			self.sons[0].printTime()
			self.sons[1].printTime()
		print "%d, %d" %(self.nr, self.time)

	def nodeById(self, nr):
		if self.nr == nr: return self
		if not self.leaf:
			node = None
			c = 0
			while not node and c < len(self.sons):
				node = self.sons[c].nodeById(nr)
				c += 1
		else: return None
		return node

	def nodeByName(self, label, fath=None):
		if self.label == label and ( not fath or self.fath.label == fath): return self
		if not self.leaf:
			node = None
			c = 0
			while not node and c < len(self.sons):
				node = self.sons[c].nodeByName(label, fath)
				c += 1
		else: return None
		return node

		

	def create_graph(self):
		def getnum(n):
			if n.hgtEnd: return n.hgtLink.nr
			return n.nr
		e=[]
		for n in self.nodes():
			a=getnum(n)
			for c in n.kids:
				b=getnum(c)
				if a<>b: e.append((a,b))
		return e
	

	def check_cycles(self):
		e = self.create_graph()
		l=[]
		d1={}
		d2={}
		for x,y in e:
			d1.setdefault(x,[]).append(y)
			d2.setdefault(y,[]).append(x)
		roots=list(set(d1.keys()).difference(d2.keys()))
		while roots:
			n=roots.pop()
			l.append(n)
			if n not in d1: continue
			ms=d1[n][:]
			for m in ms:
				d1[n].remove(m)
				d2[m].remove(n)
				if not d2[m]: 
					roots.append(m)
					d2.pop(m)
				if not d1[n]: d1.pop(n)
		if d1 or d2: raise ValueError('Cycle')
		return l



