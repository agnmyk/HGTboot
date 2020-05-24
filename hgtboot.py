from hgttree import Node
import hgtrecon
import sys
import itertools
import getopt
import copy
from glob import glob
from collections import defaultdict
import re
import random
import math

def print_help():
	print "\nUsage:"
	print "hgtboot.py -s <file with a species tree> -g <file with a gene tree> -g <file with bootstrap trees> [-d .., -l .., -t..]"
	print "\nOptions:"
	print "\t-s <file with a species tree>"
	print "\t-g <file with a gene tree>"
	print "\t-g <file with bootstrap trees> -- one tree in each line."
	print "\t-d <cost of duplication event>  -- default value = 2"
	print "\t-l <cost of loss event> -- default value = 1"
	print "\t-t <cost of horizontal gene transfer event> -- default value = 3"
	print "\t-h -- print help"
	sys.exit()

def read_data():

	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:g:b:t:d:l:i:h")
	except getopt.GetoptError as err:
		print "Wrong option"
		print_help()

	

	for o, a in opts:
		if o == "-s": stree = a
		elif o == "-g": gtree = a
		elif o == "-b": bootstraptrees = a
		elif o == "-d": 
			global DUP
			DUP = int(a)
		elif o == "-t": 
			global HGT
			HGT = int(a)
		elif o == "-l": 
			global LOSS
			LOSS = int(a)
		elif o == "-i": 
			global iterations_num 
			iterations_num = int(a)
		elif o == "-h": 
			print_help()
		
	if not any("-s" in a for a in opts): 
		print "No species tree given\n\n"
		print_help()

	if not any("-g" in a for a in opts): 
		print "No gene tree given\n\n"
		print_help()

	if not any("-b" in a for a in opts): 
		print "Bootstrap trees file not provided\n\n"
		print_help()

	streeFin = open(stree)
	S = Node.readTree(streeFin.readline().rstrip())
	gtreeFin = open(gtree)
	G = Node.readTree(gtreeFin.readline().rstrip())
	G.lca(S)

	boottrees = []
	for i in open(bootstraptrees):
		if i != "": boottrees.append(Node.readTree(i))

	return G, S, boottrees



def hgtboot(G, S, boottrees):
	new_HGT_name = "P"+str(S.HGTs_num())
	newtrees = S.generate_HGTs(new_HGT_name)

	transfer_boot_supp = {k:defaultdict(int) for k in xrange(len(newtrees))}
	transfer_supp = {k:defaultdict(int) for k in xrange(len(newtrees))}
	optsols = []

	c = 0


	solutions_dict = {k:defaultdict(str) for k in xrange(iterations_num)}
	transfer_boot_supp = {k:defaultdict(str) for k in xrange(iterations_num)}
	transfer_supp = {k:defaultdict(str) for k in xrange(iterations_num)}
	optsols = {k:{} for k in xrange(iterations_num)}
	solnum = {k:{} for k in xrange(iterations_num)}
	solnum_boot = {k:defaultdict(str) for k in xrange(iterations_num)}
	optGTree = None

	def boot(S, S_id, itr):
		best_newtrees = []
		best_optsol = float("inf")

		new_HGT_name = "P"+str(S.HGTs_num())
		newtrees = S.generate_HGTs(new_HGT_name)

		c=0
		for newS in newtrees:
			solutions, optsol = hgtrecon.opt_reconstruct(G,newS,DUP,HGT,LOSS)
			optsols[itr][newS] = optsol
			trfsols = False
			trfsols_dict = {k:0 for k in [h.label for h in newS.hgtEnds()]}

		
			solnum[itr][newS]= len(solutions)
		
			for sol in solutions:
				trf_counted = {k:0 for k in [h.label for h in newS.hgtEnds()]}
				trfs = re.findall("\^([^)(,]*)", sol)
				x = new_HGT_name in trfs
				if x: trfsols = True
				for tt in trfs:
					if not trf_counted[tt]:
						trfsols_dict[tt] += 1
					trf_counted[tt] = 1

			if trfsols: 
				if newS not in transfer_supp[itr]: transfer_supp[itr][newS] = None

				transfer_supp[itr][newS] = (S_id, c, copy.deepcopy(trfsols_dict))
				if optsol < best_optsol:
					best_optsol = optsol
					best_newtrees = [c]
					solutions_dict[itr][newS] = solutions
				elif optsol == best_optsol:
					best_newtrees.append(c)
					if newS in solutions_dict[itr]: 
						solutions_dict[itr][newS].append(solutions[0])
					else: 
						solutions_dict[itr][newS] = solutions
			
			c+=1	
		
		for c in best_newtrees:
		
			newS = newtrees[c]
			for bG in boottrees:
				bootsolutions, optsol_boot = hgtrecon.opt_reconstruct(bG,newS,DUP,HGT,LOSS)
				trfopt = False
				trfsols_dict = {k:0 for k in [h.label for h in newS.hgtEnds()]}
				for sol in bootsolutions:
					trfs = re.findall("\^([^)(,]*)", sol)
					x = new_HGT_name in trfs
					if x: trfopt = True
					for tt in trfs:
						trfsols_dict[tt] = 1

				if trfopt: 
					if newS not in transfer_boot_supp[itr]: transfer_boot_supp[itr][newS] = defaultdict(int)
					for tt in trfsols_dict:
						transfer_boot_supp[itr][newS][tt] += trfsols_dict[tt]
		
			if itr < iterations_num-1: boot(newS, c, itr+1)

	boot(S,0,0)

	print
	for i in transfer_supp:
		print "=== ITERATION", i, "==="
		for t in transfer_supp[i]:
			if transfer_supp[i][t] and transfer_boot_supp[i][t]:
				print "Species tree S':", t.printTreeNewick([])
				print "Optimal cost:", optsols[i][t]
				print "Number of optimal scenarios for S' and gene tree G:", solnum[i][t]
				print "Transfers usage:", " ".join(k+":"+str(transfer_supp[i][t][-1][k]) for k in transfer_supp[i][t][-1])
				print "Transfers support:", " ".join(k+":"+str(round(transfer_boot_supp[i][t][k]/float(len(boottrees)),2)) for k in sorted(transfer_boot_supp[i][t]))
				print "Optimal scenarios:", "\n".join(sol for sol in solutions_dict[i][t])

				print
		print
		print



iterations_num = 2
DUP, HGT, LOSS = 2, 3, 1
G, S, boottrees = read_data()
hgtboot(G, S, boottrees)

