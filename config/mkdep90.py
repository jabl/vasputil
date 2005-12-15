#!/usr/bin/env python
#
# $Id$
#

import sys, string, re

use=re.compile('\s*use\s+(\w+)', re.I)
mod=re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)

class depfile:
	def __init__(self, name):
		self.name=name
		self.oname=name[:-3]+'o'
		self.uses={}
	
	def adduses(self, u):
		self.uses[u]=''
		
def main():
	allmods={}
	dfiles=[]
	for ff in sys.argv[1:]:
		fd=open(ff, 'r')
		buf=fd.readlines()
		fd.close()
		dfiles.append(depfile(ff))

		for ln in buf:
			m=mod.match(ln)
			if m is not None:
				allmods[m.group(1)]=ff
			else:
				m=use.match(ln)
				if m is not None:
					dfiles[-1].adduses(m.group(1))
	for df in dfiles:
		deps=[]
		for dd in df.uses.keys():
			try:
				if (allmods[dd] != df.name):
					omod=allmods[dd][:-3]+'o'
					deps.append(omod)
			except:
				print >> sys.stderr, 'Missing dependency for', dd, 'in',\
				df.name
		if deps:
			dstr=df.oname+': '
			for i in deps:
				dstr=dstr+i+' '
			print dstr


if __name__ == '__main__':
	main()
