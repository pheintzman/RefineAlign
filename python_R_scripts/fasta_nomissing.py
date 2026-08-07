# Script by Pontus Skoglund (https://github.com/pontussk/fasta_nomissing.py) distributed under a GNU General Public License v3.0

import sys
from optparse import OptionParser
import random
import math

usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--fastafile", action="store", type="string", dest="fastafile",help="fastafile")
parser.add_option("--maxmissing","--maxmissing_snp", action="store", type="float", dest="maxmissing",help="",default=1.0)
parser.add_option("--maxmissing_ind", action="store", type="float", dest="maxmissing_ind",help="",default=1.0)
parser.add_option("--polymorphic", action="store_true", dest="polymorphic",help="Keep only polymorphic sites",default=False)
parser.add_option("--notriallelic", action="store_true", dest="notriallelic",help="Exclude sites with 3 or more alleles",default=False)
parser.add_option("--notransitions", action="store_true", dest="notransitions",help="Exclude transition SNPs",default=False)

parser.add_option("--verbose", action="store_true", dest="verbose",help="",default=False)

(options, args) = parser.parse_args()


nucleotides=['A','T','G','C']
symbols=['A','T','G','C','N']

seqs=[]
IDs=[]
theseq=''
for line in open(options.fastafile):
	if '>' in line:
		add=True
		IDs.append(line.rstrip('\n').lstrip('>').split()[0])
		if len(IDs) > 1:
			seqs.append(theseq)
		theseq=''
	else:
		theseq += line.rstrip('\n').rstrip('\r')
seqs.append(theseq)


nsamples=len(seqs)
seqlength=min([len(seq) for seq in seqs])

print >>sys.stderr,'number of samples:',nsamples
print >>sys.stderr,'length before filters:',seqlength

mask=[]

newseqs=[]
for n in seqs:
	newseqs.append('')

print >>sys.stderr,'individual missing rates before filters:',
for ID,seq in zip(IDs,seqs):
	missingrate=1.0*seq.upper().count('N')/len(seq)
	print >>sys.stderr,ID+':'+str(round(missingrate,2)),
print >>sys.stderr,''

for x in xrange(0,seqlength):
	pileup=''
	for s in seqs:
		pileup += s[x].upper()
	if options.verbose: print pileup
	missingrate=1.0*pileup.count('N')/len(pileup)
	
	alleles=[a for a in pileup if a in nucleotides]
	numalleles=len(set(list(alleles)))
	masking=False
	
	
	notallowed=[a for a in pileup if a not in symbols]
	if len(notallowed) >0: 
		print >>sys.stderr,'unexpected character(s)',notallowed
		exit(0)
	
	if missingrate > options.maxmissing:
		masking=True
	if options.polymorphic and numalleles <2:
		masking=True
	if options.notriallelic and numalleles > 2:
		masking=True
	if options.notransitions:
		if 'C' in pileup and 'T' in pileup:continue
		if 'G' in pileup and 'A' in pileup:continue
	
	if masking==False:
		for p in xrange(0,nsamples):
			newseqs[p] += pileup[p]
	elif masking:
		mask.append(x)
		if options.verbose: print 'masked'
		continue
		

print >>sys.stderr,'filtering sites:',len(mask)
print >>sys.stderr,'keeping sites:',len(newseqs[0])

if len(newseqs[0]) == 0:
	print >>sys.stderr,'No sites remaining in alignment, exiting'

excluded_ind=0
for ID,newseq in zip(IDs,newseqs):
	ind_missingrate=1.0*newseq.count('N')/len(newseq)
	if ind_missingrate > options.maxmissing_ind: 
		excluded_ind +=1
		continue
	print '>'+ID
	print newseq
	
print >>sys.stderr,'individual missing rates after filters:',
for ID,newseq in zip(IDs,newseqs):
	missingrate=1.0*newseq.count('N')/len(newseq)
	print >>sys.stderr,ID+':'+str(round(missingrate,2))
print >>sys.stderr,''

print >>sys.stderr,'excluded individuals:',excluded_ind
exit(0)







for ID,seq in zip(IDs,seqs):
	print '>'+ID
	newseq=''
	for x,a in zip(xrange(0,len(seq)),seq):
		if x in mask:
			continue
		else:
			newseq += a
	print newseq

print >>sys.stderr,'length after filters:',len(newseq)
