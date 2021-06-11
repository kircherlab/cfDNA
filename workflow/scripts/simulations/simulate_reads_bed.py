#!/usr/bin/env python

"""
:Author: Martin Kircher, Sebastian Roener
:Contact: martin.kircher@bihealth.de
:Date: *20.10.2020
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
from sortedcontainers import SortedDict
import random
import mmap
import gzip
import pysam

def read_genome_fastaindex(faifile):
  #READ FASTA INDEX TO MEMORY
  fastaindex = {}
  if os.path.exists(faifile):
    infile = open(faifile)
    for line in infile:
      fields = line.split()
      if len(fields) == 5:
        cname,length,start,line,cline = fields
        fastaindex[cname]=int(length),int(start),int(line),int(cline)
      else:
        sys.stderr.write('Error: Unexpected line in fasta index file: %s\n'%(line.strip()))
        sys.exit()
    infile.close()
  else:
    sys.stderr.write('Error: Fasta index file not available: %s\n'%faifile)
    sys.exit()
  return fastaindex

def get_DNA_file(chrom,start,end,fastafile,fastaindex):
  #try:
    length,bstart,bline,cline = fastaindex[chrom]
    bases = ""
    for pos in range(start,end+1):
      if pos <= length and pos >= 1:
        hpos = pos -1
        npos = (hpos // bline)*cline+(hpos % bline)
        fastafile.seek(bstart+npos)
        bases+=fastafile.read(1)
    return bases
  #except:
    #return None

def DNAiterator(chrom,start,fastamap,fastaindex):
  try:
    length,bstart,bline,cline = fastaindex[chrom]
    pos = start
    while (pos <= length and pos >= 1):
      hpos = pos -1
      npos = (hpos // bline)*cline+(hpos % bline)
      fastamap.seek(bstart+npos)
      yield fastamap.read(1)
      pos += 1
  except:
    pass
  raise StopIteration

def writeBAMentry(chrom,pos,seq,length,strand):
  global BAMoutfile,rcounter
  global options

  chromTID = BAMoutfile.gettid(chrom)
  if chromTID == -1: return
  
  qual = "I"*len(seq)
  
  if length > options.readlength:
    forward = pysam.AlignedRead()
    forward.qname = "SIM-%s-%09d"%(chrom,rcounter)
    forward.is_paired = True
    forward.is_proper_pair = True
    forward.is_read1 = True
    forward.tid = chromTID
    forward.rnext = chromTID

    reverse = pysam.AlignedRead()
    reverse.qname = "SIM-%s-%09d"%(chrom,rcounter)
    reverse.is_paired = True
    reverse.is_proper_pair = True
    reverse.is_read2 = True
    reverse.tid = chromTID
    reverse.rnext = chromTID
    
    if strand:
      forward.seq = seq[:options.readlength]
      forward.qual = qual[:options.readlength]
      forward.pos = pos-1
      forward.mpos = pos+length-options.readlength-1
      forward.isize = -length
      forward.mate_is_reverse = True

      reverse.seq = seq[-options.readlength:]
      reverse.qual = qual[-options.readlength:]
      reverse.is_reverse = True
      reverse.pos = pos+length-options.readlength-1
      reverse.mpos = pos-1
      reverse.isize = length
    
    else:
      forward.seq = seq[-options.readlength:]
      forward.qual = qual[-options.readlength:]
      forward.is_reverse = True
      forward.pos = pos-options.readlength
      forward.mpos = pos-length
      forward.isize = -length
      
      reverse.seq = seq[:options.readlength]
      reverse.qual = qual[:options.readlength]
      reverse.pos = pos-length
      reverse.mpos = pos-options.readlength
      reverse.isize = length
      reverse.mate_is_reverse = True

    forward.mapq = 255
    forward.cigar = [(0,options.readlength)]
    reverse.mapq = 255
    reverse.cigar = [(0,options.readlength)]
    BAMoutfile.write(forward)
    BAMoutfile.write(reverse)
    rcounter += 1
  else:
    forward = pysam.AlignedRead()
    forward.is_reverse = not strand
    forward.qname = "M_SIM-%s-%09d"%(chrom,rcounter)
    forward.seq = seq
    forward.qual = qual
    forward.tid = chromTID
    if strand: forward.pos = pos-1
    else: forward.pos = pos-length
    forward.mapq = 255
    forward.cigar = [(0,length)]
    forward.mrnm = -1
    forward.mpos = -1
    
    BAMoutfile.write(forward)
    rcounter += 1

  #print chrom,pos,length,seq,len(seq)


def readKmerGenome(filename):
  global options
  infile = open(filename)
  total = 0 
  minCount,maxCount = None,None
  cLength = None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
      if cLength == None:
        cLength = len(seq)
      total += count
  infile.close()
  if options.verbose: sys.stderr.write("total: %d/min: %d/max: %d\n"%(total,minCount,maxCount))
  total = float(maxCount)
  infile = open(filename)
  res = defaultdict(int)
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      res[seq] = count/total
  infile.close()
  return cLength,res


def readKMers(filename,reference):
  global options
  infile = open(filename)
  total = 0
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
  infile.close()
  if options.verbose: sys.stderr.write("total: %d/min: %d/max: %d\n"%(total,minCount,maxCount))
  total = float(maxCount)
  infile = open(filename)
  maxCount = None
  res = defaultdict(int)
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      freq = (count/total)/reference[seq]
      if maxCount == None: 
        maxCount = freq
      if maxCount < freq: maxCount = freq
      res[seq] = freq
  infile.close()
  if options.verbose: sys.stderr.write("maxCount: %d\n"%(maxCount))
  maxCount = float(maxCount)
  for key,value in res.items():
    res[key]=value/maxCount
  return res


def get_coords(bedfile):
    # BEDCOLS = ["chrom", "chromStart", "chromEnd", "name"]
    res = []
    infile = open(bedfile)
    for line in infile:
      if line.startswith("#"): continue
      fields = line.strip().split("\t")
      if len(fields) > 2:
        res.append( (fields[0].replace("chr", ""), int(fields[1]), int(fields[2]) ) )
    return res


parser = OptionParser()
parser.add_option(
    "-o",
    "--outfile",
    dest="outfile",
    help="Name of output file (def output.bam)",
    default="output.bam",
)
parser.add_option(
    "-r",
    "--region",
    dest="region",
    help="Regions to simulate (def regions.bed)",
    default="regions.bed",
)
parser.add_option(
    "-f",
    "--fasta",
    dest="reference",
    help="Fasta indexed reference genome (default hg19.fa)",
    default="hg19.fa",
)
parser.add_option(
    "-p",
    "--pipe",
    dest="pipe",
    help="Do not fetch, stream data (def Off)",
    default=False,
    action="store_true",
)
parser.add_option(
    "-d",
    "--lengthDist",
    dest="lengthDist",
    help="Length distribution (default CH01_lenDist.tsv )",
    default="CH01_lenDist.tsv",
)
parser.add_option(
    "--fwdKMerGenome",
    dest="fwdKMerGenome",
    help="Frequency distribution of KMers for left reads ends in the genome (default GRCh37_regChroms_2mers.tsv)",
    default="GRCh37_regChroms_2mers.tsv",
)
parser.add_option(
    "--revKMerGenome",
    dest="revKMerGenome",
    help="Frequency distribution of KMers for right reads ends in the genome (default GRCh37_regChroms_2mers.tsv)",
    default="GRCh37_regChroms_2mers.tsv",
)
parser.add_option(
    "--fwdPKMers",
    dest="fwdPKMers",
    help="Frequency distribution of KMers at left reads ends on plus strand (default CH01_left_2mer_f.tsv)",
    default="CH01_left_2mer_f.tsv",
)
parser.add_option(
    "--fwdMKMers",
    dest="fwdMKMers",
    help="Frequency distribution of KMers at left reads ends on minus strand (default CH01_left_2mer_r.tsv)",
    default="CH01_left_2mer_r.tsv",
)
parser.add_option(
    "--revPKMers",
    dest="revPKMers",
    help="Frequency distribution of KMers at right reads ends on plus strand (default CH01_right_2mer_f.tsv)",
    default="CH01_right_2mer_f.tsv",
)
parser.add_option(
    "--revMKMers",
    dest="revMKMers",
    help="Frequency distribution of KMers at right reads ends on minus strand (default CH01_right_2mer_r.tsv)",
    default="CH01_right_2mer_r.tsv",
)
parser.add_option(
    "-l",
    "--readlength",
    dest="readlength",
    help="Read length (default 45)",
    default=45,
    type="int",
)
parser.add_option(
    "-s",
    "--sample",
    dest="sample",
    help="Number of times to sample from a position (default 30)",
    default=30,
    type="int",
)
parser.add_option(
    "-c",
    "--correction",
    dest="correction",
    help="Likelihood correction factor to increase chances of piking a putative alignment (default 1.0)",
    default=1.0,
    type="float",
)
parser.add_option(
    "-v",
    "--verbose",
    dest="verbose",
    help="Turn debug output on",
    default=False,
    action="store_true",
)
(options, args) = parser.parse_args()

#############################
# READ LENGTH DISTRIBUTION
#############################

lengthDist = SortedDict()
if os.path.exists(options.lengthDist):
  infile = open(options.lengthDist)
  line = infile.readline() # Skip header
  total = 0
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      length,count = map(int,fields)
      total += count
  infile.close()
  infile = open(options.lengthDist)
  line = infile.readline() # Skip header
  total = float(total)
  rsum = 0
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      length,count = map(int,fields)
      rsum+=count
      lengthDist[rsum/total] = length
  infile.close()
else:
  sys.stderr.write("Error: Could not open length distribution!\n")
  sys.exit()

minLength = min(lengthDist.values())
maxLength = max(lengthDist.values())

#############################
# READ KMER DISTRIBUTIONS
#############################

forwardGenome = defaultdict(int)
reverseGenome = defaultdict(int)
forwardLength = None
reverseLength = None

if os.path.exists(options.fwdKMerGenome):
  forwardLength, forwardGenome = readKmerGenome(options.fwdKMerGenome)
else:
  sys.stderr.write("Error: Could not open left read end kmer distribution of genome!\n")
  sys.exit()

if os.path.exists(options.revKMerGenome):
  reverseLength, reverseGenome = readKmerGenome(options.revKMerGenome)
else:
  sys.stderr.write("Error: Could not open right read end kmer distribution of genome!\n")
  sys.exit()

forward = { True:defaultdict(int), False:defaultdict(int) }
reverse = { True:defaultdict(int), False:defaultdict(int) }

if os.path.exists(options.fwdPKMers):
  forward[True]=readKMers(options.fwdPKMers,forwardGenome)
else:
  sys.stderr.write("Error: Could not open left read end kmer [plus strand] distribution!\n")
  sys.exit()
  
if os.path.exists(options.fwdMKMers):
  forward[False]=readKMers(options.fwdMKMers,forwardGenome)
else:
  sys.stderr.write("Error: Could not open left read end kmer [minus strand] distribution!\n")
  sys.exit()

if os.path.exists(options.revPKMers):
  reverse[True]=readKMers(options.revPKMers,reverseGenome)
else:
  sys.stderr.write("Error: Could not open right read end kmer [plus strand] distribution!\n")
  sys.exit()
  
if os.path.exists(options.revMKMers):
  reverse[False]=readKMers(options.revMKMers,reverseGenome)
else:
  sys.stderr.write("Error: Could not open right read end kmer [minus strand] distribution!\n")
  sys.exit()

if options.verbose:
  sys.stderr.write("""
  Forward kmer length: %d
  Reverse kmer length: %d
  
  Initialized kmer sampling rates:"
  FOW+ %s %d
  FOW- %s %d
  REV+ %s %d
  REV- %s %d\n"""%(forwardLength,
                 reverseLength,
                 str(forward[True])[:100],sum(forward[True].values()),
                 str(forward[False])[:100],sum(forward[False].values()),
                 str(reverse[True])[:100],sum(reverse[True].values()),
                 str(reverse[False])[:100],sum(reverse[False].values())))

if not os.path.exists(options.reference) or not os.path.exists(options.reference+".fai"):
  sys.stderr.write("Error: Could not open reference fasta file or corresponding fasta index!\n")
  sys.exit()
  
genome_index = read_genome_fastaindex(options.reference+".fai")
genome_map = open(options.reference,'r')
#genome_map = get_fasta_mmap(options.reference)
#genome = pysam.Fastafile(options.reference)
#print(genome_index)

#################
# PARSE REGIONS
#################

options.region = options.region.replace('"', "").replace("'", "").strip()
chromosomes = []
if os.path.exists(options.region):
  chromosomes = get_coords(options.region)
else:
  sys.stderr.write("Error: Could not read genomic regions (%s)! Using all (major) chromosomes in reference fasta.\n"%(options.region))
  for chrom,(length,pos,lenLine1,lenLine2) in genome_index.items():
    if chrom.startswith("GL") or chrom.startswith("NC_") or chrom in ["MT","hs37d5"]: continue
    chromosomes.append((chrom,1,length))
chromosomes.sort()

helper = []
chrPrefix = ""
for chrom,(length,pos,lenLine1,lenLine2) in genome_index.items():
  if chrom.startswith("chr"): chrPrefix = "chr"
  helper.append({'LN': length, 'SN': chrom})
  
rcounter = 1
BAMoutfile = pysam.Samfile(options.outfile, "wb", header={'SQ': helper, 'HD': {'SO': 'unknown', 'VN': '1.4'}})

# Empirically identified factor for required iterations to achieve required coverage
targetCov = abs(options.sample)
#if options.correction < 1.0 and abs(int(options.sample*options.correction)) > 1:
 #targetCov = abs(int(options.sample*options.correction)+1)
 #if targetCov < options.sample*options.correction:
   #options.correction = 1.0

correctForward = options.correction
if ((4**forwardLength)/float(4**reverseLength) > options.correction and options.correction > 1):
  correctReverse = 1.0
else:
  correctReverse = options.correction/((4**forwardLength)/float(4**reverseLength))

if options.verbose:
  sys.stderr.write("Sampling, FwdCorrect, RevCorrect: %d / %d / %d\n"%(targetCov,correctForward,correctReverse))

count = 0
for chrom,start,end in chromosomes:
  if (chrom not in genome_index) and (chrPrefix+chrom in genome_index):
    chrom=chrPrefix+chrom
  if chrom not in genome_index:
    sys.stderr.write("Warning: Chromosomes unavailable, skipping region %s:%d-%d...\n"%(chrom,start,end))
  if options.verbose:
    sys.stderr.write("Reading %s:%d-%d...\n"%(chrom,start,end))
  
  chromLength=genome_index[chrom][0]
  context = get_DNA_file(chrom,max(1,start-maxLength),min(chromLength,start+maxLength),genome_map,genome_index).upper()
  posInWindow = maxLength
  if (start-maxLength) < 1:
    posInWindow = maxLength-(maxLength-start)-1
  #print posInWindow,len(context),context[posInWindow-25:posInWindow+25]
  for pos in range(start,end+1):
    #print chrom,pos,context[posInWindow],posInWindow,get_DNA_file(chrom,pos,pos,genome_map,genome_index)

    posStrand = []
    # CAN ONLY EXTRACT PLUS STRAND
    if (posInWindow < maxLength) and (posInWindow+maxLength <= len(context)): 
      posStrand = [True]
    # CAN ONLY EXTRACT MINUS STRAND
    elif (posInWindow >= maxLength) and (posInWindow+maxLength > len(context)): 
      posStrand = [False]
    # CAN THEORETICALLY EXTRACT FROM BOTH STRANDS
    else: posStrand = [True,False]
    
    for strand in posStrand:
      kmer = ""
      if strand: kmer = context[posInWindow:posInWindow+forwardLength]
      else: kmer = context[posInWindow-forwardLength+1:posInWindow+1]
      
      for i in range(targetCov//len(posStrand)):
        # Pick a threshold
        rval = random.random()
        if rval <= forward[strand][kmer]*correctForward:
          # Pick a length
          lval = random.random()
          selLength = lengthDist.bisect_left(lval)
        
          # Check the kmer at the other end
          rkmer = ""
          if strand: rkmer = context[posInWindow+selLength-reverseLength:posInWindow+selLength]
          else: rkmer = context[posInWindow-selLength+1:posInWindow-selLength+reverseLength+1]

          rval = random.random()
          if rval <= reverse[not strand][rkmer]*correctReverse:
            # Extract sequence and write BAM record
            selSeq = ""
            if strand: selSeq = context[posInWindow:posInWindow+selLength]
            else: selSeq = context[posInWindow-selLength+1:posInWindow+1]
            writeBAMentry(chrom,pos,selSeq,selLength,strand)
            #print strand,kmer,rkmer,selSeq,len(selSeq),selLength
      #else:
        #print rval,forward[strand][kmer]

    # ADVANCE IN WINDOW
    nbase = get_DNA_file(chrom,pos+maxLength+1,pos+maxLength+1,genome_map,genome_index).upper()
    if nbase != None:
      if len(context) == 2*maxLength+1:
        context = context[1:]+nbase
      else:
        context += nbase

    posInWindow = maxLength
    if (pos-maxLength) < 0:
      posInWindow = maxLength-(maxLength-pos)
    elif (pos+maxLength) > chromLength:
      posInWindow = max(0,maxLength-(chromLength-pos))
      
    count += 1
    #if count > 5: sys.exit()
    if options.verbose and (count % 100000 == 0): 
      sys.stderr.write(" ".join(map(str,["CurPos:",chrom,pos,context[posInWindow],posInWindow,get_DNA_file(chrom,pos+1,pos+1,genome_map,genome_index),"Sim:",rcounter]))+"\n")

if BAMoutfile != None:
  BAMoutfile.close()
