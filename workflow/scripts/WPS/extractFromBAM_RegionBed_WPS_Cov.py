#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *03.06.2014
"""

import sys, os
import argparse
import gzip
import pysam
import random
import math

from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def regionFileParser(infile):
  for line in infile:
    ########
    # implement proper bedfile reading
    ########
    chrom,start,end,cid,score,strand = line.split() # positions should be 0-based and end non-inclusive
    if chrom.startswith("chr"):
      chrom = chrom.replace("chr","")
    yield chrom, start, end, cid, score, strand
  return 

def parseRegion(regionstr):
  try:
    chrom=regionstr.split(":")[0]
    start,end=(regionstr.split(":")[-1]).split("-")
    start=str(int(start)-1)
    if chrom.startswith("chr"):
      chrom = chrom.replace("chr","")
    yield chrom, start, end, regionstr, "0", "+"
    return
  except:
    sys.stderr.write("Failed parsing coordinate string %s\n"%(regionstr))
  return


parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='BAM files with aligned reads.')
parser.add_argument("-i","--input", dest="input", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 76)",default=76,type=int)
parser.add_argument("-m","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_argument("-t","--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
parser.add_argument("-w","--protection", dest="protection", help="Base pair protection window size (default 120)",default=120,type=int)
#parser.add_argument("-c","--method", dest="method", help="Type of protection score to calculate. The default is: %(default)s and your choices are: %(choices)s.",default="WPSv1",type=str, choices=("WPSv1","WPSv2","WPSv3","WPSv4","WPSv5"))
parser.add_argument("-o","--outfile", dest="outfile", help="Outfile prefix (def 'block_%%s.tsv.gz')",default='block_%s.tsv.gz') # reserve atleast 6 digits 
parser.add_argument("-e","--empty", dest="empty", help="Keep files of empty blocks (def Off)",default=False,action="store_true")
parser.add_argument("--minInsert", dest="minInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type=int)
parser.add_argument("--maxInsert", dest="maxInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type=int)
parser.add_argument("--max_length", dest="max_length", help="Assumed maximum insert size (default 1000)",default=1000,type=int)
parser.add_argument("--downsample", dest="downsample", help="Ratio to down sample reads (default OFF)",default=None,type=float)
parser.add_argument("--onefile", dest="onefile", help="Print as single output to stdout (default OFF)",default=False,action="store_true")
parser.add_argument("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
options = parser.parse_args()

minInsSize,maxInsSize = None,None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
  minInsSize = options.minInsSize
  maxInsSize = options.maxInsSize
  if options.verbose: sys.stderr.write("Using min/max length cutoffs: %d/%d\n"%(minInsSize,maxInsSize))
  
options.outfile = options.outfile.strip("""\'""")
outfiles = {}
if not options.onefile:
  outfiles = { 'WPS':gzip.open(options.outfile%"WPS",'wt'), 'COV':gzip.open(options.outfile%"COV",'wt'), 'STARTS':gzip.open(options.outfile%"STARTS",'wt') }

protection = options.protection//2

#validChroms = set(map(str,range(1,23)+["X","Y"]))
#validChroms = set(map(str,range(1,23)+["X","Y"]))
#validChroms = [str(i) for i in range(1, 23)] + ["X","Y"]

if os.path.exists(options.input):
  if ".gz" in options.input:
    infile = gzip.open(options.input,"r")
  else:
    infile = open(options.input)
  regionIterator = regionFileParser(infile)
else:
  regionIterator = parseRegion(options.input)

for chrom,start,end,cid,score,strand in regionIterator:
    #if chrom not in validChroms: continue
    #regionStart,regionEnd = int(start)-300,int(end)+300
    regionStart,regionEnd = int(start),int(end)
    #if regionStart < 1: continue
    
    #outchrom = chrom.replace("chr","")
    #if outchrom.startswith('gi|'):
      #NCfield = outchrom.split("|")[-2]
      #if NCfield.startswith("NC"):
        #outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
    if options.verbose:
       sys.stderr.write("Processing region: %s:%s-%s %s %s %s\n"%(chrom,start,end,cid,score,strand))
    
    posRange = defaultdict(lambda:[0,0])
    filteredReads = Intersecter()

    for bamfile in options.files:
      if options.verbose: sys.stderr.write("Reading %s\n"%bamfile)
      bamfile = bamfile.strip("""\'""")
      if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
        input_file = pysam.Samfile( bamfile, "rb" )
        prefix = ""
        for tchrom in input_file.references:
          if tchrom.startswith("chr"):
            prefix = "chr"
            break
        if options.verbose: sys.stderr.write("Retrieving reads...\n")
        for read in input_file.fetch(prefix+chrom,max(regionStart-protection-1,0),regionEnd+protection+1):
          if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
          if isSoftClipped(read.cigar): continue
          
          if read.is_paired:
            if read.mate_is_unmapped: continue
            if read.rnext != read.tid: continue
            if read.is_read1 or (read.is_read2 and read.pnext+read.qlen < regionStart-protection-1):
              if read.isize == 0: continue
              if options.downsample != None and random.random() >= options.downsample: continue
              rstart = min(read.pos,read.pnext)+1 # 1-based
              lseq = abs(read.isize)
              rend = rstart+lseq-1 # end included
              if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)): continue

              #if options.verbose: sys.stderr.write("Adding interval for read...\n")
              filteredReads.add_interval(Interval(rstart,rend))
              #print read.qname,rstart,rend,rend-rstart,abs(read.isize)
              for i in range(rstart,rend+1):
                if i >= regionStart and i <= regionEnd:
                  posRange[i][0]+=1
              if rstart >= regionStart and rstart <= regionEnd:
                posRange[rstart][1]+=1
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
          else:
            if options.downsample != None and random.random() >= options.downsample: continue
            rstart = read.pos+1 # 1-based
            lseq = aln_length(read.cigar)
            rend = rstart+lseq-1 # end included
            if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)): continue
            
            filteredReads.add_interval(Interval(rstart,rend))
            #print read.qname,rstart,rend,rend-rstart,aln_length(read.cigar)
            for i in range(rstart,rend+1):
              if i >= regionStart and i <= regionEnd:
                posRange[i][0]+=1
            if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.lengthSR-10)):
              if (rstart >= regionStart and rstart <= regionEnd):
                posRange[rstart][1]+=1
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
            elif read.is_reverse:
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
            else:
              if (rstart >= regionStart and rstart <= regionEnd):
                posRange[rstart][1]+=1
      else:
        sys.stderr.write("File without BAM index skipping: %s\n"%(bamfile))
        if options.verbose: sys.stderr.write("Warning: without input reads, output files will be populated with zero/empty values.\n")

    if options.verbose: sys.stderr.write("Evaluating posRange vector...\n")
    #filename = options.outfile%cid
    #outfile = gzip.open(filename,'w')
    cov_sites = 0
    outLines = []
    wps_list = []
    cov_list = []
    starts_list = []

    for pos in range(regionStart,regionEnd+1):
      rstart,rend = pos-protection,pos+protection
      gcount,bcount,ecount = 0.0,0.0,0.0
      for read in filteredReads.find(rstart,rend):
        if (read.start > rstart):
          bcount += 1.0
        elif (read.end < rend): 
          ecount += 1.0
        else: 
          gcount += 1.0
      covCount,startCount = posRange[pos]
      cov_sites += covCount
      wpsValue = gcount-(bcount+ecount)
      #if (options.method != "WPSv1") and (2*gcount+bcount+ecount > 1):
      #  if options.method == "WPSv2":
      #    wpsValue = 2.0*gcount/(2.0*gcount+bcount+ecount)
      #  elif options.method == "WPSv3":
      #    wpsValue = 2.0*gcount/(2.0*gcount+bcount+ecount)
      #    wpsValue = wpsValue * (1.0 - math.sqrt( (wpsValue*(1.0-wpsValue)) / (2.0*gcount+bcount+ecount-1.0) ))
      #  elif options.method == "WPSv4":
      #    wpsValue = max(2.0*gcount-abs(bcount-ecount),0)/(2.0*gcount+bcount+ecount)
      #  elif options.method == "WPSv5":
      #    wpsValue = max(2.0*gcount-abs(bcount-ecount),0)/(2.0*gcount+bcount+ecount)
      #    wpsValue = wpsValue * (1.0 - math.sqrt( (wpsValue*(1.0-wpsValue)) / (2.0*gcount+bcount+ecount-1.0) ))
      #elif (options.method != "WPSv1"): 
      #  wpsValue = 0.0
      if options.onefile:
         outLines.append("%s\t%d\t%.4f\t%.4f\t%.4f\n"%(chrom,pos,covCount,startCount,wpsValue))
      else: 
         wps_list.append(wpsValue)
         cov_list.append(covCount)
         starts_list.append(startCount)

    if options.onefile:
       if strand == "-": outLines = outLines[::-1]
       for line in outLines: sys.stdout.write(line)
    else:
      if strand == "-": wps_list = wps_list[::-1]
      if strand == "-": cov_list = cov_list[::-1]
      if strand == "-": starts_list = starts_list[::-1]
  
      outfiles['WPS'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),wps_list))+"\n")
      outfiles['COV'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),cov_list))+"\n")
      outfiles['STARTS'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),starts_list))+"\n")

if not options.onefile:
  for name,filestream in outfiles.items():
    filestream.close()
