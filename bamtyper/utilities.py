#!/usr/bin/env python
###############################################################################
#                                                                             #
#    utilities.py                                                            #
#                                                                             #
#    Wraps coarse workflows                                                   #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#  888888b.                      88888888888                                  #
#  888  "88b                         888                                      #
#  888  .88P                         888                                      #
#  8888888K.   8888b.  88888b.d88b.  888  888  888 88888b.   .d88b.  888d888  #
#  888  "Y88b     "88b 888 "888 "88b 888  888  888 888 "88b d8P  Y8b 888P"    #
#  888    888 .d888888 888  888  888 888  888  888 888  888 88888888 888      #
#  888   d88P 888  888 888  888  888 888  Y88b 888 888 d88P Y8b.     888      #
#  8888888P"  "Y888888 888  888  888 888   "Y88888 88888P"   "Y8888  888      #
#                                              888 888                        #
#                                         Y8b d88P 888                        #
#                                          "Y88P"  888                        #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.1.2"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
import os
import re
import numpy as np

import pysam

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamParser:
    """Parse multiple bam files and write link info to screen """

    def __init__(self):

        # Read orientations
        # type 0 <--- --->
        # type 1 ---> --->
        # type 2 ---> <---
        self.OT = enum('OUT', 'SAME', 'IN', 'ERROR')
        
        # Links are stored as triples (contig1, contig2, linktype)
        # There are 4 linktypes:
        # SS  <--1--- ---2--> 
        # SE  <--1--- <--2---
        # ES  ---1--> ---2-->
        # EE  ---1--> <--2---
        self.LT = enum('SS','SE','ES','EE','ERROR')
        
#------------------------------------------------------------------------------
# Managing orientation types
    
    def OT2Str(self, OT):
        """For the humans!"""
        if OT == self.OT.OUT:
            return 'OUT'
        if OT == self.OT.SAME:
            return 'SAME'
        if OT == self.OT.IN:
            return 'IN'
        return 'ERROR'

#------------------------------------------------------------------------------
# Managing linking types

    def LT2Str(self, link, terse=False):
        """For the humans!
        
        link is a link-tuple"""

        if terse:
            if link == self.LT.SS:
                return "SS"
            if link == self.LT.SE:
                return "SE"
            if link == self.LT.ES:
                return "ES"
            if link == self.LT.EE:
                return "EE"
            return "??"           
        if link[2] == self.LT.SS:
            return str(link[0]) + " lies before " + str(link[1]) + " in the opposite direction with gap "+str(link[3])
        if link[2] == self.LT.SE:
            return str(link[0]) + " lies after " + str(link[1]) + " in the same direction with gap "+str(link[3])
        if link[2] == self.LT.ES:
            return str(link[0]) + " lies before " + str(link[1]) + " in the same direction with gap "+str(link[3])
        if link[2] == self.LT.EE:
            return str(link[0]) + " lies after " + str(link[1]) + " in the opposite direction with gap "+str(link[3])
        return 'Who knows?'
    
    def invertLT(self, LT):
        """Reverse the contigs == invert the link type"""
        if LT == self.LT.SS:
            return self.LT.SS 
        if LT == self.LT.SE:
            return self.LT.ES 
        if LT == self.LT.ES:
            return self.LT.SE 
        if LT == self.LT.EE:
            return self.LT.EE 
        return self.LT.ERROR
        

#------------------------------------------------------------------------------
# Linking reads

    def getLinks(self, bamFiles, full=False, verbose=False, doCoverage=False, minJoin=0):
        """Get all linking reads"""
        # first we need to know the types of reads we're dealing with
        bam_types = self.getTypes(bamFiles)
        num_bams = len(bamFiles)
        #print bam_types
        bam_count = 0
        all_links = []
        total_coverages = {}
        ref_lengths = {}
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                if verbose:
                    print "    Parsing BAM '%s' (%d of %d)" % (getBamDescriptor(bf), (bam_count+1), num_bams)
                if doCoverage:
                    (links, coverage, lengths) = self.classifyBamLinks(bam_file, bam_types[bam_count], doCoverage=True)
                    all_links += links
                    total_coverages[bam_count] = coverage
                    for cid in lengths:
                        if cid not in ref_lengths:
                            ref_lengths[cid] = lengths[cid]
                else:
                    all_links += self.classifyBamLinks(bam_file, bam_types[bam_count], doCoverage=False)
                bam_count += 1
                bam_file.close()
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise

        filtered_links = self.filterLinks(all_links, full=full, minJoin=minJoin)
        if doCoverage:
            for bam_cov in total_coverages:
                for cid in total_coverages[bam_cov]:
                    total_coverages[bam_cov][cid] = float(total_coverages[bam_cov][cid]) / float(ref_lengths[cid])
            return (filtered_links, ref_lengths, total_coverages)
        return filtered_links

    def classifyBamLinks(self, bamFile, bamType, doCoverage=False):
        """Parse a bam file (handle) to extract linking reads
        
        numPaired refers to the number of mapped pairs we need before
        we are confident to make a decsion
        """
        # just a sanity check here
        if bamType[0] == self.OT.ERROR:
            return []
        all_links = []
        coverages = {}
        seen_reads = {}
        ref_lengths = [int(x) for x in bamFile.lengths]
        references = bamFile.references
        for alignedRead in bamFile.fetch():
            if doCoverage:
                ar_name = bamFile.getrname(alignedRead.tid)
                if ar_name not in coverages:
                    coverages[ar_name] = 1
                else:
                    coverages[ar_name] += 1
                    
            # strip off any trailing ".1, _1, /1" which may be at the end of the read id
            query = re.sub("[_/\.].$", '', alignedRead.qname)
            if query in seen_reads:
                # we have a pair!
                # check to see they're on the same strand
                # rname is deprecated in pysam >= 0.4 use tid instead!
                if alignedRead.tid != seen_reads[query].tid:
                    (gap, LT) = self.determineOTDifferentRefs(alignedRead, seen_reads[query], ref_lengths, bamType)
                    if LT != self.LT.ERROR:
                        #print self.LT2Str(self.makeLink(alignedRead.tid, seen_reads[query].tid, LT))
                        all_links.append(self.makeLink(bamFile.getrname(alignedRead.tid), bamFile.getrname(seen_reads[query].tid), LT, gap))
            else:
                seen_reads[query] = alignedRead
    
        if(doCoverage):
            return (all_links, coverages, dict(zip(references,ref_lengths)))
        else:
            return all_links
    
    def filterLinks(self, links, full=False, minJoin=20):
        """Take a list of link triples and filter it down to a set up unique links
        
        links is a list of link tuples (c1, c2, LT, gap)
        minJoin denotes how many links must be present to confirm
        """
        # first we must munge the info into something useful
        gaps_hash = {}
        for link in links:
            if full:
                print self.LT2Str(link)
            # make sure the data structure is there
            if link[0] not in gaps_hash:
                gaps_hash[link[0]] = {}
            if link[1] not in gaps_hash[link[0]]:
                gaps_hash[link[0]][link[1]] = {0:[],1:[],2:[],3:[],4:[]}
            gaps_hash[link[0]][link[1]][link[2]].append(link[3])
    
        # now go through and refine it further
        filtered_links = {}
        for c1 in gaps_hash:
            for c2 in gaps_hash[c1]:
                most_common_lt_count = 0
                most_common_lt = self.LT.ERROR
                for lt in gaps_hash[c1][c2]:
                    if len(gaps_hash[c1][c2][lt]) >= minJoin:
                        # we only want one type of link per contig pair
                        if len(gaps_hash[c1][c2][lt]) > most_common_lt_count:
                            most_common_lt_count = len(gaps_hash[c1][c2][lt])
                            most_common_lt = lt
                if most_common_lt_count != 0:
                    # we accept the link as real
                    if c1 not in filtered_links:
                        filtered_links[c1] = []
                    if c2 not in filtered_links:
                        filtered_links[c2] = []
                    # store linking contig + number of links
                    filtered_links[c1].append((c2, most_common_lt_count, most_common_lt,np.mean(gaps_hash[c1][c2][most_common_lt])))
                    filtered_links[c2].append((c1, most_common_lt_count, self.invertLT(most_common_lt),np.mean(gaps_hash[c1][c2][most_common_lt])))
        return filtered_links
    
    def determineOTDifferentRefs(self, ar1, ar2, lengths, bamType):
        """Determine the orientation type and insert size of two reads
        
        We assume that:
           1. both reads are on different references
           2. ar1 is the first read in the pair
           
        I swear this is the last time I write this code!
        """
        # first check to see if the start read lies in the right position
        max_ins = bamType[1] + (3 * bamType[2]) # mean + 3 * stdev
        adjusted_ins = bamType[1]
        rl = ar1.rlen
        if ar1.pos <= ( max_ins - 2 * rl ):
            # read 1 lies at the start of its contig
            r1_at_start = True
            adjusted_ins -= (ar1.pos + rl) 
            max_ins -= (ar1.pos + rl)
        elif ar1.pos >= (lengths[ar1.tid] - max_ins + rl):
            # read 1 lies at the end of its contig
            r1_at_start = False
            adjusted_ins -= (lengths[ar1.tid] - ar1.pos)
            max_ins -= (lengths[ar1.tid] - ar1.pos + 1)
        else:
            # read 1 lies in the middle of its contig
            return (0, self.LT.ERROR)
        # now check read 2
        if ar2.pos <= ( max_ins - 2 * rl ):
            # read 2 lies at the start of its contig
            r2_at_start = True
            adjusted_ins -= (ar2.pos + rl) 
        elif ar2.pos >= (lengths[ar2.tid] - max_ins + rl):
            # read 2 lies at the end of its contig            
            r2_at_start = False
            adjusted_ins -= (lengths[ar2.tid] - ar2.pos)
        else:
            # read 2 lies in the middle of its contig
            return (0, self.LT.ERROR)

        # now put it all together!
        # print r1_at_start, ar1.pos, ar1.is_reverse, "|", r2_at_start, ar2.pos, ar2.is_reverse, "|",     
        if r1_at_start:
            # -x-1-->
            if ar1.is_reverse:
                # -<-1-->
                if r2_at_start:
                    # <--1->- -x-2-->
                    if ar2.is_reverse:
                        # <--1->- -<-2--> ==> IN + SS
                        if bamType[0] == self.OT.IN:
                            # print "0 IN + SS"
                            return (adjusted_ins, self.LT.SS)
                    else:
                        # <--1->- ->-2--> ==> SAME + SS
                        if bamType[0] == self.OT.SAME:
                            # print "1 SAME + SS"
                            return (adjusted_ins, self.LT.SS)
                else:
                    # <--1->- <x-2---
                    if ar2.is_reverse:
                        # <--1->- <>-2--- ==> SAME + SE
                        if bamType[0] == self.OT.SAME:
                            # print "2 SAME + SE"
                            return (adjusted_ins, self.LT.SE)
                    else:
                        # <--1->- <<-2--- ==> IN + SE
                        if bamType[0] == self.OT.IN:
                            # print "3 IN + SE"
                            return (adjusted_ins, self.LT.SE)
            else: # r1 agrees
                # ->-1-->
                if r2_at_start:
                    # <--2-x- ->-1-->
                    if ar2.is_reverse:
                        # <--2->- ->-1--> ==> SAME + SS
                        if bamType[0] == self.OT.SAME:
                            # print "4 SAME + SS"
                            return (adjusted_ins, self.LT.SS)
                    else:
                        # <--2-<- ->-1--> ==> OUT + SS
                        if bamType[0] == self.OT.OUT:
                            # print "5 OUT + SS"
                            return (adjusted_ins, self.LT.SS)
                else:
                    # ---2-x> ->-1-->
                    if ar2.is_reverse:
                        # ---2-<> ->-1--> ==> OUT + SE
                        if bamType[0] == self.OT.OUT:
                            # print "6 OUT + SE"
                            return (adjusted_ins, self.LT.SE)
                    else:
                        # ---2->> ->-1--> ==> SAME + SE
                        if bamType[0] == self.OT.SAME:
                            # print "7 SAME + SE"
                            return (adjusted_ins, self.LT.SE)
        else: # r1 at end
            # ---1-x>
            if ar1.is_reverse:
                # ---1-<>
                if r2_at_start:
                    # ---1-<>- -x-2-->
                    if ar2.is_reverse:
                        # ---1-<> -<-2--> ==> SAME + ES
                        if bamType[0] == self.OT.SAME:
                            # print "8 SAME + ES"
                            return (adjusted_ins, self.LT.ES)
                    else:
                        # ---1-<> ->-2--> ==> OUT + ES
                        if bamType[0] == self.OT.OUT:
                            # print "9 OUT + ES"
                            return (adjusted_ins, self.LT.ES)
                else:
                    # ---1-<> <x-2---
                    if ar2.is_reverse:
                        # ---1-<> <>-2--- ==> OUT + EE
                        if bamType[0] == self.OT.OUT:
                            # print "a OUT + EE"
                            return (adjusted_ins, self.LT.EE)
                    else:
                        # ---1-<> <<-2--- ==> SAME + EE 
                        if bamType[0] == self.OT.SAME:
                            # print "b SAME + EE"
                            return (adjusted_ins, self.LT.EE)
            else: # r1 agrees
                # ---1->>
                if r2_at_start:
                    # ---1->> -x-2-->
                    if ar2.is_reverse:
                        # ---1->> -<-2--> ==> IN + ES
                        if bamType[0] == self.OT.IN:
                            # print "c IN + SS"
                            return (adjusted_ins, self.LT.ES)
                    else:
                        # ---1->> ->-2--> ==> SAME + ES
                        if bamType[0] == self.OT.SAME:
                            # print "d SAME + ES"
                            return (adjusted_ins, self.LT.ES)
                else:
                    # ---1->> <x-2---
                    if ar2.is_reverse:
                        # ---1->> <>-2--- ==> SAME + EE
                        if bamType[0] == self.OT.SAME:
                            # print "e SAME + EE"
                            return (adjusted_ins, self.LT.EE)
                    else:
                        # ---1->> <<-2--- ==> IN + EE
                        if bamType[0] == self.OT.IN:
                            # print "f IN + EE"
                            return (adjusted_ins, self.LT.EE)
        # print 
        return (0, self.LT.ERROR)
    
    def makeLink(self, rname1, rname2, LT, gap):
        """Consistent interface for naming links"""
        if rname1 < rname2:
            return (rname1, rname2, LT, gap)
        else:
            return (rname2, rname1, self.invertLT(LT), gap)
        
#------------------------------------------------------------------------------
# Working out read types

    def getTypes(self, bamFiles, verbose=False):
        """Parse multiple bam files, keep track of links"""
        # parse the BAMs
        tmp_storage = {}
        num_bams = len(bamFiles)

        bam_count = 0
        bam_types = {}
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                if verbose:
                    print "Determining OT for BAM '%s'" % (getBamDescriptor(bf))
                (OT,ins,std) = self.classifyBamType(bam_file)
                if verbose:
                    print "Orientation: %s Insert: %d, Stdev: %d" % (self.OT2Str(OT), ins, std) 
                bam_types[bam_count] = (OT,ins,std)                
                bam_count += 1
                bam_file.close()
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise
        return bam_types
    
    def classifyBamType(self, bamFile, numPaired=10000):
        """Parse a bam file (handle) to determine the read orientation type
        
        numPaired refers to the number of mapped pairs we need before
        we are confident to make a decsion
        """
        seen_reads = {}
        OTs = { self.OT.OUT : [], self.OT.SAME : [], self.OT.IN : [] } # store OT's vs insert lengths
        num_stored = 0
        for alignedRead in bamFile.fetch():
            if(num_stored > numPaired):
                break                               # we have enough
            # strip off any trailing ".1, _1, /1" which may be at the end of the read id
            query = re.sub("[_/\.].$", '', alignedRead.qname)
            if query in seen_reads:
                # we have a pair!
                # check to see they're on the same strand
                # rname is deprecated in pysam >= 0.4 use tid instead!
                if alignedRead.tid == seen_reads[query].tid:
                    if alignedRead.pos < seen_reads[query].pos:
                        # alignedRead is mapped BEFORE the saved one
                        (OT,ins) = self.determineOTSameRef(alignedRead, seen_reads[query])
                    else:
                        # otherwise it comes after
                        (OT,ins) = self.determineOTSameRef(seen_reads[query], alignedRead)
                    OTs[OT].append(ins)
                    num_stored += 1
            else:
                seen_reads[query] = alignedRead
        
        # now determine the stats for this bamfile
        # we'd like the vast majority of the OT's to agree
        cutoff = int(num_stored * 0.9)
        for OT in OTs:
            if len(OTs[OT]) > cutoff:
                return (OT, int(np.mean(OTs[OT])), int(np.std(OTs[OT])))
        
        # return garbage
        return (self.OT.ERROR, 0)
     
    def determineOTSameRef(self, ar1, ar2):
        """Determine the orientation type and insert size of two reads
        
        We assume that:
           1. both reads are on the same reference
           2. ar1 comes before ar2
        """
        isize = ar2.pos - ar1.pos + ar1.rlen
        if ar1.is_reverse:
            # <-1--
            if ar2.is_reverse:
                # <-1-- <-2--
                return (self.OT.SAME, isize)
            else:
                # <-1-- -2-->
                return (self.OT.OUT, isize)
        else:
            # -->
            if ar2.is_reverse:
                # --1-> <-2--
                return (self.OT.IN, isize)
            else:
                # --1-> --2->
                return (self.OT.SAME, isize)
    
def getBamDescriptor(fullPath):
    """AUX: Reduce a full path to just the file name minus extension"""
    return os.path.splitext(os.path.basename(fullPath))[0]

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
