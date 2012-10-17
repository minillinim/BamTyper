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
#                                         "Y88P"  888                         #
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
import os
import re
import numpy

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
        # STARTSEND   <--1--- <--2---
        # ES    ---1--> ---2-->
        # EE      ---1--> <--2---
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

    def LT2Str(self, link):
        """For the humans!
        
        link is a link-triple"""
        if link[2] == self.LT.SS:
            return "contig %d lies before contig %d in the opposite direction" % (link[0], link[1])
        if link[2] == self.LT.SE:
            return "contig %d lies after contig %d in the same direction" % (link[0], link[1])
        if link[2] == self.LT.ES:
            return "contig %d lies before contig %d in the same direction" % (link[0], link[1])
        if link[2] == self.LT.EE:
            return "contig %d lies after contig %d in the opposite direction" % (link[0], link[1])
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

    def getLinks(self, bamFiles, full=False):
        """Get all linking reads"""
        # first we need to know the types of reads we're dealing with
        bam_types = self.getTypes(bamFiles)
        bam_count = 0
        all_links = []
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                all_links += self.classifyBamLinks(bam_file, bam_types[bam_count])
                bam_count += 1
                bam_file.close()
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise

        print len(all_links)
        self.filterLinks(all_links)

    def classifyBamLinks(self, bamFile, bamType):
        """Parse a bam file (handle) to extract linking reads
        
        numPaired refers to the number of mapped pairs we need before
        we are confident to make a decsion
        """
        # just a sanity check here
        if bamType[0] == self.OT.ERROR:
            return []
        all_links = []
        seen_reads = {}
        ref_lengths = bamFile.lengths
        for alignedRead in bamFile.fetch():
            # strip off any trailing ".1, _1, /1" which may be at the end of the read id
            query = re.sub("[_/\.].$", '', alignedRead.qname)
            if query in seen_reads:
                # we have a pair!
                # check to see they're on the same strand
                # rname is deprecated in pysam >= 0.4 use tid instead!
                if alignedRead.tid != seen_reads[query].tid:
                    fix here
                    if alignedRead.is_read1:
                        print "K", seen_reads[query].is_read1 
                        LT = self.determineOTDifferentRefs(alignedRead, seen_reads[query], ref_lengths, bamType)
                        if LT != self.LT.ERROR:
                            #print self.LT2Str(self.makeLink(alignedRead.tid, seen_reads[query].tid, LT))
                            all_links.append(self.makeLink(alignedRead.tid, seen_reads[query].tid, LT))
                    else:
                        print "J", seen_reads[query].is_read1 
                        LT = self.determineOTDifferentRefs(seen_reads[query], alignedRead, ref_lengths, bamType)
                        if LT != self.LT.ERROR:
                            #print self.LT2Str(self.makeLink(seen_reads[query].tid, alignedRead.tid, LT))
                            all_links.append(self.makeLink(seen_reads[query].tid, alignedRead.tid, LT))
            else:
                seen_reads[query] = alignedRead
    
        return all_links
    
    def filterLinks(self, links, minJoin=20):
        """Take a list of link triples and filter it down to a set up unique links
        
        links is a list of link triples (c1, c2, LT)
        minJoin denotes how many links must be present to confirm
        """
        # first we must munge the info into something useful
        links_hash = {}
        for link in links:
            # make sure the data structure is there
            if link[0] not in links_hash:
                links_hash[link[0]] = {}
            if link[1] not in links_hash[link[0]]:
                links_hash[link[0]][link[1]] = {0:0,1:0,2:0,3:0,4:0} # 'SS':0,'SE':0,'ES':0,'EE':0,'ERROR':0
            
            links_hash[link[0]][link[1]][link[2]] += 1
            print self.LT2Str(link), link
    
        # now go through and refine it further
        filtered_links = {}
        for c1 in links_hash:
            for c2 in links_hash[c1]:
                for lt in links_hash[c1][c2]:
                    if links_hash[c1][c2][lt] >= minJoin:
                        # we accept the link as real
                        if c1 not in filtered_links:
                            filtered_links[c1] = []
                        if c2 not in filtered_links:
                            filtered_links[c2] = []
                        filtered_links[c1].append(c2)
                        filtered_links[c2].append(c1)
    
    def determineOTDifferentRefs(self, ar1, ar2, lengths, bamType):
        """Determine the orientation type and insert size of two reads
        
        We assume that:
           1. both reads are on different references
           2. ar1 is the first read in the pair
           
        I swear this is the last time I write this code!
        """
        # first check to see if the start read lies in the right position
        adjusted_ins = 0
        if ar1.pos <= bamType[1]:
            # read 1 lies at the start of its contig
            r1_at_start = True
            adjusted_ins = bamType[1] - r1_at_start 
        elif ar1.pos >= lengths[ar1.tid] - bamType[1]:
            # read 1 lies at the end of its contig
            r1_at_start = False
            adjusted_ins = bamType[1] - (lengths[ar1.tid] - ar1.pos)
        else:
            # read 1 lies in the middle of its contig
             return self.LT.ERROR
        # now check read 2
        if ar2.pos <= adjusted_ins:
            # read 2 lies at the start of its contig
            r2_at_start = True
        elif ar2.pos >= lengths[ar2.tid] - adjusted_ins:
            # read 2 lies at the end of its contig            
            r2_at_start = False
        else:
            # read 2 lies in the middle of its contig
            return self.LT.ERROR

        # no put it all together!        
        if r1_at_start:
            # -x-1-->
            if ar1.is_reverse:
                # -<-1-->
                if r2_at_start:
                    # <--1->- -x-2-->
                    if ar2.is_reverse:
                        # <--1->- -<-2--> ==> IN + SS
                        if bamType[0] == self.OT.IN:
                            return self.LT.SS
                    else:
                        # <--1->- ->-2--> ==> SAME + SS
                        if bamType[0] == self.OT.SAME:
                            return self.LT.SS
                else:
                    # <--1->- <x-2---
                    if ar2.is_reverse:
                        # <--1->- <>-2--- ==> SAME + SE
                        if bamType[0] == self.OT.SAME:
                            return self.LT.SE
                    else:
                        # <--1->- <<-2--- ==> IN + SE
                        if bamType[0] == self.OT.IN:
                            return self.LT.SE
            else: # r1 agrees
                # ->-1-->
                if r2_at_start:
                    # <--2-x- ->-1-->
                    if ar2.is_reverse:
                        # <--2->- ->-1--> ==> SAME + SS
                        if bamType[0] == self.OT.SAME:
                            return self.LT.SS
                    else:
                        # <--2-<- ->-1--> ==> OUT + SS
                        if bamType[0] == self.OT.OUT:
#                            print ar1.tid, ar1.pos, ar1.is_unmapped 
#                            print ar2.tid, ar2.pos, ar2.is_unmapped
#                            print "=======4"
                            return self.LT.SS
                else:
                    # ---2-x> ->-1-->
                    if ar2.is_reverse:
                        # ---2-<> ->-1--> ==> OUT + SE
                        if bamType[0] == self.OT.OUT:
                            return self.LT.SE
                    else:
                        # ---2->> ->-1--> ==> SAME + SE
                        if bamType[0] == self.OT.SAME:
                            return self.LT.SE
        else: # r1 at end
            # ---1-x>
            if ar1.is_reverse:
                # ---1-<>
                if r2_at_start:
                    # ---1-<>- -x-2-->
                    if ar2.is_reverse:
                        # ---1-<> -<-2--> ==> SAME + ES
                        if bamType[0] == self.OT.SAME:
                            return self.LT.ES
                    else:
                        # ---1-<> ->-2--> ==> OUT + ES
                        if bamType[0] == self.OT.OUT:
                            return self.LT.ES
                else:
                    # ---1-<> <x-2---
                    if ar2.is_reverse:
                        # ---1-<> <>-2--- ==> OUT + EE
                        if bamType[0] == self.OT.OUT:
                            return self.LT.EE
                    else:
                        # ---1-<> <<-2--- ==> SAME + EE 
                        if bamType[0] == self.OT.SAME:
                            return self.LT.EE
            else: # r1 agrees
                # ---1->>
                if r2_at_start:
                    # ---1->> -x-2-->
                    if ar2.is_reverse:
                        # ---1->> -<-2--> ==> IN + ES
                        if bamType[0] == self.OT.IN:
                            return self.LT.ES
                    else:
                        # ---1->> ->-2--> ==> SAME + ES
                        if bamType[0] == self.OT.SAME:
                            return self.LT.ES
                else:
                    # ---1->> <x-2---
                    if ar2.is_reverse:
                        # ---1->> <>-2--- ==> SAME + EE
                        if bamType[0] == self.OT.SAME:
                            return self.LT.EE
                    else:
                        # ---1->> <<-2--- ==> IN + EE
                        if bamType[0] == self.OT.IN:
                            return self.LT.EE
        return self.LT.ERROR
    
    def makeLink(self, rname1, rname2, LT):
        """Consistent interface for naming links"""
        if rname1 < rname2:
            return (rname1, rname2, LT)
        else:
            return (rname2, rname1, self.invertLT(LT))
        
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
                (OT,ins) = self.classifyBamType(bam_file)
                if verbose:
                    print "Orientation: %s Insert: %d" % (self.OT2Str(OT), ins) 
                bam_types[bam_count] = (OT,ins)                
                bam_count += 1
                bam_file.close()
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise
        return bam_types
    
    def classifyBamType(self, bamFile, numPaired=1000):
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
                return (OT, numpy.median(OTs[OT]))
        
        # return garbage
        return (self.OT.ERROR, 0)
     
    def determineOTSameRef(self, ar1, ar2):
        """Determine the orientation type and insert size of two reads
        
        We assume that:
           1. both reads are on the same reference
           2. ar1 comes before ar2
        """
        isize = ar2.pos - ar1.pos
        if ar1.is_reverse:
            # <--
            if ar2.is_reverse:
                # <-- <--
                return (self.OT.SAME, isize)
            else:
                # <-- -->
                return (self.OT.OUT, isize)
        else:
            # -->
            if ar2.is_reverse:
                # --> <--
                return (self.OT.IN, isize)
            else:
                # --> -->
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
