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

    def __init__(self, full=False):
        self.doFullRecords = full
        
        # Read orientations
        # type 0 <--- --->
        # type 1 ---> --->
        # type 2 ---> <---
        self.OT = self.enum('OUT', 'SAME', 'IN', 'ERROR')
    
    def enum(*sequential, **named):
        enums = dict(zip(sequential, range(len(sequential))), **named)
        return type('Enum', (), enums)
    
    def OT2Str(self, OT):
        """For the humans!"""
        if OT == self.OT.OUT:
            return 'OUT'
        if OT == self.OT.SAME:
            return 'SAME'
        if OT == self.OT.IN:
            return 'IN'
        return 'ERROR'
    
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
                (OT,ins) = self.classifyBam(bam_file)
                if verbose:
                    print "Orientation: %s Insert: %d" % (self.OT2Str(OT), ins) 
                bam_types[bam_count] = (OT,ins)                
                bam_count += 1
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise
        return bam_types
    
    def classifyBam(self, bamFile, numPaired=1000):
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
            if(alignedRead.is_paired and            # basic QA so we're not wasting our time
               alignedRead.is_proper_pair and 
               not alignedRead.is_duplicate and 
               not alignedRead.mate_is_unmapped and 
               not alignedRead.is_unmapped
               ):
                # strip off any trailing ".1, _1, /1" which may be at the end of the read id
                query = re.sub("[_/\.].$", '', alignedRead.qname)
                if query in seen_reads:
                    # we have a pair!
                    # check to see they're on the same strand
                    # rname is deprecated in pysam >= 0.4 use tid instead!
                    if alignedRead.rname == seen_reads[query].rname:
                        if alignedRead.qstart < seen_reads[query].qstart:
                            # alignedRead is mapped BEFORE the saved one
                            (OT,ins) = self.determineOT(alignedRead, seen_reads[query])
                        else:
                            # otherwise it comes after
                            (OT,ins) = self.determineOT(seen_reads[query], alignedRead)
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
     
    def determineOT(self, ar1, ar2):
        """Determine the orientation type and insert size of two reads"""
        isize = abs(ar1.isize)
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################
