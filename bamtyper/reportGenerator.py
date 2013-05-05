#!/usr/bin/env python
###############################################################################
#                                                                             #
#    reportGenerator.py                                                       #
#                                                                             #
#    Contig / mapping reports                                                 #
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
from sys import exc_info, exit
import pysam
import errno
import gzip
from os.path import splitext as op_splitext, basename as op_basename
from string import maketrans as s_maketrans
from os.path import splitext as osp_splitext, basename as osp_basename, dirname as osp_dirname, join as osp_join
from os import makedirs as os_makedirs
from re import sub as re_sub

import numpy as np
np.seterr(all='raise')     

# BamTyper imports
from utilities import BamParser

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ReportGenerator:
    """main wrapper for producing reports on assembly mappings"""
    def __init__(self, contigsFileName, bamFileNames):
        self.bamFileNames = bamFileNames
        self.contigsFile = contigsFileName
        
        print contigsFileName, bamFileNames
 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self):
        self.compl = s_maketrans('ACGTacgtnN', '0110011000')

    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break
                
    def parse(self, contigFile, kse=None):
        """Do the heavy lifting of parsing"""
        if kse is None:
            kse = KmerSigEngine()
        tmp_storage = {} # save everything here first so we can sort accordingly
        for cid,seq,qual in self.readfq(contigFile):
            tmp_storage[cid] = (kse.getKSig(seq.upper()), len(seq))

        # get a list of contig names
        # This array is used like everywhere dude...
        con_names = {}
        for cid in sorted(tmp_storage.keys()):     
            con_names[cid] = tmp_storage[cid][1]
        return (tmp_storage, con_names)

    def getGC(self, seq):
        """Get the GC of a sequence"""
        Ns = seq.count('N') + seq.count('n')
        return sum([float(x) for x in list(seq.translate(self.compl))])/float(len(seq) - Ns)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
