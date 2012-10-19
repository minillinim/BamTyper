#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bamtyper.py                                                              #
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
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
import os

# BamTyper imports
import utilities
 
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamTyperOptionsParser():
    def __init__(self): pass
    
    def parseOptions(self, options ):

        if(options.subparser_name == 'type'):
            # Work out the read orientation type and return
            BamParser = utilities.BamParser()
            BamParser.getTypes(options.bamfiles, verbose=True)

        elif(options.subparser_name == 'links'):
            # Calculate linking reads and possibly coverages
            BamParser = utilities.BamParser()
            if options.coverage:
                (filtered_links, ref_lengths, total_coverages) = BamParser.getLinks(options.bamfiles, full=options.verbose, doCoverage=True, minJoin=options.min)
                for cid in filtered_links:
                    print cid, 
                    for fl in filtered_links[cid]: 
                        print ", [",fl[0],",",fl[1],",",BamParser.LT2Str(fl[2], terse=True),", %0.0f ]" % fl[3],
                    print 
                for cid in ref_lengths:
                    line_vals = [cid]
                    for i in range(len(options.bamfiles)):
                        if cid in total_coverages[i]:
                            line_vals.append("%0.4f" %total_coverages[i][cid])
                        else:
                            line_vals.append("0.0")
                    print "\t".join(line_vals)
            else:
                filtered_links = BamParser.getLinks(options.bamfiles, full=options.verbose, doCoverage=False, minJoin=options.min)
                for cid in filtered_links:
                    print cid, 
                    for fl in filtered_links[cid]: 
                        print ", [",fl[0],",",fl[1],",",BamParser.LT2Str(fl[2], terse=True),", %0.0f ]" % fl[3],
                    print 
        return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################
