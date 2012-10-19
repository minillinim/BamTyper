                                                                               
    888888b.                      88888888888                                  
    888  "88b                         888                                      
    888  .88P                         888                                      
    8888888K.   8888b.  88888b.d88b.  888  888  888 88888b.   .d88b.  888d888  
    888  "Y88b     "88b 888 "888 "88b 888  888  888 888 "88b d8P  Y8b 888P"    
    888    888 .d888888 888  888  888 888  888  888 888  888 88888888 888      
    888   d88P 888  888 888  888  888 888  Y88b 888 888 d88P Y8b.     888      
    8888888P"  "Y888888 888  888  888 888   "Y88888 88888P"   "Y8888  888      
                                                888 888                        
                                           Y8b d88P 888                        
                                            "Y88P"  888                         
                                                                             

# Overview

Ability to work out the orientation and insert size of a paired read data file
Can estimate relative orientation and gap between pairs of contigs in the bam file (Useful for scaffolding) 

# Installation

use pip:

```sh
$ pip install BamTyper
```
# Usage
##As a library:
###Get the type of the reads:
    #!/usr/bin/env python
    from bamtyper.utilities import BamParser as BTBP
    BP = BTBP()
    bam_types = BP.getTypes(bamFiles)

    Where:
        bamFiles - a list of BAM filenames
    
    Returns:
        bam_types - a dict containing information about the insert size and 
                    relative orientation of reads in the bam file
                      { bam1 : (type, ins, stdev), ... }
                    Where:
                        type - orientation type:
                            0 : OUT  : <--- --->
                            1 : SAME : ---> --->
                            2 : IN   : ---> <---
                        ins - estimated insert size (of original DNA fragment)
                        stdev - standard deviation of insert size      

###Get linking pairs:
    #!/usr/bin/env python
    from bamtyper.utilities import BamParser as BTBP
    BP = BTBP()
    (links, ref_lengths, coverages) = BP.getLinks(bamFiles, doCoverage=True)
    
    Where:
        bamFiles - a list of BAM filenames
    
    Returns:
        links - a dict containing information about links between two contigs
                  {c1: (c2, num_links, link_type, gap), ... }
                Where:
                    num_links - Number of paired reads confirming the link
                    link_type - Relative orientation of the two contigs (Start and End)
                      SS, SE, ES, EE or ERROR
                    gap - Estimated gap between te two contigs
                    
        coverages - a dict containing the FRAGMENT coverage of each contig n the bam file(s)
                  {c1 : (cov1, cov2, ..., covN), ... }
                  
        ref_lengths - a dict containing the lengths of all contigs
                  {c1 : len, ... }
    
    Notes:
        bamtyper will automatically work out the orientation and insert size of the 
        reads in each bam file and base it's estimations of link_type and gap on this 

##On the command line:

    bamtyper type - Parse BAM files and determine reads type
    
        $ bamtyper type bamfile.bam
        
        Determining OT for BAM 'bamfile'
        Orientation: IN Insert: 301, Stdev: 29
 
    bamtyper links - Parse BAM files and get linking reads

    Usage 1:

        $ bamtyper links bamfile.bam 

      1.
        contig2 , [ contig1 , 39 , SE , 69 ]
        contig1 , [ contig2 , 39 , ES , 69 ]
        
        implies a layout which looks like:
        
        ---1--> 69bp ---2-->

      2.
        contig3 , [ contig2 , 3 , SS , 58 ]
        contig2 , [ contig3 , 3 , SS , 58 ] , [ contig1 , 4 , EE , 45 ]
        contig1 , [ contig2 , 4 , EE , 45 ]

        implies a layout which looks like:
        
        ---1--> 45bp <--2--- 58bp ---3-->

    Usage 2: report FRAGMENT coverage too!

        $ bamtyper links bamfile1.bam bamfile2.bam -c

        contig3 , [ contig2 , 3 , SS , 58 ]
        contig2 , [ contig3 , 3 , SS , 58 ] , [ contig1 , 4 , EE , 45 ]
        contig1 , [ contig2 , 4 , EE , 45 ]
        contig3 0.6206 0.5234
        contig2 0.6558 0.0123
        contig1 0.6523 0.5634

        Where:
            contig3 0.6206 0.5234
            
            Reports 0.6206 fragments per base in bamfile1 and 0.5234 in bamfile2
            If these were 100bp reads then this would imply
            coverages of 62x and 52x respectively
            
# Administration

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/BamTyper

This software is currently unpublished.

Copyright © 2012 Michael Imelfort. See LICENSE.txt for further details.

