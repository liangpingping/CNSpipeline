#!/usr/bin/env python
#-*- coding:utf-8 -*-

##------------------------------------------------------------------------------------------------
##
## Usages:
##
##        python *.py query_genome.fa number_of_thread
##
##------------------------------------------------------------------------------------------------

from __future__ import division
import os
from Bio import SeqIO
import sys
import time
import subprocess
import math

query=sys.argv[1]	
q_abspath=os.path.abspath(query)
q_fullname=os.path.basename(query)
q_name=os.path.splitext(q_fullname)[0]
threads=sys.argv[2]	
	
os.system('mkdir -p ./%s/chain ./%s/maf  ./%s/Nib_q '%(q_name,q_name,q_name))
	


if __name__=='__main__': 
    main()
