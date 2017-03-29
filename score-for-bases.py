##===============================================================================================================
##     This script used to counting the score of every base from maf files,  we regard the bases of query species
## which mismatch to the reference species as match spot. That is,once the sites in query species are not bases(-),
## then we think it is aligned to the reference species.
 
##  Usage:
##  python script.py *.maf num_species 
##  It will generates the results named *-score.bed,*-score.wig
##===============================================================================================================

#!/usr/bin/env python
#-*- coding:utf-8 -*-

from __future__ import division
import sys
import os

target=sys.argv[1]
q_name=os.path.basename(target)
name=os.path.splitext(q_name)[0]
num_species=sys.argv[2]

os.system('maf-sort %s >sort-%s'%(target,name))
os.system("sed -i '1d' sort-%s "%(name))
os.system("sed -i '$d' sort-%s "%(name))
os.system("sed -i '/^#/d' sort-%s "%(name))

def main():
    f1=open('%s-score.bed'%(name),'w')
    f2=open('sort-%s'%(name),'r')
    f3=open('%s-score.wig'%(name),'w')
    f22=f2.read().split('\n\n')

    for i in range(len(f22)):
        list=f22[i].split('\n')
        list1=','.join(list[1].split())
        chr_os=list1.split(',')[1]
        start_os=int(list1.split(',')[2])
        match=list1.split(',')[3]
        end_os=str(int(start_os)+int(match))
        seq_os=list1.split(',')[6]
  
    f1.close()
    f2.close()
    f3.close()
    
if __name__=='__main__':
    main()
