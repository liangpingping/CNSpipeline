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
        f3.write('fixedStep'+'\t'+'chrom='+chr_os+'\t'+'start='+str(start_os+1)+'\t'+'step=1'+'\n')
        for j in range(len(seq_os)):
            num=0
            if seq_os[j]!='-':
                start_os+=1
                for a in list[2:]:
                    a=','.join(a.split())
                    if a != '':
                        # name1=a.split(',')[1]
                        # beg=a.split(',')[2]
                        # num=a.split(',')[3]
                        seq1=a.split(',')[6]
                        if seq1[j]!='-':
                            num+=1
            else:
                continue
            score=num/int(num_species)
            f3.write(str(score)+'\n')
            f1.write(chr_os+'\t'+str(start_os-1)+'\t'+str(start_os)+'\t'+str(score)+'\t'+str(num)+'\n')

    f1.close()
    f2.close()
    f3.close()
    
if __name__=='__main__':
    main()
