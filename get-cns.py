##===============================================================================================================
 
##  Usage:
##  python get-cns.py score-file CDS-file cutoff

##===============================================================================================================

#!/usr/bin/env python
#-*- coding:utf-8 -*-

from __future__ import division
import sys
import os

target=sys.argv[1]
q_name=os.path.basename(target)
name=os.path.splitext(q_name)[0]
cds_file=sys.argv[2]
cutoff=sys.argv[3]


def get_cns_fragments():
    f1=open(target,'r')
    f11=f1.readlines()
    f2=open('%s-larger'%(name),'w')

    for i in f11:
        score=i.split('\t')[3]
        if score >= cutoff:
            f2.write(i)

    f1.close()
    f2.close()
    os.system('bedtools subtract -a %s-larger -b %s>%s-cns-bases.bed'%(name,cds_file,name))
    os.system('bedtools merge -d 3 -i %s-cns-bases.bed>%s-merge3.bed'%(name,name))
    os.system('bedtools subtract -a %s-merge3.bed -b %s >%s-merge3-suubtract-cds.bed'%(name,cds_file,name))

def get_cns():
    f1=open('%s-merge3-suubtract-cds.bed'%(name),'r')
    f11=f1.readlines()
    f2=open('%s-cns.bed'%(name),'w')

    for i in f11:
        chro=i.split('\t')[0]
        start=int(i.split('\t')[1])
        end=int(i.split('\t')[2])
        if (end-start) >5:
            f2.write(i)

   # os.system('''awk '{if (($3-$2)>5);{print $0}}' %s-merge3-suubtract-cds.bed >%s-cns.bed'''%(name,name)) 
    os.system('rm %s-larger %s-cns-bases.bed %s-merge3.bed %s-merge3-suubtract-cds.bed'%(name,name,name,name))
    

if __name__=='__main__':
    get_cns_fragments()
    get_cns()
