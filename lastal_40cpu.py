#!/usr/bin/env python
#-*- coding:utf-8 -*-
##------------------------------------------------------------------------------------------------
##
## Usages:
##
##        python *.py query_genome.fa
##
## NB: this script uses 40cpu.
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
	
os.system('mkdir -p ./%s/chain ./%s/maf  ./%s/Nib_q '%(q_name,q_name,q_name))
	
def deal_query():
    tf3=open(query,'r')
    tf4=open('%s/q1.fa'%(q_name),'w')
    tf44=open('./reference/%s-multi'%(q_name),'w')
    for seq_record in SeqIO.parse(tf3,'fasta'):
	id=str(seq_record.id)
        seq=str(seq_record.seq)  
        length=str(len(seq))
        tf4.write('>'+ q_name +'.'+ id +'\n'+ seq +'\n')
        tf44.write('>'+ q_name + ':'+ id + ':' + '1' + ':' + '+' + ':' + length + '\n'+ seq +'\n')
    tf3.close()
    tf4.close()
    tf44.close()
    os.system('fold -w 70 ./%s/q1.fa >./%s/q.fa'%(q_name,q_name))
    os.system('fold -w 70 ./reference/%s-multi> ./reference/%s'%(q_name,q_name))
    os.system('rm ./%s/q1.fa ./reference/%s-multi'%(q_name,q_name))
	
def faSplit():	
    os.system('faSplit byname ./%s/q.fa ./%s/Nib_q/'%(q_name,q_name))
    os.system('faToTwoBit ./%s/q.fa ./%s/q.2bit'%(q_name,q_name))
    os.system('faSize ./%s/q.fa -detailed > ./%s/q.sizes'%(q_name,q_name))

def faToNib():
    tf5=open('%s/fatonib.sh'%(q_name),'w')
    tf55=os.listdir('./%s/Nib_q/'%(q_name))
    for i in tf55:
        faname=os.path.splitext(i)[0]
        tf5.write('faToNib %s/Nib_q/%s %s/Nib_q/%s.nib \n'%(q_name,i,q_name,faname))
	
def lastal():
    list2=open('%s/list_of_lastal.sh'%(q_name),'w')
    files2=os.listdir('./%s/Nib_q/'%(q_name))
    length=len(files2)
    lines=math.ceil(length/40)
    if lines <2:
        lines=2
    elif int(lines)%2 ==0:
        lines=int(lines)
    else:
        lines=int(lines)+1
    for name2 in files2:
        ext=os.path.splitext(name2)[1]
        if ext == '.fa':
            exname2=os.path.splitext(name2)[0]
            list2.write('lastal -p HOXD70 -e4000 -C2 -P2 -m100 ./reference/tdb  ./%s/Nib_q/%s > %s/maf/%s.maf \n\
            last-split -m1 %s/maf/%s.maf | maf-convert psl | axtChain stdin ./reference/Nib_t/ ./%s/Nib_q/ \
            ./%s/chain/%s.chain -linearGap=loose -psl  \n' %(q_name,name2,q_name,exname2,q_name,exname2,q_name,q_name,exname2))
    list2.close()
    os.system('cp all_lastal_40.sh all_sub_40.sh others %s/'%(q_name))
    os.system('sed -i "s/name/%s/g" %s/all_lastal_40.sh'%(q_name,q_name))
    os.system('sed -i "s/name/%s/g" %s/all_sub_40.sh'%(q_name,q_name))
    os.system('split -l %s ./%s/list_of_lastal.sh -d -a 2 ./%s/lastal_'%(lines,q_name,q_name))
    os.system('split -l 10 %s/all_lastal_40.sh -d -a 2 ./%s/sub'%(q_name,q_name))
    no_cpu=int(math.ceil(length/lines))
    if no_cpu <40:
        os.system('head -n %s %s/all_sub_40.sh >%s/sub.sh'%(no_cpu,q_name,q_name))
        os.system('parallel -j 40 < %s/sub.sh'%(q_name))
    elif os.path.isfile('%s/lastal_40'%(q_name)):
        os.system('cat ./%s/lastal_40 >> %s/lastal_39'%(q_name,q_name))
        os.system('parallel -j 40 < %s/all_sub_40.sh'%(q_name))
    else:
        os.system('parallel -j 40 < %s/all_sub_40.sh'%(q_name))

def others():    
    time.sleep(60)            
    os.system('sed -i "s/q_name/%s/g" %s/others'%(q_name,q_name))
    ID=subprocess.Popen('cat %s/jobID.txt | tr "\n" ":" |tr -d ".melon" | sed s/.$//'%(q_name),shell=True,stdout=subprocess.PIPE)
    allID=ID.communicate()[0]
    os.system('sed -i "s/all_jobID/%s/g" %s/others'%(allID,q_name))
    os.system('qsub ./%s/others'%(q_name))

def main():
    deal_query()
    faSplit()
    faToNib()
    os.system('parallel -k -j 20 <%s/fatonib.sh'%(q_name))
    lastal()
    others()

if __name__=='__main__': 
    main()
