#!/usr/bin/env python
#-*- coding:utf-8 -*-

##------------------------------------------------------------------------------------------------
##
## Usages:
##
##        python *.py query_genome.fa number_of_thread reference_genome_name
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
reference_genome_name=sys.argv[3]	
	
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
    list3=open('%s/list_of_chain.sh'%(q_name),'w')
    files2=os.listdir('./%s/Nib_q/'%(q_name))
    length=len(files2)
    for name2 in files2:
        ext=os.path.splitext(name2)[1]
        if ext == '.fa':
            exname2=os.path.splitext(name2)[0]
            list2.write('lastal -p HOXD70 -e4000 -C2 -P2 -m100 ./reference/tdb  ./%s/Nib_q/%s > %s/maf/%s.maf \n' %(q_name,name2,q_name,exname2))
            list3.write('last-split -m1 %s/maf/%s.maf | maf-convert psl | axtChain stdin ./reference/Nib_t/ ./%s/Nib_q/ \
            ./%s/chain/%s.chain -linearGap=loose -psl  \n' %(q_name,exname2,q_name,q_name,exname2))
    list2.close()
    list3.close()

def others():    
    # Chaining
    os.system('find ./%s/chain/ -name "*.chain" | chainMergeSort -inputList=stdin | gzip -c >./%s/all.chain.gz'%(q_name,q_name))
    os.system('chainPreNet ./%s/all.chain.gz ./reference/t.sizes ./%s/q.sizes ./%s/all.pre.chain'%(q_name,q_name,q_name))
    # Netting
    os.system('chainNet ./%s/all.pre.chain -minSpace=1 ./reference/t.sizes ./%s/q.sizes stdout /dev/null | netSyntenic stdin  ./%s/noClass.net'%(q_name,q_name,q_name))
    # format transformation
    os.system('netToAxt ./%s/noClass.net ./%s/all.pre.chain ./reference/t.2bit ./%s/q.2bit stdout | axtSort stdin %s/sorted_t-q.axt'%(q_name,q_name,q_name,q_name))
    os.system('axtToMaf ./%s/sorted_t-q.axt ./reference/t.sizes ./%s/q.sizes ./%s/%s.maf'%(q_name,q_name,q_name,q_name))

    os.system('maf-convert sam ./%s/%s.maf >./%s/%s.sam'%(q_name,q_name,q_name,q_name))
    os.system('samtools view -bt ./reference/t.fa.fai ./%s/%s.sam>./%s/%s.bam'%(q_name,q_name,q_name,q_name))
    os.system('samtools sort ./%s/%s.bam  ./%s/sort-%s'%(q_name,q_name,q_name,q_name))
    os.system('samtools index ./%s/sort-%s.bam'%(q_name,q_name))

    os.system('bamToBed -i ./%s/sort-%s.bam  > ./%s/%s.bed'%(q_name,q_name,q_name,q_name))
    os.system('sort -k1,1 -k2,2n ./%s/%s.bed > ./%s/sort-%s.bed'%(q_name,q_name,q_name,q_name))
    os.system('bedtools merge -i ./%s/sort-%s.bed > ./%s/merge-sort-%s.bed'%(q_name,q_name,q_name,q_name))
    
    # Screen out
    os.system('maf_sort ./%s/%s.maf  Osativa >./%s/%s.orig.maf'%(q_name,q_name,q_name,q_name))
    os.system('single_cov2 ./%s/%s.orig.maf R=%s >./reference/tba/Osativa.%s.sing.maf'%(q_name,q_name,reference_genome_name,q_name))

def main():
    deal_query()
    faSplit()
    faToNib()
    os.system('parallel -k -j %s <%s/fatonib.sh'%(threads,q_name))
    lastal()
    os.system('parallel -k -j %s <%s/list_of_lastal.sh'%(threads,q_name))
    os.system('parallel -k -j %s <%s/list_of_chain.sh'%(threads,q_name))
    others()

if __name__=='__main__': 
    main()
