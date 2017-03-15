#!/usr/lib/python
#-*- coding:utf-8 -*-

##-------------------------------------------------------------------------------------------------

##    This script refers to dealing with reference genome

## NB:based on last-759
##-------------------------------------------------------------------------------------------------
import sys
sys.path.append('/share/workplace/software/python-packages/biopython')
from Bio import SeqIO
import os
from datetime import *


target=sys.argv[1]
t_abspath=os.path.abspath(target)
t_fullname=os.path.basename(target)
t_name=os.path.splitext(t_fullname)[0]

os.system('mkdir -p reference/Nib_t reference/tba')

def deal_target():
    tf1=open(target,'r')
    tf2=open('./reference/t1.fa','w')
    tf22=open('./reference/tba/%s.fa'%(t_name),'w')
    for seq_record in SeqIO.parse(tf1,'fasta'):
	id=str(seq_record.id)
        seq=str(seq_record.seq)
        length=str(len(seq))
        tf2.write('>'+ t_name + '.'+ id +'\n'+ seq +'\n')
        tf22.write('>'+ t_name + ':'+ id + ':' + '1' + ':' + '+' + ':' + length + '\n'+ seq +'\n')
    tf1.close()
    tf2.close()
    tf22.close()	
    os.system('fold -w 70 ./reference/t1.fa > ./reference/t.fa') 
    os.system('fold -w 70 ./reference/tba/%s.fa >./reference/tba/%s'%(t_name,t_name)) 
    os.system('rm ./reference/t1.fa ./reference/tba/%s.fa'%(t_name)) 
	
def lastdb():
    os.system('/share/workplace/home/ping/software/last-759/bin/lastdb -uMAM8 -P10 ./reference/tdb  ./reference/t.fa \n' )
    os.system('samtools faidx ./reference/t.fa \n' )

def trans():
    os.system('faSplit byname ./reference/t.fa ./reference/Nib_t/')
    os.system('faToTwoBit ./reference/t.fa ./reference/t.2bit') 
    os.system('for i in ./reference/Nib_t/*; do faToNib $i `echo $i | sed -e s/.fa/.nib/`; done')  
    os.system('faSize ./reference/t.fa -detailed > ./reference/t.sizes')

def main():
    start1=datetime.now()
    deal_target()
    lastdb()
    trans()
    end1=datetime.now()
    print 'lastdb uses'+ str(end1-start1) 

if __name__=='__main__': 
    start=datetime.now()
    print 'whole program starts at' + str(start)
    main()
    end=datetime.now()
    print 'whole program ends at' + str(end)
    print 'totally uses ' + str(end-start)	
	





