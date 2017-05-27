##=======================================================================================================
##  This script used to extract shared cds and shared noncoding sequences from maf file.
##  Usage:
##  python script.py *.maf cdsfile 
##  It will generates the results named *.cds_shared,*.noncoding_shared 
##=======================================================================================================

#!/usr/lib/python
#-*- coding:utf-8 -*-

import sys
import os

target=sys.argv[1]
cdsfile=sys.argv[2]
q_name=os.path.basename(target)
name=os.path.splitext(q_name)[0]

os.system('maf-sort %s >sort-%s'%(target,name))
os.system('''grep 'Osativa' sort-%s|awk '{print $2"\t"$3"\t"$3+$4}' >%s-Osativa.bed '''%(name,name))
os.system('bedtools intersect -a %s-Osativa.bed -b %s >%s-intersect.bed'%(name,cdsfile,name))
os.system('sort -k1,1 -k2,2n %s-intersect.bed >sort-%s-intersect.bed'%(name,name))
os.system("sed -i '1d' sort-%s "%(name))
os.system("sed -i '$d' sort-%s "%(name))
os.system("sed -i '/^#/d' sort-%s "%(name))


def extract():
    f1=open('%s-Osativa.bed'%(name),'r')
    f2=open('sort-%s-intersect.bed'%(name),'r')
    f3=open('%s.index'%(name),'w')
    f11=f1.readlines()
    f22=f2.readlines()
    for i in f11:
        chri=i.rstrip().split('\t')[0]
        starti=i.rstrip().split('\t')[1]
        endi=i.rstrip().split('\t')[2]
        a=''
        for j in f22:
            chrj=j.rstrip().split('\t')[0]
            startj=j.rstrip().split('\t')[1]
            endj=j.rstrip().split('\t')[2]
            if chrj==chri and int(starti) <= int(startj) and int(endj)<=int(endi):  ## because cds file intersect with os.bed,so endj will less than endi if starti < startj <= endi
                a+=startj+','+endj+';'
        f3.write(i.rstrip()+'\t'+a+'\n')
    f1.close()
    f2.close()
    f3.close()

def get_sequences():
    f1=open('%s.index'%(name),'r')
    f2=open('sort-%s'%(name),'r')
    f3=open('%s.cds_shared'%(name),'w')
    f4=open('%s.noncoding_shared'%(name),'w')
    f11=f1.readlines()
    f22=f2.read().split('\n\n')
    lines=len(f11)
    for i in range(lines):
        chri=f11[i].rstrip().split('\t')[0]
        starti=f11[i].rstrip().split('\t')[1]
        endi=f11[i].rstrip().split('\t')[2]
        length=len(f11[i].rstrip().split('\t'))

        list=f22[i].split('\n')
        list1=','.join(list[1].split())
        chr_os=list1.split(',')[1]
        start_os=list1.split(',')[2]
        match=list1.split(',')[3]
        end_os=str(int(start_os)+int(match))
        seq_os=list1.split(',')[6]
        if length==3:
            for a in list[1:]:
                a=','.join(a.split())
                if a != '':
                    name1=a.split(',')[1]
                    beg=a.split(',')[2]
                    num=a.split(',')[3]
                    seq1=a.split(',')[6]
                    f4.write(name1 +'\t'+beg+'\t'+str(int(beg)+int(num))+'\t'+seq1+'\n')

        if length==4:
            array=(f11[i].rstrip().split('\t')[3]+end_os+',').split(';')
	    if len(array)>2:
                kong=len(seq_os)-int(match)
                seq_os=list1.split(',')[6]
                bstart=array[0].split(',')[0]
                ran0=int(bstart)-int(start_os)
                b0=ran0
                l0=len(seq_os[:b0].replace('-',''))
                while l0<ran0:
                    b0+=1
                    l0=len(seq_os[:b0].replace('-',''))
                s_non0=seq_os[:b0]
                if s_non0 != '':
                    f4.write(chr_os+'\t'+start_os+'\t'+str(int(start_os)+ran0)+'\t'+s_non0+'\n')
                    for c in list[2:]:
                        c=','.join(c.split())
                        name2=c.split(',')[1]
                        beg2=c.split(',')[2]
                        num2=c.split(',')[3]
                        seq2=c.split(',')[6][:b0]
                        f4.write(name2+'\t'+beg2+'\t'+str(int(beg2)+ran0)+'\t'+seq2+'\n')

                for b in range(len(array)-2):
                    startb=array[b].split(',')[0]
                    endb=array[b].split(',')[1]
                    startb1=array[b+1].split(',')[0]
                    ran1=int(endb)-int(start_os)
                    b1=ran1
		    l1=len(seq_os[:b1].replace('-',''))
                    while l1<ran1:
                        b1+=1
                        l1=len(seq_os[:b1].replace('-',''))
                                        
		    ranb1=int(startb1)-int(start_os)
                    b1b=ranb1
                    l1b=len(seq_os[:b1b].replace('-',''))
		    while l1b<ranb1:
                        b1b+=1
                        l1b=len(seq_os[:b1b].replace('-',''))
                    s_non=seq_os[b1:b1b]
          	    s_cds=seq_os[b0:b1]
                    if s_non != '':
                        f4.write(chr_os+'\t'+str(int(start_os)+ran1)+'\t'+str(int(start_os)+ranb1)+'\t'+s_non+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b1:b1b]
				f4.write(name2+'\t'+str(int(beg2)+ran1)+'\t'+str(int(beg2)+ranb1)+'\t'+seq2+'\n')
                    if s_cds != '':
	      	        f3.write(chr_os+'\t'+str(int(start_os)+ran0)+'\t'+str(int(start_os)+ran1)+'\t'+s_cds+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
				name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b0:b1]
                                f3.write(name2+'\t'+str(int(beg2)+ran0)+'\t'+str(int(beg2)+ran1)+'\t'+seq2+'\n')
                    b0=b1b
                last_start=array[-2].split(',')[0]
                last_end=array[-2].split(',')[1]
                last_start1=array[-1].split(',')[0]
                if int(last_end)==int(last_start1):
                    last_cds=seq_os[b1b:]
                    if last_cds != '':
                        f3.write(chr_os+'\t'+str(int(start_os)+ranb1)+'\t'+str(int(start_os)+int(match))+'\t'+last_cds+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b1b:]
                                f3.write(name2+'\t'+str(int(beg2)+ranb1)+'\t'+str(int(beg2)+int(num2))+'\t'+seq2+'\n')
                else:
                    jiange1=int(last_end)-int(start_os)
                    bc1=jiange1
                    l1bc1=len(seq_os[:bc1].replace('-',''))
                    while l1bc1<jiange1:
                        bc1+=1
                        l1bc1=len(seq_os[:bc1].replace('-',''))
                    last_cds=seq_os[b0:bc1]
                    last_non=seq_os[bc1:]
                    if last_non != '':
                        f4.write(chr_os+'\t'+str(int(start_os)+int(jiange1))+'\t'+str(int(start_os)+int(match))+'\t'+last_non+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][bc1:]
                                f4.write(name2+'\t'+str(int(beg2)+int(jiange1))+'\t'+str(int(beg2)+int(num2))+'\t'+seq2+'\n')
                    if last_cds != '':
                        f3.write(chr_os+'\t'+str(int(start_os)+ranb1)+'\t'+str(int(start_os)+int(jiange1))+'\t'+last_cds+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b0:bc1]
                                f3.write(name2+'\t'+str(int(beg2)+ranb1)+'\t'+str(int(beg2)+int(jiange1))+'\t'+seq2+'\n')
			
			
	    if len(array)==2:
                kong=len(seq_os)-int(match)
                seq_os=list1.split(',')[6]
                bstart=array[0].split(',')[0]
      	        bend=array[0].split(',')[1]
	        start2=array[1].split(',')[0]
                ran0=int(bstart)-int(start_os)
                b0=ran0
                l0=len(seq_os[:b0].replace('-',''))
                while l0<ran0:
                    b0+=1
                    l0=len(seq_os[:b0].replace('-',''))
                s_non0=seq_os[:b0]
                if s_non0 != '':
                    f4.write(chr_os+'\t'+start_os+'\t'+str(int(start_os)+ran0)+'\t'+s_non0+'\n')
                    for c in list[2:]:
                        c=','.join(c.split())
                        if c != '':
                            name2=c.split(',')[1]
                            beg2=c.split(',')[2]
                            num2=c.split(',')[3]
                            seq2=c.split(',')[6][:b0]
                            f4.write(name2+'\t'+beg2+'\t'+str(int(beg2)+ran0)+'\t'+seq2+'\n')

                ran1=int(bend)-int(start_os)
                b1=ran1
                l1=len(seq_os[:b1].replace('-',''))
                while l1<ran1:
                    b1+=1
                    l1=len(seq_os[:b1].replace('-',''))
                if int(bend)==int(start2):
                    s_cds0=seq_os[b0:]
                    if s_cds0 != '':
                        f3.write(chr_os+'\t'+str(int(start_os)+ran0)+'\t'+str(int(start_os)+int(match))+'\t'+s_cds0+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b0:]
                                f3.write(name2+'\t'+str(int(beg2)+ran0)+'\t'+str(int(beg2)+int(num2))+'\t'+seq2+'\n')

                if int(bend)!=int(start2):
                    s_cds0=seq_os[b0:b1]
                    if s_cds0 != '':
                        f3.write(chr_os+'\t'+str(int(start_os)+ran0)+'\t'+str(int(start_os)+ran1)+'\t'+s_cds0+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b0:b1]
                                f3.write(name2+'\t'+str(int(beg2)+ran0)+'\t'+str(int(beg2)+ran1)+'\t'+seq2+'\n')
                    non1=seq_os[b1:]
                    if non1 != '':
                        f4.write(chr_os+'\t'+str(int(start_os)+ran1)+'\t'+str(int(start_os)+int(match))+'\t'+non1+'\n')
                        for c in list[2:]:
                            c=','.join(c.split())
                            if c != '':
                                name2=c.split(',')[1]
                                beg2=c.split(',')[2]
                                num2=c.split(',')[3]
                                seq2=c.split(',')[6][b1:]
                                f4.write(name2+'\t'+str(int(beg2)+ran1)+'\t'+str(int(num2)+int(beg2))+'\t'+seq2+'\n')

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    

if __name__=='__main__':
    extract()
    get_sequences()
