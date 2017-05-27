#!usr/lib/python
# -*- coding:utf-8 -*-

##------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------
## Usage:
##    python *.py gene_CDS.bed score.bed chromosome
##------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------

from __future__ import division
import sys
import os
import subprocess


db1=os.path.abspath(sys.argv[1]) 
fullname=os.path.basename(db1)
name1=os.path.splitext(fullname)[0]
bed=sys.argv[2] 
keyword=sys.argv[3] 

def head():
    f1=open(db1,'r')
    f3=open('%s.gene_dengfen'%(name1),'w')
    f4=open('%s.gene_ave_median'%(name1),'w')
    f11=f1.readlines()
    os.system('sort -k1,1 -k2,2n %s>sort-k2-%s'%(db1,name1))
    min0=subprocess.Popen('''head -n1 sort-k2-%s|awk '{print $2}' '''%(name1),shell=True,stdout=subprocess.PIPE)
    mini=min0.communicate()[0].strip('\n')
    os.system('sort -k1,1 -k3,3n %s>sort-k3-%s'%(db1,name1))
    maxi0=subprocess.Popen('''tail -n1 sort-k3-%s|awk '{print $3}' '''%(name1),shell=True,stdout=subprocess.PIPE)
    maxi=maxi0.communicate()[0].strip('\n')
    os.system('''grep %s$"\t" %s|awk -F '\t' '{if(%s<=$2&&$3<=%s)print $0}' >%s.part.bed'''%(keyword,bed,mini,maxi,name1))
    f2=open('%s.part.bed'%(name1),'r')
    f22=f2.readlines()
    lines=len(f11)
    ac=''
    for i in range(lines):
        #chr_os=f11[i].split()[0]
        #start_os=int(f11[i].split()[1])
        #end_os=int(f11[i].split()[2])
        gene_ID=f11[i].split()[3]
        cds=f11[i].split()[5][:-1]
        dire=f11[i].split()[4]
        num=0
        list1=[]
        for k in cds.split(';'):
            begin=int(k.split(':')[0])
            end=int(k.split(':')[1])
            cha=end-begin
            num+=cha
            for a in range(cha):
                begin+=1
                list1.append(begin)
        if len(list1)<20:
            gene_score=[]
            for b in f22:
                position=b.strip().split('\t')[2]
                if int(position) <=list1[-1] and int(position) in list1:
                    score=float(b.strip().split('\t')[3])
                    gene_score.append(score)
            gene_average=sum(gene_score)/num
            for j in range(20):
                f3.write('gene'+str(j+1)+'\t'+str(gene_average)+'\n')
        else:    
            yu=len(list1)%20
            bases=(len(list1)-yu)/20
            l_qian=[bases+1]
            l_hou=[bases]
            l=l_qian*yu+l_hou*(20-yu)
            start_sites=0
            gene_score=[]
            for j in range(20):
                gene=list1[start_sites:start_sites+int(l[j])]
                maxi=gene[-1]
                list_of_score=[]
                start_sites=start_sites+int(l[j])
                for b in f22:
                    position=b.strip().split('\t')[2]
                    if int(position) <=maxi and int(position) in gene:
                        score=float(b.strip().split('\t')[3])
                        list_of_score.append(score)
                average=sum(list_of_score)/len(gene)
                gene_score+=list_of_score
                if dire=='+':
                    ac='gene'
                else:
                    ac='-gene'
                f3.write(ac+str(j+1)+'\t'+str(average)+'\n')

            gene_average=sum(gene_score)/len(list1)

        if len(gene_score)<=num:
            list11=gene_score+[0]*(num-len(gene_score))
            list2=sorted(list11)
        else:
            list2=sorted(gene_score)
        if len(list2)%2==1:
            median=list2[int((len(list2)+1)/2)]
        else:
            median=(list2[int(len(list2)/2)-1]+list2[int(len(list2)/2)])/2
        f4.write(gene_ID+'\t'+str(gene_average)+'\t'+str(median)+'\n')


    os.system('rm %s.part.bed sort-k2-%s sort-k3-%s'%(name1,name1,name1))
    f1.close()
    f2.close()
    f3.close()

if __name__=="__main__":
    head()
