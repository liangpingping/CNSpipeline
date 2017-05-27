#!/usr/lib/python
#-*- coding:utf-8 -*-

import os
import sys

target=sys.argv[1]
q_name=os.path.basename(target)
name=os.path.splitext(q_name)[0]
os.system(''' grep -v ^$ %s|grep -v 'score'|awk '{print $2"\t"$7}' >%s_qu_kong''' %(target,name))
 
def hebing():
    f1=open('%s_qu_kong'%(name),'r')
    f2=open('%s.fa'%(name),'w')
    f11=f1.readlines()
    os=''
    oi=''
    ob=''
    tae=''
    tu=''
    ata=''
    hv=''
    bd=''
    phe=''
    sv=''
    si=''
    pha=''
    pv=''
    sb=''
    zm=''
    et=''
    oth=''
    aco=''
    pda=''
    eg=''
    ma=''
    sp=''
    atr=''

    for i in f11:
        name1=i.split('\t')[0]
        seq=i.rstrip().split('\t')[1]
        if 'Osativa' in name1:
            os+=seq
        if 'Oindica' in name1:
            oi+=seq
        if 'Obrachyantha' in name1:
            ob+=seq
        if 'Taestivum' in name1:
            tae+=seq
        if 'Turartu' in name1:
            tu+=seq
        if 'Atauschii' in name1:
            ata+=seq
        if 'Hvulgare' in name1:
            hv+=seq
        if 'Bdistachyon' in name1:
            bd+=seq
        if 'Pheterocycla' in name1:
            phe+=seq
        if 'Sviridis' in name1:
            sv+=seq
        if 'Sitalica' in name1:
            si+=seq
        if 'Phallii' in name1:
            pha+=seq
        if 'Pvirgatum' in name1:
            pv+=seq
        if 'Sbicolor' in name1:
            sb+=seq
        if 'Zmays' in name1:
            zm+=seq
        if 'Etef' in name1:
            et+=seq
        if 'Othomaeum' in name1:
            oth+=seq
        if 'Acomosus' in name1:
            aco+=seq
        if 'Pdactylifera' in name1:
            pda+=seq
        if 'Eguineensis' in name1:
            eg+=seq
        if 'Macuminata' in name1:
            ma+=seq
        if 'Spolyrhiza' in name1:
            sp+=seq
        if 'Atrichopoda' in name1:
            atr+=seq
    print len(tae.replace('-','')),len(tu.replace('-','')),len(ata.replace('-','')),len(hv.replace('-','')),len(bd.replace('-','')),len(phe.replace('-','')),len(ob.replace('-','')),len(oi.replace('-','')),len(os.replace('-','')),len(et.replace('-','')),len(oth.replace('-','')),len(sb.replace('-','')),len(zm.replace('-','')),len(sv.replace('-','')),len(si.replace('-','')),len(pha.replace('-','')),len(pv.replace('-','')),len(aco.replace('-','')),len(ma.replace('-','')),len(eg.replace('-','')),len(pda.replace('-','')),len(sp.replace('-','')),len(atr.replace('-',''))
    print 'tae,tu,ata,hv,bd,phe,ob,oi,os,et,oth,sb,zm,sv,si,pha,pv,aco,ma,eg,pda,sp,atr'
    f2.write('>T.aestivum\n'+tae+'\n'+'>T.urartu\n'+tu+'\n'+'>A.tauschii\n'+ata+'\n'+'>H.vulgare\n'+hv+'\n'+'>B.distachyon\n'+bd+'\n'+'>P.heterocycla\n'+phe+'\n'+'>O.brachyantha\n'+ob+'\n'+'>O.indica\n'+oi+'\n'+'>O.sativa\n'+os+'\n'+'>E.tef\n'+et+'\n'+'>O.thomaeum\n'+oth+'\n'+'>S.bicolor\n'+sb+'\n'+'>Z.mays\n'+zm+'\n'+'>S.viridis\n'+sv+'\n'+'>S.italica\n'+si+'\n'+'>P.hallii\n'+pha+'\n'+'>P.virgatum\n'+pv+'\n'+'>A.comosus\n'+aco+'\n'+'>M.acuminata\n'+ma+'\n'+'>E.guineensis\n'+eg+'\n'+'>P.dactylifera\n'+pda+'\n'+'>S.polyrhiza\n'+sp+'\n'+'>A.trichopoda\n'+atr+'\n')
    f1.close()
    f2.close()

if __name__=='__main__':
    hebing()
    os.system('rm %s_qu_kong'%(name))
                                           
