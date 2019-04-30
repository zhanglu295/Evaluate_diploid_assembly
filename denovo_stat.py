#!/usr/bin/env python
# coding: utf-8

# ### Import essential packages

# In[1]:
'''
from collections import defaultdict
import matplotlib
import re
import numpy
import operator
from scipy import stats as sta
import seaborn as sns
import math
from collections import defaultdict
from scipy.interpolate import spline
import operator
from itertools import groupby
import sys
import gzip
import getopt
'''
# ## Input Contig alignment

# In[2]:

class aligncontig(object):#class 'aligncontig' to store alignment information for contigs,
          def __init__(self):
            self.start=0 # start of contig alignment ,
            self.end=0   # end of contig alignment ,
            self.chrid='' # total length including gap,
            self.alignlength=0 #total length without misassemblies,
            self.forward=1
            self.misstart=[] #list for the starts of aligned contigs, which may include the breakpoints of misassemblies,
            self.misend=[]   #list for the ends of aligned contigs, which may include the breakpoints of misassemblies,
def alignment_contig_coordinate(tsv_file,chr_id,chr_len):
        f_tsv = open(tsv_file,'r')
        count = 0
        contig_dict = defaultdict(aligncontig)
        start_align=[]
        end_align=[]
        start_contig=[]
        end_contig=[]
        offset_array=[]
        offset_array.append(0)
        offset=0
        chrid=''
        chridlist=[]
        for length_chr in chr_len:
            offset=offset+length_chr
            offset_array.append(offset)
        alignlength=0
        processid=0
        for line in f_tsv:
            data=line.rsplit()
            if count==0:
                count+=1
                continue
            elif count%2 == 1:
                data=line.rsplit()
                if len(data)==4:
                    continue
                index_offset=chr_id.index(data[4])
                offset_add=offset_array[index_offset]
                if processid==0:
                    processid=int(data[5])
                if processid==0 or processid==int(data[5]):
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                    chridlist.append(data[4])
                else:
                    candlen=defaultdict(int)
                    for index in range(len(chridlist)):
                         candlen[chridlist[index]]+=abs(int(end_align[index])-int(start_align[index]))+1
                    maxchr=''
                    maxvalue=0
                    for key,value in candlen.items():
                        if value>maxvalue:
                            maxvalue=value
                            maxchr=key
                    start_align_final=[]
                    end_align_final=[]
                    alignlength_final=0
                    chrid_final=''
                    start_contig_final=[]
                    end_contig_final=[]
                    chrid_final=maxchr
                    for index in range(len(chridlist)):
                        if chridlist[index]==maxchr:
                            start_align_final.append(start_align[index])
                            end_align_final.append(end_align[index])
                            alignlength_final=alignlength_final+end_align[index]-start_align[index]+1
                            start_contig_final.append(start_contig[index])
                            end_contig_final.append(end_contig[index])
                    
                    Contig=aligncontig()
                    Contig.start=min(start_align_final)
                    Contig.end=max(end_align_final)
                    Contig.misstart=start_align_final
                    Contig.misend=end_align_final
                    Contig.alignlength=alignlength_final
                    Contig.chrid=chrid_final
                    length_forward=0
                    length_reverse=0
                    for index in range(len(start_contig_final)):
                        if start_contig_final[index]<end_contig_final[index]:
                            length_forward+=end_contig_final[index]-start_contig_final[index]
                        else:
                            length_reverse+=start_contig_final[index]-end_contig_final[index]
                    if length_forward>=length_reverse:
                        Contig.forward=1
                    else:
                        Contig.forward=0
                    contig_dict[processid]=Contig
                    start_align=[]
                    end_align=[]
                    start_contig=[]
                    end_contig=[]
                    chridlist=[]
                    alignlength=0
                    processid=int(data[5])
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[1])-int(data[0]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[0])-int(data[1]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3])) 
                    chridlist.append(data[4])
            count=count+1

        candlen=defaultdict(int)
        for index in range(len(chridlist)):
            candlen[chridlist[index]]+=abs(int(end_align[index])-int(start_align[index]))+1
            maxchr=''
            maxvalue=0
            for key,value in candlen.items():
                if value>maxvalue:
                    maxvalue=value
                    maxchr=key
            start_align_final=[]
            end_align_final=[]
            alignlength_final=0
            chrid_final=''
            start_contig_final=[]
            end_contig_final=[]
            chrid_final=maxchr
            for index in range(len(chridlist)):
                if chridlist[index]==maxchr:
                    start_align_final.append(start_align[index])
                    end_align_final.append(end_align[index])
                    alignlength_final=alignlength_final+end_align[index]-start_align[index]+1
                    start_contig_final.append(start_contig[index])
                    end_contig_final.append(end_contig[index])
        Contig=aligncontig()
        Contig.start=min(start_align_final)
        Contig.end=max(end_align_final)
        Contig.misstart=start_align_final
        Contig.misend=end_align_final
        Contig.alignlength=alignlength_final
        Contig.chrid=chrid_final
        length_forward=0
        length_reverse=0
        for index in range(len(start_contig_final)):
            if start_contig_final[index]<end_contig_final[index]:
                length_forward+=end_contig_final[index]-start_contig_final[index]
            else:
                length_reverse+=start_contig_final[index]-end_contig_final[index]
        if length_forward>=length_reverse:
            Contig.forward=1
        else:
            Contig.forward=0
        contig_dict[processid]=Contig
        return contig_dict


# Input contig alignment breaked by misassembly and indel >5 

def alignment_miss_contig(tsv_file,chr_id,chr_len):
        f_tsv = open(tsv_file,'r')
        count = 0
        contig_dict = defaultdict(aligncontig)
        start_align=[]
        end_align=[]
        start_contig=[]
        end_contig=[]
        alignlength=0
        processid=0
        processid_mis=0
        contig_status=1
        indexid=0
        offset_array=[]
        offset_array.append(0)
        offset=0
        for length_chr in chr_len:
            offset=offset+length_chr
            offset_array.append(offset)
        for line in f_tsv:
            if count==0:
                count+=1
                continue
            elif count%2==0:
                data=line.rsplit()
                if data[0]=='relocation,' or line=='indel: indel (> 5bp)' or data[0]=='inversion' or data[0]=='translocation':
                    contig_status=0            
            elif count%2 == 1:
                data=line.rsplit()
                if len(data)==4:
                    continue
                index_offset=chr_id.index(data[4])
                offset_add=offset_array[index_offset]
                if processid==0:
                    processid=int(data[5])
                if contig_status==1 and processid==int(data[5]):
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[1])-int(data[0]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[0])-int(data[1]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                if processid!=int(data[5]) or contig_status==0:
                    contig_status=1
                    Contig=aligncontig()
                    Contig.start=min(start_align)
                    Contig.end=max(end_align)
                    Contig.alignlength=alignlength
                    if start_contig[0]<end_contig[0]:
                        Contig.forward=1
                    else:
                        Contig.forward=0
                    Contig.misstart=start_align
                    Contig.misend=end_align
                    contig_dict[processid_mis]=Contig
                    processid_mis+=1
                    start_align=[]
                    end_align=[]
                    start_contig=[]
                    end_contig=[]
                    alignlength=0
                    processid=int(data[5])
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[1])-int(data[0]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[0])-int(data[1]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
            count=count+1
        contig=aligncontig()
        Contig.start=min(start_align)
        Contig.end=max(end_align)
        if start_contig[0]<end_contig[0]:
               Contig.forward=1
        else:
               Contig.forward=0
        Contig.misstart=start_align
        Contig.misend=end_align
        Contig.alignlength=alignlength
        contig_dict[processid_mis]=Contig
        return contig_dict


# ## Calculate Contig Statistics
def CalculateNx_contig_align(contig_dict):
    NX_align_contig=[]
    contig_align_list=[]
    contig_align_list_new=[]
    for k,v in contig_dict.items():
       contig_align_list.append(v.alignlength)
    lastlength=0
    contig_align_list_new=sorted(contig_align_list,reverse=True)
    totallen=sum(contig_align_list_new)
    accu_len=0
    for contiglen in contig_align_list_new:
        lastlength=accu_len
        accu_len=accu_len+contiglen
        if accu_len>=0.1*totallen and lastlength<0.1*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.2*totallen and lastlength<0.2*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.3*totallen and lastlength<0.3*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.4*totallen and lastlength<0.4*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.5*totallen and lastlength<0.5*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.6*totallen and lastlength<0.6*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.7*totallen and lastlength<0.7*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.8*totallen and lastlength<0.8*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.9*totallen and lastlength<0.9*totallen:
            NX_align_contig.append(contiglen/1000)
    return NX_align_contig

def CalculateNx_contig(infile,threshold):
    NX_contig=[]
    lastlength=0
    if '.gz' in infile:
         contig_fasta=gzip.open(infile,"r")
    else:
        contig_fasta=open(infile,"r")
    contig_list=[]
    lengthlist=defaultdict(int)
    contigid=''
    contigseq=0
    index=0
    for contigx in contig_fasta:
        if '.gz' in infile:
            contig=contigx.decode()
        else:
            contig=contigx 
        if contig[0]!='>':
            contig_clean=contig.strip('\n')
            contigseq=contigseq+len(contig_clean)
        else:
            if index!=0:
                lengthlist[contigid]=contigseq
            if contigseq>threshold*1000 and index!=0:
                contig_list.append(contigseq)    
            contigseq=0
            contigid=contig.strip('\n')[1:]
        index+=1
    if contigseq>threshold*1000:
        contig_list.append(contigseq)
    lengthlist[contigid]=contigseq
    contig_list_new=sorted(contig_list,reverse=True)
    totallen=sum(contig_list_new)
    accu_len=0
    #out1=open('Raw_len','w')
    for contiglen in contig_list_new:
        #out1.write(str(contiglen)+'\n')
        lastlength=accu_len
        accu_len=accu_len+contiglen
        if accu_len>=0.1*totallen and lastlength<0.1*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.2*totallen and lastlength<0.2*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.3*totallen and lastlength<0.3*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.4*totallen and lastlength<0.4*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.5*totallen and lastlength<0.5*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.6*totallen and lastlength<0.6*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.7*totallen and lastlength<0.7*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.8*totallen and lastlength<0.8*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.9*totallen and lastlength<0.9*totallen:
            NX_contig.append(contiglen/1000)
    #out1.close()
    return NX_contig,lengthlist
def CalculateNGX_contig(lengthlist,reflen):
    contig_sort=sorted(list(lengthlist.values()),reverse=True)
    accu_len=0
    NGX_contig=[]
    for contiglen in contig_sort:
        lastlength=accu_len
        accu_len=accu_len+contiglen
        if accu_len>=0.1*reflen and lastlength<0.1*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.2*reflen and lastlength<0.2*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.3*reflen and lastlength<0.3*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.4*reflen and lastlength<0.4*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.5*reflen and lastlength<0.5*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.6*reflen and lastlength<0.6*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.7*reflen and lastlength<0.7*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.8*reflen and lastlength<0.8*reflen:
            NGX_contig.append(contiglen/1000)
        if accu_len>=0.9*reflen and lastlength<0.9*reflen:
            NGX_contig.append(contiglen/1000)
    return NGX_contig
def CalculateNAx_contig(misassemblist):
    lengthlist=[]
    NAX_contig=[]
    lastlength=0
    #out2=open('cutlength','w')
    for key,value in misassemblist.items():
        lengthlist.append(value.alignlength)
        #out2.write(str(key)+'\t'+str(value.alignlength)+'\n')
    AA=sorted(lengthlist,reverse=True)
    totallength=sum(AA)
    length_agg=0
    #out2=open('cutlength','w')
    for line in AA:
        #out2.write(str(line)+'\n')
        lastlength=length_agg
        length_agg=length_agg+line
        if length_agg>=0.1*totallength and lastlength<0.1*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.2*totallength and lastlength<0.2*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.3*totallength and lastlength<0.3*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.4*totallength and lastlength<0.4*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.5*totallength and lastlength<0.5*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.6*totallength and lastlength<0.6*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.7*totallength and lastlength<0.7*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.8*totallength and lastlength<0.8*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.9*totallength and lastlength<0.9*totallength:
            NAX_contig.append(line/1000)
    #out2.close()
    return NAX_contig
# ## Get chrid and chrlen

class contigcov:#class to store contig length and coverage,
        def __init__(self,length,coverage):
            self.length=length
            self.coverage=coverage
class uncovered:#class for uncovered genomic regions,
        def __init__(self, start, end, length):
            self.start=start
            self.end=end
            self.length=length
# In[20]:

def CalculateNCx_contig(alignment,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        alignlen=[]
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        for key,value in alignment.items():
            alignlen.append(value)
        alignlen.sort(key=operator.attrgetter('alignlength'),reverse=True)
        index=0
        for contig in alignlen:
            index+=1
            start=contig.misstart
            end=contig.misend
            if contig.start==laststart and contig.end==lastend:
                continue
            else:
                laststart=contig.start
                lastend=contig.end
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold: 
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)                                                                                                           
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)   
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover 
            list_contig_cov.append(contigcov(contig.alignlength/1000,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(contig.alignlength/1000)
        return NCx,list_contig_cov,round(nowlength/ref_len,4)


# ## Input Scaffold alignment

# In[21]:

class scaffolding(object):
          def __init__(self):
            self.scaffold_len=0
            self.contig_num=0
            self.pro_correct=0
            self.pro_error=0
            self.pro_skip=0
            self.contigalign_start=[]
            self.contigalign_end=[]
def analyze_scaffold_info(info_file,align_contig):
        f_info = open(info_file,"r")
        scaffold_dict1 = defaultdict(list) #contig id,
        scaffold_dict2 = defaultdict(list) #the order of contigs in scaffold based on assembly,
        scaffold_dict3 = defaultdict(list) #contig end position,
        scaffold_dict4 = defaultdict(list) #the order of contigs in scaffold based on alignment,
        scaffold_dict5 = defaultdict(list) #contig alignment length,
        scaffold_dict6 = defaultdict(list) #contig start list,
        scaffold_dict7 = defaultdict(list) #contig end list,
        total=0 
        for line in f_info:
            data = line.rsplit()
            scaffold_id = int(data[0])
            contig_id = int(data[1])
            contig_order = int(data[2])
            if contig_id in list(align_contig.keys()):
                scaffold_dict1[scaffold_id].append(contig_id)
                scaffold_dict2[scaffold_id].append(contig_order)
                scaffold_dict3[scaffold_id].append(align_contig[contig_id].end)
                scaffold_dict5[scaffold_id].append(align_contig[contig_id].alignlength)
                scaffold_dict6[scaffold_id].extend(align_contig[contig_id].misstart)
                scaffold_dict7[scaffold_id].extend(align_contig[contig_id].misend)
        scaffoldlist=[]
        for key,value in scaffold_dict3.items():
            if len(value) > 1:
                scaffold_dict4[key] = sorted(value)
                scaffold_info=scaffolding()
                scaffold_info.contig_num=len(scaffold_dict3[key])
                scaffold_info.contigalign_start=scaffold_dict6[key]
                scaffold_info.contigalign_end=scaffold_dict7[key]
                for i in range(len(scaffold_dict3[key])): 
                    scaffold_info.scaffold_len=scaffold_info.scaffold_len+scaffold_dict5[key][i]/1000000
                total+=scaffold_info.scaffold_len
            else:
                scaffold_info=scaffolding()
                scaffold_info.contigalign_start=scaffold_dict6[key]
                scaffold_info.contigalign_end=scaffold_dict7[key]
                scaffold_info.contig_num=len(scaffold_dict3[key])
                scaffold_info.scaffold_len=scaffold_info.scaffold_len+scaffold_dict5[key][-1]/1000000
                total+=scaffold_info.scaffold_len
            scaffoldlist.append(scaffold_info)
        return scaffoldlist


# ### block basic information (break scaffold by misscaffolding)

class block(object):
          def __init__(self):
            self.block_num=0
            self.block_alignlen=0
            self.block_startlist=[]
            self.block_endlist=[]
def analyze_block_info(info_file,alignment,lengthlist):
        f_info = open(info_file,"r")
        scaffold_dict1 = defaultdict(list)
        scaffold_dict2 = defaultdict(list)
        scaffold_dict3 = defaultdict(list)
        scaffold_dict4 = defaultdict(list)
        scaffold_dict5 = defaultdict(list)
        scaffold_dict6 = defaultdict(list)
        scaffold_dict7 = defaultdict(list)
        scaffold_dict8 = defaultdict(list)
        for line in f_info:
            data = line.rsplit()
            scaffold_id = int(data[0])
            contig_id = int(data[1])
            gap_size=int(data[3])
            scaffold_dict8[scaffold_id].append(int(gap_size))
            if contig_id in list(alignment.keys()):
                contig_order = int(data[2])
                scaffold_dict1[scaffold_id].append(contig_id)
                scaffold_dict2[scaffold_id].append(contig_order)
                scaffold_dict3[scaffold_id].append(alignment[contig_id].end)
                scaffold_dict5[scaffold_id].append(alignment[contig_id].alignlength)
                scaffold_dict6[scaffold_id].extend(alignment[contig_id].misstart)
                scaffold_dict7[scaffold_id].extend(alignment[contig_id].misend)
        block_list=[]
        #out=open('breakpoint_all_new','w')
        #out1=open('breakpoint_pos_new','w')
        #out2=open('breakpoint_stat_new','w')
        num_translocation=0
        num_relocation=0
        num_indel=0
        num_inversion=0
        num_overlap=0
        for key,value in scaffold_dict3.items():
                count_correct = 1
                contig_id_a = scaffold_dict1[key][0]
                block_len=alignment[contig_id_a].alignlength/1000000
                block_startlist=alignment[contig_id_a].misstart[:]
                block_endlist=alignment[contig_id_a].misend[:]
                scaffold_dict4[key] = sorted(value)
                if len(value)==1:
                    BL=block()
                    BL.block_num=count_correct
                    BL.block_alignlen=block_len
                    BL.block_startlist=block_startlist
                    BL.block_endlist=block_endlist
                    block_list.append(BL)
                    continue      
                repeat_length_list=[]   
                repeat_end_list=[]      
                for i in range(len(scaffold_dict3[key])-1):
                    contig_id_a = scaffold_dict1[key][i]
                    contig_id_b = scaffold_dict1[key][i+1]
                    a = scaffold_dict3[key][i]
                    b = scaffold_dict3[key][i+1]
                    idx_a = scaffold_dict4[key].index(a)
                    idx_b = scaffold_dict4[key].index(b)
                    contig_order_a=scaffold_dict2[key][i]
                    contig_order_b=scaffold_dict2[key][i+1]
                    include=0
                    #out.write(str(key)+'\t'+str(contig_id_a)+'\t'+str(contig_id_b)+'\t'+str(contig_order_a)+'\t'+str(contig_order_b)+'\t'+str(idx_a)+'\t'+str(idx_b)+'\n')
                    overlap=0
                    threshold1=alignment[contig_id_a].alignlength*0.95
                    threshold2=alignment[contig_id_b].alignlength*0.95
                    Start1=alignment[contig_id_a].misstart
                    End1=alignment[contig_id_a].misend
                    Start2=alignment[contig_id_b].misstart
                    End2=alignment[contig_id_b].misend
                    for pos1 in range(len(Start1)):
                        for pos2 in range(len(Start2)):
                           if Start2[pos2]<=Start1[pos1] and End2[pos2]>=Start1[pos1] and End2[pos2]<=End1[pos1]:
                               overlap+=End2[pos2]-Start1[pos1]+1
                           elif Start2[pos2]<=Start1[pos1] and End2[pos2]>=End1[pos1]:
                               overlap+=End1[pos1]-Start1[pos1]+1
                           elif Start2[pos2]>=Start1[pos1] and End2[pos2]<=End1[pos1]:
                               overlap+=End2[pos2]-Start2[pos2]+1
                           elif Start2[pos2]>=Start1[pos1] and Start2[pos2]<=End1[pos1] and End2[pos2]>=End1[pos1]:
                               overlap+=End1[pos1]-Start2[pos2]+1
                    order_contig=0
                    if alignment[contig_id_a].chrid!=alignment[contig_id_b].chrid:
                        num_translocation+=1
                    elif alignment[contig_id_a].forward!=alignment[contig_id_b].forward:
                        num_inversion+=1
                    if alignment[contig_id_a].chrid==alignment[contig_id_b].chrid and alignment[contig_id_a].forward==1 and alignment[contig_id_b].forward==1:
                        if idx_a>idx_b or idx_a+1<idx_b:
                            num_relocation+=1
                        if idx_a<idx_b:
                            if overlap>threshold1 or overlap>threshold2:
                                num_overlap+=1
                            if idx_a+1==idx_b and overlap<threshold1 and overlap<threshold2:
                                order_contig=1  
                            if order_contig==1:
                                if contig_order_a+1!=contig_order_b:
                                    sumdist_contig=0
                                    sumdist_gap=0
                                    indel=0
                                    for x in range(contig_id_b-contig_id_a-1):
                                        sumdist_contig=sumdist_contig+lengthlist[str(contig_id_a + x +1)]
                                        if lengthlist[str(contig_id_a + x +1)]>500:
                                            indel=1
                                    sumdist_gap=alignment[contig_id_b].start-alignment[contig_id_a].end+1
                                    if indel==1: 
                                        if sumdist_gap<sumdist_contig-2500 or sumdist_gap>sumdist_contig+2500:
                                            num_indel+=1
                    if alignment[contig_id_a].chrid==alignment[contig_id_b].chrid and alignment[contig_id_a].forward==0 and alignment[contig_id_b].forward==0:
                        if idx_a<idx_b or idx_a-1>idx_b:
                            num_relocation+=1
                        if idx_a>idx_b:
                            if overlap>threshold1 or overlap>threshold2:
                                num_overlap+=1
                        if idx_a-1==idx_b and overlap<threshold1 and overlap<threshold2:
                            order_contig=1
                        if order_contig==1:
                            if contig_order_a+1!=contig_order_b:
                                sumdist_contig=0
                                sumdist_gap=0
                                indel=0
                                for x in range(contig_id_b-contig_id_a-1):
                                    sumdist_contig=sumdist_contig+lengthlist[str(contig_id_a + x +1)]
                                    if lengthlist[str(contig_id_a + x +1)]>500:
                                        indel=1
                                sumdist_gap=alignment[contig_id_a].start-alignment[contig_id_b].end+1
                                if indel==1:
                                    if sumdist_gap<sumdist_contig-2500 or sumdist_gap>sumdist_contig+2500:
                                        num_indel+=1
                    order_contig=0
                    if alignment[contig_id_a].chrid==alignment[contig_id_b].chrid and alignment[contig_id_a].forward==1 and alignment[contig_id_b].forward==1:
                       if idx_a<idx_b:
                            if idx_a+1==idx_b and overlap<threshold1 and overlap<threshold2:
                                order_contig=1 
                            elif idx_a+1<idx_b:
                                 if len(set(scaffold_dict4[key][idx_a:idx_b]))==1:
                                     order_contig=1
                                     repeat_end_list.append(alignment[contig_id_a].end)
                                     repeat_length_list.append(alignment[contig_id_a].alignlength)
                       if order_contig==1:
                            if contig_order_a+1==contig_order_b:
                                include=1
                            else: 
                                sumdist_contig=0
                                sumdist_gap=0
                                indel=0
                                for x in range(contig_id_b-contig_id_a-1):
                                    sumdist_contig=sumdist_contig+lengthlist[str(contig_id_a + x +1)]
                                    if lengthlist[str(contig_id_a + x +1)]>500:
                                        indel=1
                                sumdist_gap=alignment[contig_id_b].start-alignment[contig_id_a].end+1
                                if indel==0:
                                    include=1
                                elif sumdist_gap>sumdist_contig-2500 and sumdist_gap<sumdist_contig+2500:
                                    include=1
                    elif alignment[contig_id_a].chrid==alignment[contig_id_b].chrid and alignment[contig_id_a].forward==0 and alignment[contig_id_b].forward==0:
                         if idx_a>idx_b:
                            if idx_a-1==idx_b and overlap<threshold1 and overlap<threshold2:
                                order_contig=1
                            elif idx_a-1>idx_b:
                                 if len(set(scaffold_dict4[key][idx_a:idx_b]))==1:
                                     order_contig=1
                                     repeat_end_list.append(alignment[contig_id_a].end)
                                     repeat_length_list.append(alignment[contig_id_a].alignlength)
                         if order_contig==1:
                            if contig_order_a+1==contig_order_b:
                                include=1
                            else:
                                sumdist_contig=0
                                sumdist_gap=0
                                indel=0
                                for x in range(contig_id_b-contig_id_a-1):
                                    sumdist_contig=sumdist_contig+lengthlist[str(contig_id_a + x +1)]
                                    if lengthlist[str(contig_id_a + x +1)]>500:
                                        indel=1
                                sumdist_gap=alignment[contig_id_a].start-alignment[contig_id_b].end+1
                                if indel==0:
                                    include=1
                                elif sumdist_gap>sumdist_contig-2500 and sumdist_gap<sumdist_contig+2500:
                                    include=1
                    if include==1:
                        count_correct += 1
                        block_len+=alignment[contig_id_b].alignlength/1000000
                        block_startlist.extend(alignment[contig_id_b].misstart)
                        block_endlist.extend(alignment[contig_id_b].misend) 
                    else:
                        #out1.write(str(contig_id_a)+'\t'+str(alignment[contig_id_a].chrid)+'\t'+str(a)+'\n')
                        #out.write('break\n')
                        BL=block()
                        BL.block_num=count_correct
                        BL.block_alignlen=block_len
                        BL.block_startlist=block_startlist
                        BL.block_endlist=block_endlist   
                        block_list.append(BL)
                        if contig_id_b in alignment:
                            count_correct=1
                            block_len=alignment[contig_id_b].alignlength/1000000
                            block_startlist=alignment[contig_id_b].misstart[:]
                            block_endlist=alignment[contig_id_b].misend[:]
                BL=block()
                BL.block_num=count_correct
                BL.block_alignlen=block_len
                BL.block_startlist=block_startlist
                BL.block_endlist=block_endlist  
                block_list.append(BL)
        #out2.write('translocation: ')
        #out2.write(str(num_translocation))
        #out2.write('\n')
        #out2.write('relocation: ')
        #out2.write(str(num_relocation))
        #out2.write('\n')
        #out2.write('indel: ')
        #out2.write(str(num_indel))
        #out2.write('\n')
        #out2.write('inversion: ')
        #out2.write(str(num_inversion))
        #out2.write('\n')
        #out2.write('overlap: ')
        #out2.write(str(num_overlap))
        #out2.write('\n')
        #out.close()
        #out1.close()
        #out2.close()
        return block_list


# ## Evaluation for assemble scaffold

# ### Calculate NX for Block (Mb)

def CalculateNx_block(blocklist):
        Nx=[]
        total_length=0
        subtotal_length=0
        blocklist.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        for block in blocklist:
            total_length+=block.block_alignlen
        for block in blocklist:
            subtotal_length+=block.block_alignlen
            temp=subtotal_length-block.block_alignlen
            if temp<total_length*0.1 and subtotal_length>=total_length*0.1:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.2 and subtotal_length>=total_length*0.2:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.3 and subtotal_length>=total_length*0.3:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.4 and subtotal_length>=total_length*0.4:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.5 and subtotal_length>=total_length*0.5:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.6 and subtotal_length>=total_length*0.6:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.7 and subtotal_length>=total_length*0.7:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.8 and subtotal_length>=total_length*0.8:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.9 and subtotal_length>=total_length*0.9:
                Nx.append(block.block_alignlen)
            if temp<total_length and subtotal_length>=total_length:
                Nx.append(block.block_alignlen)
        return Nx


# ### Calculate NX for Scaffold (Mb)

# In[24]:

def CalculateNx_scaffold_align(scaffoldlist):
        Nx=[]
        total_length=0
        subtotal_length=0
        scaffoldlist.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
        for scaffold in scaffoldlist:
            total_length+=scaffold.scaffold_len
        for scaffold in scaffoldlist:
            subtotal_length+=scaffold.scaffold_len
            temp=subtotal_length-scaffold.scaffold_len
            if temp<total_length*0.1 and subtotal_length>=total_length*0.1:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.2 and subtotal_length>=total_length*0.2:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.3 and subtotal_length>=total_length*0.3:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.4 and subtotal_length>=total_length*0.4:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.5 and subtotal_length>=total_length*0.5:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.6 and subtotal_length>=total_length*0.6:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.7 and subtotal_length>=total_length*0.7:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.8 and subtotal_length>=total_length*0.8:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.9 and subtotal_length>=total_length*0.9:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length and subtotal_length>=total_length:
                Nx.append(scaffold.scaffold_len)
        return Nx


# ### Calculate NGX for Scaffold

# In[25]:
'''
def CalculateNGx_scaffold(scaffoldlist,ref_len):
    NGx=[]
    subtotal_length=0
    scaffoldlist.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
    for scaffold in scaffoldlist:
        subtotal_length+=scaffold.scaffold_len
        temp=subtotal_length-scaffold.scaffold_len
        if temp<ref_len*0.1 and subtotal_length>=ref_len*0.1:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.2 and subtotal_length>=ref_len*0.2:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.3 and subtotal_length>=ref_len*0.3:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.4 and subtotal_length>=ref_len*0.4:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.5 and subtotal_length>=ref_len*0.5:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.6 and subtotal_length>=ref_len*0.6:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.7 and subtotal_length>=ref_len*0.7:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.8 and subtotal_length>=ref_len*0.8:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.9 and subtotal_length>=ref_len*0.9:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len and subtotal_length>=ref_len:
            NGx.append(scaffold.scaffold_len)
    return NGx
'''

# ### Calculate NGX for block

# In[26]:
'''
def CalculateNGx_block(blocklist,ref_len):
        NGx=[]
        total_length=0
        subtotal_length=0
        blocklist.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        for block in blocklist:
            subtotal_length+=block.block_alignlen
            temp=subtotal_length-block.block_alignlen
            if temp<ref_len*0.1 and subtotal_length>=ref_len*0.1:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.2 and subtotal_length>=ref_len*0.2:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.3 and subtotal_length>=ref_len*0.3:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.4 and subtotal_length>=ref_len*0.4:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.5 and subtotal_length>=ref_len*0.5:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.6 and subtotal_length>=ref_len*0.6:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.7 and subtotal_length>=ref_len*0.7:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.8 and subtotal_length>=ref_len*0.8:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.9 and subtotal_length>=ref_len*0.9:
                NGx.append(block.block_alignlen)
            if temp<ref_len and subtotal_length>=ref_len:
                NGx.append(block.block_alignlen)
        return NGx

'''
# ### Calculate NCX for scaffold

# In[27]:

def CalculateNCx_scaffold(align_scaffold,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        align_scaffold.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
        index=0
        A=0
        for scaffold in align_scaffold:
            index+=1
            start=scaffold.contigalign_start
            end=scaffold.contigalign_end
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold:
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)    
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover
            list_contig_cov.append(contigcov(scaffold.scaffold_len,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(scaffold.scaffold_len)
        return NCx,list_contig_cov


# ### Calculate NCX for block

# In[28]:

def CalculateNCx_block(align_block,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        align_block.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        index=0
        A=0
        for block in align_block:
            index+=1
            start=block.block_startlist
            end=block.block_endlist
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold:
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)  
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)   
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover 
            list_contig_cov.append(contigcov(block.block_alignlen,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(block.block_alignlen)
        return NCx,list_contig_cov
## reference id and length
def reflen_id(in_path_ref):
    gzornot=0
    if "gz" in in_path_ref:
        infile=gzip.open(in_path_ref,'r')
        gzornot=1
    else:
        infile=open(in_path_ref,'r')
    chrid=[]
    chrlen=[]
    singlen=0
    for line in infile:
        lineinfo=""
        if gzornot==1:
           lineinfo=line.decode()
        else:
           lineinfo=line
        if lineinfo[0]=='>':
            A_1=lineinfo[1:].replace('  ','__').strip('\n')
            A_2=A_1.replace(' ','_')
            B=A_2.replace(':','_')
            chrid.append(B)
            if singlen!=0:
                chrlen.append(singlen)
                singlen=0
        else:
            singlen=singlen+len(lineinfo[:])
    chrlen.append(singlen)
    return chrid,chrlen,sum(chrlen)/1000000
def helpinfo():
    helpinfo=\
    '''
        Calculate NX,NAX,NCX and NCAX for contigs and scaffolds
        Version: 1.0.0
        Dependents: Python (>=3.0)
        Last Updated Date: 2017-07-26
        Contact: zhanglu295@gmail.com

        Usage: python Cal_Denovo_Stat <options>

        Options:
                -t --tsv, tsv files located in ./contigs_reports/all_alignments_*.tsv from QUAST(multiple inputs separated by comma)
                -c --contig, compressed or uncompressed contig fasta files (multiple inputs separated by comma, only used to calculate Contig NX)
                -s --scaffold, conpressed or uncompressed scaffold fasta files (multiple inputs separated by comma, only used to calculate Scaffold NX)
                -r --reference, compressed or uncompressed reference genome fasta, if reference is missing, only alignment free statistics would be caluclated.
                -i --info, .info file(three columns: scaffold id, contig id, order of contig in the scaffold)
                -a --min_contig, the length of minimum contigs (kb)(only for calculating N50, the short contigs for calculating NA50,NC50 and NCA50 should be eliminated first by -m in QUAST)
                -b --min_scaffold, the length of minimum scaffold (kb)
                -p --prefix, the prefix of output file(default:Assembly_summary)
                -l --label, Human-readable assembly names (multiple labels separated by comma).

                -o --out, the path to output
                -h --help, help info
    '''
    print(helpinfo)
def main():
    tsvlist=[]
    contiglist=[]
    scaffoldlist=[]
    reference=None
    infolist=[]
    labellist=[]
    outpath=None
    min_contig=0
    min_scaffold=0
    prefix="Assembly_summary"
    opts, args = getopt.gnu_getopt(sys.argv[1:], 't:c:s:r:i:l:a:b:p:o:h', ['tsv', 'contig','scaffold' ,'reference', 'info' ,'label','min_contig','min_scaffold','prefix','out','help'])
    for o, a in opts:
        if o == '-t' or o == '--tsv':
                tsvlist=a.split(',') 
        if o == '-c' or o == '--contig':
                contiglist = a.split(',')
        if o == '-s' or o == '--scaffold':
                scaffoldlist = a.split(',')
        if o == '-r' or o == '--reference':
                reference = a
        if o == '-i' or o == '--info':
                infolist = a.split(',')
        if o == '-l' or o == '--label':
                labellist = a.split(',')
        if o =='-a' or o =='--min_contig':
                min_contig = float(a)
        if o =='-b' or o =='--min_scaffold':
                min_scaffold = float(a)
        if o =='-p' or o =='--prefix':
                prefix = a
        if o == '-o' or o == '--out':
                outpath = a
        if o == '-h' or o == '--help':
                helpinfo()
                sys.exit(-1)
    outfile=open(outpath+'/'+prefix+'.txt','w')
    if reference==None:
           outfile.write('Statistics')
           outfile.write('\t')
           if len(labellist)==0:
               for index in range(len(contiglist)):
                   outfile.write('Sample'+str(index+1)+'\t')
           else:
               for label in labellist:
                   outfile.write(label+'\t')
           outfile.write('\n')
           outfile.write('Contig N50(kb)'+'\t')
           No_contig=[]
           No_scaffold=[]
           mean_contig=[]
           mean_scaffold=[]
           max_contig=[]
           max_scaffold=[]    
           for contig in contiglist:
               Contig_NX,lengthlist_contig=CalculateNx_contig(contig,min_contig)              
               No_contig.append(len(list(lengthlist_contig.values())))
               mean_contig.append(np.mean(list(lengthlist_contig.values())))
               max_contig.append(max(list(lengthlist_contig.values())))
               outfile.write(str(round(Contig_NX[4],1))+'\t')
           outfile.write('\n')
           outfile.write('Number of contigs'+'\t')
           for contig in No_contig:
               outfile.write(str(contig)+'\t')
           outfile.write('\n')
           outfile.write('Mean contigs length(kb)'+'\t')
           for contig in mean_contig:
               outfile.write(str(round(contig/1000,1))+'\t')
           outfile.write('\n')
           outfile.write('Maximum contigs length(kb)'+'\t')
           for contig in max_contig:
               outfile.write(str(round(contig/1000,1))+'\t')
           outfile.write('\n')
           outfile.write('Scaffold N50(kb)'+'\t')
           for scaffold in scaffoldlist:
               Scaffold_NX,lengthlist_scaffold=CalculateNx_contig(scaffold,min_scaffold)
               No_scaffold.append(len(list(lengthlist_scaffold.values())))
               mean_scaffold.append(np.mean(list(lengthlist_scaffold.values())))
               max_scaffold.append(max(list(lengthlist_scaffold.values())))
               outfile.write(str(round(Scaffold_NX[4],1))+'\t')
           outfile.write('\n')
           outfile.write('Number of scaffolds'+'\t')
           for scaffold in No_scaffold:
               outfile.write(str(scaffold)+'\t')
           outfile.write('\n')
           outfile.write('Mean scaffolds length(kb)'+'\t')
           for scaffold in mean_scaffold:
               outfile.write(str(round(scaffold/1000,1))+'\t')
           outfile.write('\n')
           outfile.write('Maximum scaffolds length(kb)'+'\t')
           for scaffold in max_scaffold:
               outfile.write(str(round(scaffold/1000,1))+'\t')
           outfile.write('\n')
           outfile.close()
    else:
           No_contig=[]
           No_scaffold=[]
           mean_contig=[]
           mean_scaffold=[]
           max_contig=[]
           max_scaffold=[]
           NGX_contig=[]
           NGX_scaffold=[]
           outfile.write('Statistics')
           outfile.write('\t')
           if len(labellist)==0:
               for index in range(len(contiglist)):
                    outfile.write('Sample'+str(index+1)+'\t')
           else:
               for label in labellist:
                   outfile.write(label+'\t')
           outfile.write('\n')
           print('Reading refernece genome...')
           (chrid,chrlen,singlen)=reflen_id(reference)
           print('Calculationg Contig N50...')
           outfile.write('Contig N50(kb)'+'\t')
           for contig in contiglist:
               Contig_NX,lengthlist_contig=CalculateNx_contig(contig,min_contig)
               No_contig.append(len(list(lengthlist_contig.values())))
               mean_contig.append(np.mean(list(lengthlist_contig.values())))
               max_contig.append(max(list(lengthlist_contig.values())))
               NGX_contig.append(CalculateNGX_contig(lengthlist_contig,singlen*1000000))
               Contig_NX,lengthlist_contig=CalculateNx_contig(contig,min_contig)
               outfile.write(str(round(Contig_NX[4],1))+'\t')
           outfile.write('\n')
           print('Calculationg Contig NG50...')
           outfile.write('Contig NG50(kb)'+'\t')
           for NGX in NGX_contig:
               if len(NGX)<5:
                   outfile.write('0'+'\t')
               else:
                   outfile.write(str(round(NGX[4],1))+'\t')
           outfile.write('\n')
           outfile.write('Number of contigs'+'\t')
           for contig in No_contig:
               outfile.write(str(contig)+'\t')
           outfile.write('\n')
           outfile.write('Mean contigs length(kb)'+'\t')
           for contig in mean_contig:
               outfile.write(str(round(contig/1000,1))+'\t')
           outfile.write('\n')
           outfile.write('Maximum contigs length(kb)'+'\t')
           for contig in max_contig:
               outfile.write(str(round(contig/1000,1))+'\t')
           outfile.write('\n')
           print('Calculationg Scaffold N50...')
           outfile.write('Scaffold N50(kb)'+'\t')
           for scaffold in scaffoldlist:
               Scaffold_NX,lengthlist_scaffold=CalculateNx_contig(scaffold,min_scaffold)
               No_scaffold.append(len(list(lengthlist_scaffold.values())))
               mean_scaffold.append(np.mean(list(lengthlist_scaffold.values())))
               max_scaffold.append(max(list(lengthlist_scaffold.values())))
               NGX_scaffold.append(CalculateNGX_contig(lengthlist_scaffold,singlen*1000000))
               outfile.write(str(round(Scaffold_NX[4],1))+'\t')
           outfile.write('\n')
           print('Calculationg Scaffold NG50...')
           outfile.write('Scaffold NG50(kb)'+'\t')
           for NGX in NGX_scaffold:
               if len(NGX)<5:
                   outfile.write('0'+'\t')
               else:
                   outfile.write(str(round(NGX[4],1))+'\t')
           outfile.write('\n')
           outfile.write('Number of scaffolds'+'\t')
           for scaffold in No_scaffold:
               outfile.write(str(scaffold)+'\t')
           outfile.write('\n')
           outfile.write('Mean scaffolds length(kb)'+'\t')
           for scaffold in mean_scaffold:
               outfile.write(str(round(scaffold/1000,1))+'\t')
           outfile.write('\n')
           outfile.write('Maximum scaffolds length(kb)'+'\t')
           for scaffold in max_scaffold:
               outfile.write(str(round(scaffold/1000,1))+'\t')
           outfile.write('\n')
           index=1
           #print('Reading refernece genome...')
           #(chrid,chrlen,singlen)=reflen_id(reference)
           list_align_contig=[]
           list_misassemble_contig=[]
           list_align_scaffold=[]
           list_align_block=[]
           print('Reading Alignment files...')
           for tsvid in tsvlist:
               align_contig=alignment_contig_coordinate(tsvid,chrid,chrlen)
               misassemble_contig=alignment_miss_contig(tsvid,chrid,chrlen)
               list_align_contig.append(align_contig)
               list_misassemble_contig.append(misassemble_contig)
           for infoid in range(len(infolist)):
               align_scaffold=analyze_scaffold_info(infolist[infoid],list_align_contig[infoid])
               align_block=analyze_block_info(infolist[infoid],list_align_contig[infoid],lengthlist_contig)
               list_align_scaffold.append(align_scaffold)
               list_align_block.append(align_block)
           print('Calculating Contig N50 (Just keep segments can be aligned to the reference genome)')
           outfile.write('Contig N50 (alignment,kb)'+'\t')
           index=0
           for contigid in list_align_contig:
               NX_contig_align=CalculateNx_contig_align(contigid)
               outfile.write(str(round(NX_contig_align[4],1))+'\t')
           outfile.write('\n')
           outfile.write('Scaffold N50 (alignment,kb)'+'\t')
           print('Calculating Scaffold N50 (Just keep segments can be aligned to the reference genome)')
           for scaffoldid in list_align_scaffold:
               NX_scaffold_align=CalculateNx_scaffold_align(scaffoldid)
               outfile.write(str(round(NX_scaffold_align[4]*1000,1))+'\t')
           outfile.write('\n')
           print('Calculate Contig NA50...')
           outfile.write('Contig NA50(kb)'+'\t')
           for misassemble_contig in list_misassemble_contig:
               NA_contig=CalculateNAx_contig(misassemble_contig)
               outfile.write(str(round(NA_contig[4],1))+'\t')
           outfile.write('\n')
           print('Calculate Scaffold NA50...')
           outfile.write('Scaffold NA50(kb)'+'\t')
           for align_block in list_align_block:
               NX_block=CalculateNx_block(align_block)
               outfile.write(str(round(NX_block[4]*1000,1))+'\t')
           outfile.write('\n')
           coverage_list=[]
           print('Calculate Contig NC50...')
           outfile.write('Contig NC50(kb)'+'\t')
           for contigid in list_align_contig:
               [NC_contig,list_NC_contig,coverage]=CalculateNCx_contig(contigid,singlen)
               coverage_list.append(coverage)
               if len(NC_contig)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_contig[4],1)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold NC50...')
           outfile.write('Scaffold NC50(kb)'+'\t')
           for scaffoldid in list_align_scaffold:
               [NC_scaffold,list_NC_scaffold]=CalculateNCx_scaffold(scaffoldid,singlen)
               if len(NC_scaffold)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_scaffold[4]*1000,1)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Contig NCA50...')
           outfile.write('Contig NCA50(kb)'+'\t')
           for misassemble_contig in list_misassemble_contig:
               [NCA_contig,list_NCA_contig,coverage]=CalculateNCx_contig(misassemble_contig,singlen)
               if len(NCA_contig)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NCA_contig[4],1)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold NCA50...')
           outfile.write('Scaffold NCA50(kb)'+'\t')
           for align_block in list_align_block: 
               [NC_block,list_NC_block]=CalculateNCx_block(align_block,singlen)
               if len(NC_block)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_block[4]*1000,1)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold genome coverage...')
           outfile.write('Genome Coverage(%)'+'\t')
           for coverage in coverage_list:
               outfile.write(str(round(coverage*100,2))+'%\t')
           outfile.close()
    return None
if __name__=="__main__":
    import sys
    if len(sys.argv) == 1:
         helpinfo()
    else:
        from collections import defaultdict
        import re
        import numpy as np
        import operator
        from scipy import stats as sta
        import math
        from itertools import groupby
        import gzip
        import getopt
        main()

