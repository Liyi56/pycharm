#!/usr/bin/env python
# -*- coding=utf-8 -*-
import math
promoters={}
enhancer={}
enhancer2={}
H2AZ={}
peak_H2AZ={}
enhancer_H2AZ={}
enhancer_gene={}
m3 = open('enhancer_flower.txt', 'r')
prev_chr=""
ii1=0
for line3 in m3:
    lm3 = line3.rstrip("\n").split("\t")
    chrm = lm3[0].lower()
    start=lm3[1]
    end=lm3[2]
    closed_gene=lm3[3]
    if ii1 == 0:
        ii1 += 1
        continue
    if chrm != prev_chr:
        enhancer[chrm] = {}
        enhancer_gene[chrm]={}
        prev_chr = chrm
    enhancer[chrm][start]=end
    enhancer_gene[chrm][start]=closed_gene
m3.close()

def overlap(a1,a2):
    for i in a2:
        if not float(a2[i])<float(start) and not float(i)>float(a1[start]):
            enhancer2[i] = a2[i]
            return True
        else:
            continue

f2=open('Table S1.txt')
i1=0
ii2=0
prev_chr2=""
for line1 in f2:
    lf2 = line1.rstrip("\n").split("\t")
    Chr=lf2[0].strip().lower()
    start=lf2[1].strip()
    end=lf2[2].strip()
    if ii2==0:
        ii2+=1
        continue
    if Chr != prev_chr2:
        H2AZ[Chr] = {}
        peak_H2AZ[Chr]={}
        prev_chr2 = Chr
    H2AZ[Chr][start]=end
    for i in enhancer[Chr]:
        if not float(enhancer[Chr][i])<=float(start) and not float(i)>=float(H2AZ[Chr][start]):
            i1 += 1
            enhancer2[i] = enhancer[Chr][i]
            peak_H2AZ[Chr][start] = end
            enhancer_H2AZ[i]=enhancer_gene[Chr][i]
            break #因为一个peak可能对应多个enhancer
            #return True
        else:
            continue
    # if overlap(H2AZ[Chr],enhancer[Chr]):#比较每个peak 和enhancer即{Chr1:{},Chr2:{},Chr3:{},Chr4:{},Chr5:{}}是否有overlap
    #     peak_H2AZ[start]=end
    #     i1+=1
print i1,len(peak_H2AZ['chr1'])+len(peak_H2AZ['chr2'])+len(peak_H2AZ['chr3'])+len(peak_H2AZ['chr4'])+len(peak_H2AZ['chr5']),peak_H2AZ
f2.close()
print len(enhancer2),enhancer2
enhancer3=enhancer2
enhancer2={}

peak_H3K27={}
H3K27={}
f2=open('Table S3.txt')
i2=0
ii3=0
prev_chr3=""
for line1 in f2:
    lf2 = line1.rstrip("\n").split("\t")
    Chr=lf2[0].strip('')
    start=lf2[1].strip('')
    end=lf2[2].strip('')
    if ii3==0:
        ii3+=1
        continue
    if Chr != prev_chr3:
        H3K27[Chr] = {}
        peak_H3K27[Chr]={}
        prev_chr3 = Chr
    H3K27[Chr][start]=end
    for i in enhancer[Chr]:
        if not float(enhancer[Chr][i])<=float(start) and not float(i)>=float(H3K27[Chr][start]):
            enhancer2[i] = enhancer[Chr][i]
            peak_H3K27[Chr][start] = end
            i2 += 1
            break
        else:
            continue
    # if overlap(H3K27[Chr],enhancer[Chr]):#比较每个peak 和enhancer即{Chr1:{},Chr2:{},Chr3:{},Chr4:{},Chr5:{}}是否有overlap
    #     i2+=1
    #     peak_H3K27[start]=end
print i2,len(peak_H3K27['chr1'])+len(peak_H3K27['chr2'])+len(peak_H3K27['chr3'])+len(peak_H3K27['chr4'])+len(peak_H3K27['chr5']),peak_H3K27
print len(enhancer2),enhancer2
print len(enhancer_H2AZ),enhancer_H2AZ
f2.close()
H2AZ_H3K27=[i for i in enhancer3.keys() if i in enhancer2.keys()]
print 100*(len(enhancer3)-len(H2AZ_H3K27))/len(enhancer3),len(enhancer3)-len(H2AZ_H3K27)
print len(H2AZ_H3K27)
print 100*(len(enhancer2)-len(H2AZ_H3K27))/len(enhancer2),len(enhancer2)-len(H2AZ_H3K27)

gene_id_fpkm={}
m3 = open('genes.fpkm_tracking', 'r')
for line3 in m3:
    lm3 = line3.rstrip("\n").split("\t")
    gene_id = lm3[0]
    fpkm = lm3[9]
    gene_id_fpkm[gene_id]=fpkm
m3.close()

m3 = open('enhancer_flower.txt', 'r')
m4=open('H2AZ_nonH2AZ_enhancer.txt','w')
ii1 = 0
for line3 in m3:
    lm3 = line3.rstrip("\n").split("\t")
    chrm = lm3[0].lower()
    start = lm3[1]
    end = lm3[2]
    closed_gene = lm3[3]
    if ii1 == 0:
        ii1 += 1
        m4.write(line3.strip()+'\t'+'Type'+'\t'+'group'+'\t'+'WT_FPKM'+'\t'+'log2(FPKM+1)'+'\n')
        continue
    if closed_gene not in gene_id_fpkm.keys():
        gene_id_fpkm[closed_gene]=0
    if start in enhancer3.keys() and end==enhancer3[start]:
        m4.write(line3.strip()+'\t'+'H2AZ'+'\t'+'1'+'\t'+str(gene_id_fpkm[closed_gene])+'\t'+str(math.log(float(gene_id_fpkm[closed_gene])+1,2))+'\n')
    else:
        m4.write(line3.strip()+'\t'+'non_H2AZ'+'\t'+'2'+'\t'+str(gene_id_fpkm[closed_gene])+'\t'+str(math.log(float(gene_id_fpkm[closed_gene])+1,2))+'\n')
m4.close()
m3.close()



m4=open('H2AZ_nonH2AZ_enhancer_H3K27me3.txt','w')
ii1 = 0
for i in enhancer2.keys():
    m3 = open('enhancer_flower.txt', 'r')
    for line3 in m3:
        lm3 = line3.rstrip("\n").split("\t")
        chrm = lm3[0].lower()
        start = lm3[1]
        end = lm3[2]
        closed_gene = lm3[3]
        if ii1 == 0:
            ii1 += 1
            m4.write(line3.strip()+'\t'+'Type'+'\t'+'group'+'\t'+'WT_FPKM'+'\t'+'log2(FPKM+1)'+'\n')
            continue
        if closed_gene not in gene_id_fpkm.keys():
            gene_id_fpkm[closed_gene]=0
        if i==start:
            if i in H2AZ_H3K27:
                m4.write(line3.strip() + '\t' + 'w/' + '\t' + '1' + '\t' + str(gene_id_fpkm[closed_gene]) + '\t' + str(math.log(float(gene_id_fpkm[closed_gene]) + 1, 2)) + '\n')
            else:
                m4.write(line3.strip()+'\t'+'w/o'+'\t'+'2'+'\t'+str(gene_id_fpkm[closed_gene])+'\t'+str(math.log(float(gene_id_fpkm[closed_gene])+1,2))+'\n')
        else:
            continue
    m3.close()
m4.close()
