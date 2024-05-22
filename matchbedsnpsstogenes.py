#requires interlap
#pip3 install interlap
#python3 matchpeakstogenes.py

from interlap import InterLap
import pandas as pd

#size of distance outside the peak to look
windowsize = 10000

#extract name from gff defline
def def_parse(astr):
    dd = {}
    a = astr.split(';')
    for b in a:
        c = b.split('=')
        if len(c) < 2: continue
        dd[c[0]] = c[1]
    return dd

chr_tree = {}

fh = open("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3") # file downloaded from https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/
for x in fh:
    if x[0] == '#': continue
    y = x.strip().split('\t')
    if len(y) < 2: continue
    if y[2] != 'gene': continue
    mychr = y[0]
    mystart = int(y[3])
    mystop = int(y[4])
    dd = def_parse(y[-1])
    myname = dd['ID']
    if not mychr in chr_tree: chr_tree[mychr] = InterLap()
    chr_tree[mychr].add((mystart,mystop,{"name":myname,"strand":y[6]}))

fh = open("WW-MM_bQTL.bed") # file containing bQTL
for x in fh:
    y = x.strip().split('\t')
    mychr = y[0].split('-')[1]
    if not mychr in chr_tree:
        continue
    gene_dist_dict = {}
    upstream_genes = chr_tree[mychr].find((int(y[1])-windowsize,int(y[1])))
    for g in upstream_genes:
        mystart,mystop,myinfodict = g
        if not myinfodict["strand"] == '-': continue
        gene_dist_dict[myinfodict["name"]] = int(y[1]) - mystop
        if gene_dist_dict[myinfodict["name"]] < 0:
            gene_dist_dict[myinfodict["name"]] = 0
    downstream_genes = chr_tree[mychr].find((int(y[1]),int(y[1])+windowsize))
    for g in downstream_genes:
        mystart,mystop,myinfodict = g
        if not myinfodict["strand"] == '+': continue
        gene_dist_dict[myinfodict["name"]] = mystart - int(y[1])
        if gene_dist_dict[myinfodict["name"]] < 0:
            gene_dist_dict[myinfodict["name"]] = 0
    plist = y[:]
    if len(gene_dist_dict) == 0:
        plist.extend(["intergenic","intergenic"])
    else:
        gene_list = list(gene_dist_dict)
        gene_list.sort(key=lambda a:gene_dist_dict[a])
        if gene_dist_dict[gene_list[0]] == 0:
            plist.append("intragenic")
        else:
            plist.append(str(gene_dist_dict[gene_list[0]]))
        plist.append(gene_list[0])

    print("\t".join(plist))
