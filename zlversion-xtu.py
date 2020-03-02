import pybedtools
from pyfasta import Fasta
import sys
from bx.intervals.intersection import Intersecter, Interval

##########################################################################################
################## MAKE SURE THE ORDER FOR ALL OF SPECIES SHOULD BE SAME! ################
##########################################################################################
path = sys.argv[1]

zmgff = pybedtools.BedTool("{0}/Zmays_493_RefGen_V4.gene.gff3".format(path))
sigff = pybedtools.BedTool("{0}/Sitalica_312_v2.2.gene.gff3".format(path))
sbgff = pybedtools.BedTool("{0}/Sbicolor_313_v3.1.gene.gff3".format(path))

zmfasta = Fasta("{0}/Zmays_493_APGv4-good.fa".format(path))
sifasta = Fasta("{0}/Sitalica_312_v2-good.fa".format(path))
sbfasta = Fasta("{0}/Sbicolor-good_313_v3.0.fa".format(path))

fastas = [zmfasta, sifasta, sbfasta]
gffs = [zmgff, sigff, sbgff]

fh = open('/home/zliang/Documents/genomeinfo/CNS/sorghum3_maize4_intelligent.csv', 'r')
fh.readline()
sdict = []

sh = open(sys.argv[2],'r')

only = set([])

for line in sh:
    new = line.strip()
    only.add(new)
sh.close()

syn_zm = set([])
syn_si = set([])
syn_sb = set([])

temp = open('temp-cand.csv','w')

for line in fh:
    nlist = []
    new = line.strip().split(',')
    if new[3] not in only and new[4] not in only:continue
    if new[0].startswith('N'):
        continue
    if new[3].startswith('N'):
        continue
    if new[6].startswith('N'):
        continue
    nlist.append(new[3])
    nlist.append(new[6])
    nlist.append(new[0])
    syn_zm.add(new[3])
    syn_si.add(new[6])
    syn_sb.add(new[0])
    temp.write(' '.join(nlist)+'\n')

temp.close()

syn = [syn_zm, syn_si, syn_sb]

syn_region = {}
mdict = {}
for i in gffs:
    k = gffs.index(i)
    syn_region[k] = {}
    for y in i:
        if y[2] != 'gene':
            continue
        if k == 0:
            m = y.name
            mygene = m.split('.')[0]
            if 'Zm' not in mygene:
                continue
            myg = mygene
        elif k == 1:
            m = y.name
            m1 = m.split('.')
            mygene = m1[0]+'.'+m1[1]
            myg = mygene
        elif k == 2:
            m = y.name
            m1 = m.split('.')
            mygene = m1[0]+'.'+m1[1]
            myg = mygene
        if myg not in syn[k]:
            continue
        if y[0] not in syn_region[k]:
            syn_region[k][y[0]] = Intersecter()
        syn_region[k][y[0]].add_interval(Interval(int(y[3]), int(y[4]), myg))
        if myg not in mdict:
            mdict[myg] = []
        mdict[myg].append(y.chrom)
        mdict[myg].append(int(y[3]))
        mdict[myg].append(int(y[4]))
        mdict[myg].append(y.strand)
temp.close()

cand = open('temp-cand.csv','r')
count = 0
for line in cand:
    count += 1
    new = line.strip().split(' ')
    myanns = []
    myseqs = []
    string = []
    for i in new:
        mygene = i
        mychr = mdict[i][0]
        mystart = mdict[i][1]
        mystop = mdict[i][2]
        mystrand = mdict[i][3]
        k = new.index(i)
        a = syn_region[k][mychr].find(mystart - 2000, mystart - 1) # finding intersection of 2kb region with other genes
        b = syn_region[k][mychr].find(mystop + 1, mystop + 2000)
        st = []
        sp = []
        if len(a) > 0: # if there are genes fallen into 2kb upstream region
            for c in a:
                sp.append(c.end)
            m = max(sp) + 1
            myst = m # using the end position to represent the start position of genes 
        else:
            myst = mystart - 2000
            if myst < 0:
                myst = 1
        if len(b) > 0:
            for c in b:
                st.append(c.start)
            n = min(st) - 1
            mysp = n
        else:
            mysp = mystop + 2000
        string = [mygene, mychr, mystrand, str(myst), str(mystart), str(mystop), str(mysp)]
        myanns.append("> {0}".format(' '.join(string)))
        myseqs.append(">{0}\n{1}".format(mygene, fastas[k].sequence({"chr": mychr, 'start': myst, 'stop': mysp})))
    out = open("outputs/gene_extract_{0}.fasta".format(str(count).zfill(5)), 'w')
    for q in myanns:
        out.write(q + '\n')
    for d in myseqs:
        out.write(d + '\n')
    out.close()
fh.close()
cand.close()