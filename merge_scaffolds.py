#! /usr/bin/env python

# Sort by longest VGW scaffold, then go through vgw.gff looking for other
# transcripts with a highish coverage that map to overlapping exon positions.
# Check flanking genome sequence matches, then assemble together with cap3.
# Re order if necessary/add in other info to produce new vgw scaffold with
# information from multiple transcripts

# Requires CAP3, BLAST and MAFFT


import os
import re
import argparse
import subprocess
import tempfile
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo

parser = argparse.ArgumentParser(description="Find duplicate genome scaffolds and merge")
parser.add_argument("genome", help="Parsed VGW scaffolds")
parser.add_argument("gff", help="GFF file with transcript mapped to VGW scaffolds")
parser.add_argument("transcripts", help="fasta file of transcript sequences")
parser.add_argument('-o', '--outfile', nargs='?', default = "merged.fasta", help="Outfile name (default %(default)s)")


args = parser.parse_args()

genomeorder = []
genomelens = dict()
genomeseq = dict()

for record in SeqIO.parse(args.genome, 'fasta'):
    if len(str(record.seq)) not in genomelens:
        genomelens[len(str(record.seq))] = []
    genomelens[len(str(record.seq))].append(record.name)
    genomeseq[record.name] = str(record.seq)


transseq = dict()

for record in SeqIO.parse(args.transcripts, 'fasta'):
    transseq[record.name] = str(record.seq)


for key, value in sorted(genomelens.iteritems()):
    #print key, value
    for v in value:
        genomeorder.append(v)

fout = open(args.outfile, "w")

alldone = []

for genome in reversed(genomeorder):

    print genome

    if genome in alldone:
        continue

    selflocs = dict()
    overlaps = dict()

    alldone.append(genome)

    try:
        res = subprocess.check_output("grep '^"+genome+"[[:space:]]' "+args.gff, shell=True)
    except:
        print "No results in gff file"
        continue

    interesting = 0
    skipgene = []

    for r in res.split("\n"):
        if len(r) < 2:
            continue
        #print r
        data = r.split("\t")
        if data[2] == "gene":
            continue
        trans = data[8][3:]
        trans = re.sub("\.mrna.*","",trans)
        #print trans
        if trans == genome and data[2] == "exon":
            selflocs[int(data[3])] = int(data[4])

    for r in res.split("\n"):
        if len(r) < 2:
            continue
        #print r
        data = r.split("\t")
        if data[2] == "gene":
            continue
        trans = data[8][3:]
        trans = re.sub("\.mrna.*","",trans)

        if trans != genome and trans not in skipgene:
            print r
            coverage = None
            if data[2] == "mRNA":
                coverage = re.sub(".*coverage=", "", data[8])
                coverage = re.sub(";identity.*", "", coverage)
                print coverage
                if float(coverage) < 20:
                    skipgene.append(trans)
            if data[2] == "exon":
                for s,e in sorted(selflocs.iteritems()):
                    #print s, e
                    if int(data[3]) < e and int(data[4]) > s:
                        if trans not in overlaps:
                            overlaps[trans] = dict()
                        overlaps[trans][s] = e
                        print "Overlaps!"

            #print r

    if len(overlaps) > 0:
        print overlaps
        merge = []
        winning = 0
        for t in overlaps:

            if t not in genomeseq:
                alldone.append(t)
                print "No VGW scaffold to merge with"
                continue

            if t in alldone:
                print "Already merged this scaffold elsewhere"
                continue

            tmp1 = tempfile.NamedTemporaryFile(delete=False)
            tmp1.write(">"+genome+"\n"+genomeseq[genome]+"\n")
            tmp1.close()
            tmp2 = tempfile.NamedTemporaryFile(delete = False)
            tmp2.write(">"+t+"\n"+genomeseq[t]+"\n")
            tmp2.close()

            hits = dict()
            lesshits = dict()

            res  = subprocess.check_output("blastn -query "+tmp1.name+" -subject "+tmp2.name+" -outfmt 6 ", shell=True)
            for r in res.split("\n"):
                if len(r) < 2:
                    continue
                print r
                data = r.split("\t")
                if float(data[2]) < 99:
                    if float(data[2]) > 95:
                        lesshits[int(data[6])] = int(data[7])
                    continue
                hits[int(data[6])] = int(data[7])
                #print r

            fail = 0

            for s, e in overlaps[t].iteritems():
                print s, e
                gotit = 0
                for gs, ge in hits.iteritems():
                    if s > gs and e < ge:
                        print "In BLAST!"
                        gotit += 1
                        break
                if gotit == 0:
                    lessgotit = 0
                    for gs, ge in lesshits.iteritems():
                        if s > gs and e < ge:
                            print "In lower identity BLAST!"
                            lessgotit += 1
                            break
                    if lessgotit > 0:
                        fail += 1
                    else:
                        fail += len(overlaps[t])

            os.unlink(tmp1.name)
            os.unlink(tmp2.name)

            if fail >= len(overlaps[t]):
                print "They all failed!!"
            else:
                merge.append(t)

        if len(merge) == 0:
            print "Nothing to merge"
            fout.write(">"+genome+"\n"+genomeseq[genome]+"\n")
            continue


        tmpname = "tmpmerge.fa"
        tmp = open(tmpname, 'w')
        contigs = genomeseq[genome].split("N"*500)
        for c in range(len(contigs)):
            tmp.write(">"+genome+"_contig"+str(c+1)+"\n"+contigs[c]+"\n")
        for t in merge:
            contigs = genomeseq[t].split("N"*500)
            for c in range(len(contigs)):
                tmp.write(">"+t+"_contig"+str(c+1)+"\n"+contigs[c]+"\n")
        tmp.close()

        subprocess.call("cap3 "+tmpname+" -f 300 -h 50 > "+tmpname+"_info ; cat "+tmpname+".cap.singlets "+tmpname+".cap.contigs > "+tmpname+"_cap3", shell=True)

        cap3seq = dict()
        extra = 0
        extraids = []

        for record in SeqIO.parse(tmpname+"_cap3", "fasta"):
            cap3seq[record.name] = record.seq
            #print record.name
            if genome+"_contig" not in record.name and not record.name.startswith("Contig"):
                extra += 1
                extraids.append(record.name)
                print record.name

        if extra > 0:
            #print "WILL CHECK EXTRA CANNOT BE COMBINED WITH BLAST"
            print extraids

            subprocess.call("makeblastdb -dbtype nucl -in "+tmpname+"_cap3 -parse_seqids", shell=True)

            groups = dict()
            remove = []
            beenreversed = []

            for i in extraids:

                tmp1 = tempfile.NamedTemporaryFile(delete = False)
                tmp1.write(">"+i+"\n"+str(cap3seq[i])+"\n")
                tmp1.close()

                done = []

                blastres = subprocess.check_output("blastn -db "+tmpname+"_cap3 -query "+tmp1.name+" -outfmt 6", shell=True)
                for r in blastres.split("\n"):
                    if len(r) == 0:
                        continue
                    data = r.split("\t")
                    if data[0] == data[1]:
                        continue
                    gene1 = re.sub("_contig.*", "", data[0])
                    gene2 = re.sub("_contig.*", "", data[1])
                    if gene1 == gene2 or float(data[2]) < 95:
                        continue
                    if data[0] in done:
                        continue
                    print r
                    if int(data[8]) > int(data[9]):
                        print "Need to revcom "+data[0]
                        if data[0] not in beenreversed:
                            cap3seq[data[0]] = cap3seq[data[0]].reverse_complement()
                            beenreversed.append(data[0])
                    # elif gene2 == genome and int(data[8]) > int(data[7]):
                    #     print "Need to revcom "+data[0]
                    #     if data[0] not in beenreversed:
                    #         cap3seq[data[0]] = cap3seq[data[0]].reverse_complement()
                    #         beenreversed.append(data[0])
                    group = 0
                    for g in groups:
                        for c in groups[g]:
                            if c == data[0] or c == data[1]:
                                group = g
                                break
                    if group == 0:
                        group = len(groups)+1
                    if group not in groups:
                        groups[group] = []
                    if data[0] not in groups[group]:
                        groups[group].append(data[0])
                    if data[1] not in groups[group]:
                        groups[group].append(data[1])
                    done.append(data[0])

                os.unlink(tmp1.name)

            print groups

            for g in groups:

                tmp1 = tempfile.NamedTemporaryFile(delete = False)
                for i in groups[g]:
                    print i
                    tmp1.write(">"+i+"\n"+str(cap3seq[i])+"\n")

                tmp1.close()

                p1 = subprocess.Popen("mafft --auto "+tmp1.name, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

                #muscle_cline = MuscleCommandline(input=tmpfile)
                stdout, stderr = p1.communicate()
                align = AlignIO.read(StringIO(stdout), "fasta")
                #print str(align)
                summary_align = AlignInfo.SummaryInfo(align)
                consensus = summary_align.dumb_consensus(ambiguous='N')
                #print str(consensus).upper()

                os.unlink(tmp1.name)

                n = 0
                for a in str(consensus).upper():
                    if a == "N":
                        n += 1
                pern = (float(n)/float(len(str(consensus))))*100
                print "Percent N: "+str(pern)+"%"

                if pern > 5:
                    continue

                cap3out = open(tmpname+"_cap3", 'a')
                cap3out.write('>consensus_contig'+str(g)+"\n"+str(consensus).upper()+"\n")
                cap3out.close()

                cap3seq["consensus_contig"+str(g)] = consensus

                for i in groups[g]:
                    remove.append(i)
                    if i in extraids:
                        extraids.remove(i)

            #print extraids

            cap3out2 = open(tmpname+"_cap32", 'w')
            for record in SeqIO.parse(tmpname+"_cap3", "fasta"):
                if record.name in remove:
                    continue
                cap3out2.write(">"+record.name+"\n"+str(record.seq)+"\n")

            cap3out2.close()
            #break
        else:
            subprocess.call("cp "+tmpname+"_cap3 "+tmpname+"_cap32", shell=True)

        subprocess.call("makeblastdb -dbtype nucl -in "+tmpname+"_cap32 -parse_seqids", shell=True)

        tmp = tempfile.NamedTemporaryFile(delete = False)
        tmp.write(">"+genome+"\n"+transseq[genome]+"\n")
        tmp.close()

        beenreversed = []

        orderdict = dict()
        p1 = subprocess.Popen('blastn -db '+tmpname+'_cap32 -query '+tmp.name+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
        for l in iter(p1.stdout.readline,''):
            l = l.rstrip()
            #print l
            data = l.split("\t")
            if int(data[8]) > int(data[9]):
                print "Need to reverse compliment "+data[1]
                if data[1] not in beenreversed:
                    cap3seq[data[1]] = cap3seq[data[1]].reverse_complement()
                    beenreversed.append(data[1])
            if int(data[6]) not in orderdict:
                orderdict[int(data[6])] = dict()
            if data[1] in beenreversed:
                orderdict[int(data[6])][int(data[9])] = data[1]
            else:
                orderdict[int(data[6])][int(data[8])] = data[1]

        os.unlink(tmp.name)

        laststart = 0

        order = []
        for start in sorted(orderdict):
            for qstart in sorted(orderdict[start], reverse=True):
                print orderdict[start][qstart]+"\t"+str(start)+"\t"+str(qstart)
                if orderdict[start][qstart] not in order and start != laststart:
                    order.append(orderdict[start][qstart])
                laststart = start

        print order

        if extra > 0:

            extrats = []

            for i in extraids:

                if i in order:
                    print "Already ordered!"
                    continue

                if i in remove:
                    print "Already merged"
                    continue

                t = re.sub("_contig.*", "", i)
                if t not in extrats:
                    extrats.append(t)

            for t in extrats:

                tmp = tempfile.NamedTemporaryFile(delete = False)
                tmp.write(">"+t+"\n"+transseq[t]+"\n")
                tmp.close()

                orderdict = dict()
                p1 = subprocess.Popen('blastn -db '+tmpname+'_cap32 -query '+tmp.name+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
                for l in iter(p1.stdout.readline,''):
                    l = l.rstrip()
                    print l
                    data = l.split("\t")
                    if int(data[6]) not in orderdict:
                        orderdict[int(data[6])] = dict()
                    orderdict[int(data[6])][int(data[8])] = data[1]

                prevpos = -1
                inserted = 0

                laststart = 0

                for start in sorted(orderdict):
                    for qstart in sorted(orderdict[start], reverse=True):
                        print orderdict[start][qstart]+"\t"+str(start)+"\t"+str(qstart)
                        if orderdict[start][qstart] not in order:
                            print "New contig"
                            if start != laststart:
                                print "Adding"
                                order.insert(prevpos+1, orderdict[start][qstart])
                        else:
                            prevpos = order.index(orderdict[start][qstart])
                        laststart = start

                print order

                os.unlink(tmp.name)

            #break

        finalseq = []

        for c in order:
            finalseq.append(str(cap3seq[c]))

        bridge = "N"*500

        fout.write(">"+genome+" "+" ".join(merge)+"\n"+bridge.join(finalseq)+"\n")

        for t in merge:
            alldone.append(t)

        #break

    else:

        fout.write(">"+genome+"\n"+genomeseq[genome]+"\n")

fout.close()
