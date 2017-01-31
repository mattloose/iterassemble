#! /usr/bin/python

import subprocess
import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Summarise vgw results using gmap')
parser.add_argument('genome', help='vgw output file')
parser.add_argument('transcripts', help='transcript sequences')
parser.add_argument('-o','--overwrite', action='store_true', help='overwrite gmap database')
parser.add_argument('-g','--gff3', nargs='?', default='vgw.gff', help='output gff3 filename, default (%(default)s)')
parser.add_argument('-s','--summary', nargs='?', default='vgw_summary.txt', help='output file name, default (%(default)s)')

args = parser.parse_args()

if not os.path.exists("gmapdb/"):
    subprocess.call("mkdir gmapdb")
    args.overwrite = True

if args.overwrite or not os.path.exists("gmapdb/"+args.genome):
    subprocess.call("gmap_build -d "+args.genome+" -D ./gmapdb "+args.genome, shell=True)

genomehash = dict()

for record in SeqIO.parse(args.genome, 'fasta'):
    genomehash[record.name] = record.seq

transhash = dict()

for record in SeqIO.parse(args.transcripts, 'fasta'):
    transhash[record.name] = record.seq

subprocess.call("gmap -d "+args.genome+" -D ./gmapdb -f 2 -z sense_force "+args.transcripts+" > "+args.gff3, shell=True)

gene = ""
path = ""
mrna = ""
newid = ""
last = []

resdict = dict()

with open(args.gff3, 'r') as ins:
    for l in ins:
        if l.startswith('#'):
            continue
        data = l.rstrip().split("\t")
        #print data
        if data[2] == 'mRNA':
            if len(last) > 1:
                #print "Last exon!"
                info = last[8].split(";")
                coords = info[3].split(" ")
                #print "Genome len: "+str(len(genomehash[gene]))
                resdict[newid]['3prime'] = len(genomehash[last[0]]) - int(last[4])
                #print "Trans len: "+str(len(transhash[gene]))
                if coords[3] == "+" and int(coords[2]) < len(transhash[gene]):
                    extra = len(transhash[gene]) - int(coords[2])
                    #print extra
                    resdict[newid]['bp_notmapped_ext'] += extra
                last = []
            info = data[8].split(";")
            gene = info[1].split("=")[1]
            path = info[2].split(".")[1]
            mrna = info[0].split(".")[1]
            #print "Gene: %s \t Path: %s \t mRNA: %s" % (gene, path, mrna)
            newid = gene+":"+path+":"+mrna
            if newid not in resdict:
                resdict[newid] = dict()
                resdict[newid]['exons'] = 0
                resdict[newid]['bp_notmapped_int'] = 0
                resdict[newid]['bp_notmapped_ext'] = 0
                resdict[newid]['introns_bridged'] = 0
                resdict[newid]['introns'] = []
                resdict[newid]['bridged'] = []
                resdict[newid]['5prime'] = 0
                resdict[newid]['3prime'] = 0
                resdict[newid]['genome'] = data[0]
        if data[2] == 'exon':
            resdict[newid]['exons'] += 1
            if 'exon1;' in data[8]:
                #print "First exon!"
                resdict[newid]['5prime'] = int(data[3]) - 1
                info = data[8].split(";")
                coords = info[3].split(" ")
                if coords[3] == "+":
                    resdict[newid]['bp_notmapped_ext'] += (int(coords[1])-1)
            else:
                previnfo = last[8].split(";")
                prevcoords = previnfo[3].split(" ")
                info = data[8].split(";")
                coords = info[3].split(" ")
                #print "Trans: %s - %s" % (prevcoords[2], coords[1])
                if int(prevcoords[2])+1 < int(coords[1]):
                    extra = int(coords[1]) - int(prevcoords[2]) - 1;
                    resdict[newid]['bp_notmapped_int'] += extra

                #print "Genome: %s - %s" % (last[4], data[3])
                length = int(data[3]) - int(last[4])
                resdict[newid]['introns'].append(length)
                intronseq = genomehash[data[0]][int(last[4]):int(data[3])]
                #print intronseq
                if not "N"*500 in intronseq:
                    #print "BRIDGED"
                    resdict[newid]['introns_bridged'] += 1
                    resdict[newid]['bridged'].append(length)
            last = data
        if data[2] == 'CDS' and len(last) > 1:
            #print "Last exon!"
            info = last[8].split(";")
            coords = info[3].split(" ")
            #print "Genome len: "+str(len(genomehash[gene]))
            resdict[newid]['3prime'] = len(genomehash[last[0]]) - int(last[4])
            #print "Trans len: "+str(len(transhash[gene]))
            if coords[3] == "+" and int(coords[2]) < len(transhash[gene]):
                extra = len(transhash[gene]) - int(coords[2])
                #print extra
                resdict[newid]['bp_notmapped_ext'] += extra
            last = []

ins.close()

with open(args.summary, 'w') as ins:
    ins.write("Genome\tGene\tPath\tmRNA\tTranscript.length\tvgw.length\tvgw.minus500.length\tNo.Exons\t5prime.bp\t3prime.bp\tunmapped.internal.bp\tunmapped.external.bp\tmax.intron\tmin.intron\tmean.intron\tno.bridged\tmin.bridged\tmax.bridged\tmean.bridged\n")
    for i in resdict:
        names = i.split(":")
        ins.write(resdict[i]['genome']+"\t")
        ins.write("\t".join(names)+"\t")
        ins.write("%s\t%s\t%s\t" % (len(transhash[names[0]]), len(genomehash[resdict[i]['genome']]), len(genomehash[resdict[i]['genome']]) - (500 * ((resdict[i]['exons'] - 1) - resdict[i]['introns_bridged']))))
        ins.write("%s\t%s\t%s\t%s\t%s\t" % (resdict[i]['exons'], resdict[i]['5prime'], resdict[i]['3prime'], resdict[i]['bp_notmapped_int'], resdict[i]['bp_notmapped_ext']))
        if resdict[i]['exons'] > 1:
            introns = resdict[i]['introns']
            introns = sorted(introns)
            ins.write("%s\t%s\t%s\t" % (max(introns), min(introns), sum(introns)/len(introns)))
        else:
            ins.write("0\t0\t0\t")
        ins.write(str(resdict[i]['introns_bridged'])+"\t")
        if resdict[i]['introns_bridged'] > 0:
            bridged = resdict[i]['bridged']
            bridged = sorted(bridged)
            ins.write("%s\t%s\t%s\n" % (max(bridged), min(bridged), sum(bridged)/len(bridged)))
        else:
            ins.write("0\t0\t0\n")
