#! /usr/bin/python

import subprocess
import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Parse vgw output using gmap')
parser.add_argument('genome', help='vgw output file')
parser.add_argument('transcripts', help='transcript sequences')
parser.add_argument('-o','--overwrite', action='store_true', help='overwrite gmap database')

args = parser.parse_args()

if not os.path.exists("gmapdb/"):
    subprocess.call("mkdir gmapdb")
    args.overwrite = True

if args.overwrite or not os.path.exists("gmapdb/"+args.genome):
    subprocess.call("gmap_build -d "+args.genome+" -D ./gmapdb "+args.genome, shell=True)

gmapres = dict()

p1 = subprocess.Popen("gmap -d "+args.genome+" -D ./gmapdb -f 3 "+args.transcripts, shell=True, stdout=subprocess.PIPE)

for l in iter(p1.stdout.readline,''):
    l = l.rstrip()
    data = l.split("\t")
    if len(data) < 4:
        continue
    if "Target="+data[0]+" " not in data[8]:
        continue
    #print data
    if data[0] not in gmapres:
        gmapres[data[0]] = []
    gmapres[data[0]].append(int(data[3]))


subprocess.call("makeblastdb -in "+args.genome+" -dbtype nucl", shell=True)

blastres = dict()

p2 = subprocess.Popen("blastn -db "+args.genome+" -query "+args.genome+" -outfmt 6 | sort -n -k 7", shell=True, stdout=subprocess.PIPE)

for l in iter(p2.stdout.readline,''):
    l = l.rstrip()
    data = l.split("\t")
    if data[0] != data[1] or data[2] != "100.00":
        continue
    #print data
    if data[0] not in blastres:
        blastres[data[0]] = []
    blastres[data[0]].append(int(data[6]))

outfile = args.genome+"_gmapparsed"
out = open(outfile, 'w')

for record in SeqIO.parse(args.genome, "fasta"):
    newseq = []
    if record.name not in gmapres:
        continue
    #print record.name
    #print gmapres[record.name]
    #print blastres[record.name]
    for a in range(len(blastres[record.name])-1):
        cstart = blastres[record.name][a]
        nextstart = blastres[record.name][a+1]
        #print str(cstart)+" - "+str(nextstart)
        for g in gmapres[record.name]:
            #print "\t"+str(g)
            if g >= cstart and g < nextstart:
                #print "This is a good chunk!"
                newseq.append(str(record.seq[cstart-1:nextstart-501]))
                break

    laststart = blastres[record.name][-1]
    #print laststart
    for g in gmapres[record.name]:
        if g >= laststart:
            #print "Last chunk is good!"
            newseq.append(str(record.seq[laststart-1:len(record.seq)]))
            break
    out.write(">"+record.name+"\n")
    out.write(("N"*500).join(newseq)+"\n")
