#! /usr/bin/env python

import subprocess
import os
import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description='Parse vgw output using gmap')
parser.add_argument('genome', help='vgw output file')
parser.add_argument('transcripts', help='transcript sequences')
parser.add_argument('-o','--overwrite', action='store_true', help='overwrite gmap database')
parser.add_argument('-g','--gff', nargs='?', default='parsing.gff', help='intermediate gff filename, default (%(default)s)')
parser.add_argument('--docker_vol', nargs='?', default='/data', help='shared docker volume, default (%(default)s)')

args = parser.parse_args()

args.gff = args.docker_vol + "/" + args.gff

genomename = args.genome[len(args.docker_vol):]
if genomename.startswith("/"):
    genomename = genomename[1:]

if not os.path.exists(args.docker_vol + "/gmapdb/"):
    subprocess.call("mkdir "+args.docker_vol+"/gmapdb", shell=True)
    args.overwrite = True

if args.overwrite or not os.path.exists(args.docker_vol +"/gmapdb/"+genomename):
    subprocess.call("gmap_build -d "+genomename+" -D "+args.docker_vol+"/gmapdb "+args.genome, shell=True)

gmapres = dict()

subprocess.call("gmap -d "+genomename+" -D "+args.docker_vol+"/gmapdb -t 10 -f 3 "+args.transcripts+" > "+args.gff, shell=True)

with open(args.gff, 'r') as ins:
    for l in ins:
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
ins.close()

# subprocess.call("makeblastdb -in "+args.genome+" -dbtype nucl", shell=True)
#
# blastres = dict()
#
# p2 = subprocess.Popen("blastn -db "+args.genome+" -query "+args.genome+" -outfmt 6 | sort -n -k 7", shell=True, stdout=subprocess.PIPE)
#
# for l in iter(p2.stdout.readline,''):
#     l = l.rstrip()
#     data = l.split("\t")
#     if data[0] != data[1] or data[2] != "100.00":
#         continue
#     #print data
#     if data[0] not in blastres:
#         blastres[data[0]] = []
#     blastres[data[0]].append(int(data[6]))

chunkbreak = "N"*500
#print chunkbreak

outfile = args.genome+"_gmapparsed"
out = open(outfile, 'w')

seqlens = []

for record in SeqIO.parse(args.genome, "fasta"):
    #print record.name
    newseq = []
    startlist = []
    endlist = []

    if record.name not in gmapres:
        continue

    startlist.append(0)
    for m in re.finditer(chunkbreak, str(record.seq)):
        #print "Break: %s - %s" % (m.start(), m.end())
        startlist.append(m.end())
        endlist.append(m.start())

    endlist.append(len(record.seq))
    #print record.name
    #print gmapres[record.name]
    #print blastres[record.name]
    for a in range(len(startlist)):
        #print "%s: %s - %s" % (a, startlist[a], endlist[a])
        # cstart = blastres[record.name][a]
        # nextstart = blastres[record.name][a+1]
        #print str(cstart)+" - "+str(nextstart)
        for g in gmapres[record.name]:
            #print "\t"+str(g)
            if g >= startlist[a] and g < endlist[a]:
                #print "This is a good chunk!"
                newseq.append(str(record.seq)[startlist[a]:endlist[a]])
                break

    # laststart = blastres[record.name][-1]
    # #print laststart
    # for g in gmapres[record.name]:
    #     if g >= laststart:
    #         #print "Last chunk is good!"
    #         newseq.append(str(record.seq[laststart-1:len(record.seq)]))
    #         break
    #print newseq
    seqlens.append(len("".join(newseq)))
    out.write(">"+record.name+"\n")
    out.write(("N"*500).join(newseq)+"\n")
    #break

seqlens = sorted(seqlens)
print "Num: "+str(len(seqlens))
print "Max: "+str(max(seqlens))
print "Min: "+str(min(seqlens))
print "Total: "+str(sum(seqlens))
half = float(sum(seqlens))/2.0
cum = 0
for l in seqlens:
    cum += l
    if cum >= half:
        print "N50: "+str(l)
        break
