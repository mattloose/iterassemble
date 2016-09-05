#! /usr/bin/python
import subprocess
import sys, getopt, argparse
import re
import os.path
import glob
import multiprocessing as mp
from Bio import SeqIO
import atexit
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo


def load_ids(file):

    idlist = []
    with open(file,'r') as ins:
        for line in ins:
            if re.match(">", line):
                line = re.sub("^>","",line)
                line = re.sub("\s.*$","",line)
                idlist.append(line)

    ins.close()
    return idlist


def assemble (i, id, arr1, args):
    dir = id + "_files"
    if not os.path.exists(dir):
        subprocess.call('mkdir '+dir, shell=True)

    f1 = dir + "/iter" + str(i) + "_R1.fastq"
    f2 = dir + "/iter" + str(i) + "_R2.fastq"
    fids = dir + "/iter" + str(i) + "_ids.txt"

    with open(fids, 'w') as ins:
        ins.write("\n".join(arr1))
    ins.close()

    if i > 1:
        subprocess.call("bwa index "+dir+"/iter"+str(i-1)+"_cap3_pass.fasta; bwa mem "+dir+"/iter"+str(i-1)+"_cap3_pass.fasta "+dir+"/iter"+str(i-1)+"_R1.fastq "+dir+"/iter"+str(i-1)+"_R2.fastq | samtools view -F 4 -f 8 - | awk '{print $1}' >> "+fids, shell=True)

    subprocess.call("ls "+args.d+"/seq*_R1.fastq | parallel -j 2 -k 'cat {} | fqextract "+fids+"' > "+f1, shell=True)
    subprocess.call("ls "+args.d+"/seq*_R2.fastq | parallel -j 2 -k 'cat {} | fqextract "+fids+"' > "+f2, shell=True)

    conf = dir + "/conf.txt"
    with open(conf, 'w') as ins:
        ins.write("max_rd_len=150\n[LIB]\navg_ins="+str(args.insert)+"\nreverse_seq=0\nasm_flags=3\nrd_len_cutoff=100\nrank=1\npair_num_cutoff=3\n")
        ins.write("q1="+f1+"\n")
        ins.write("q2="+f2+"\n")
    ins.close()

    soapout1 = dir + "/iter" + str(i) + "_63soap";
    soapout2 = dir + "/iter" + str(i) + "_31soap";
    if os.path.exists(soapout1+".scafSeq"):
        subprocess.call("rm "+soapout1+"* "+soapout2+"* ", shell=True)
    soapout = dir + "/iter" + str(i) + "_soap";

    subprocess.call('SOAPdenovo-63mer all -s '+conf+' -K 63 -R -F -p 2 -o '+soapout1, shell=True)

    subprocess.call('SOAPdenovo-63mer all -s '+conf+' -R -F -p 2 -o '+soapout2, shell=True)

    subprocess.call('cat ' + soapout1 + '.scafSeq ' + soapout2 + '.scafSeq > ' +  soapout + '.scafSeq', shell=True)

    if i > 1:
        #fa = dir + '/iter' + str(i) + '.fasta'
        #subprocess.call('cat '+f1+' '+f2+' | awk \'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}\' > '+fa+ ' ; cap3 '+fa, shell=True)
        subprocess.call('bwa mem '+dir+'/iter'+str(i-1)+'_cap3_pass.fasta '+f1+' '+f2+' > '+dir+'/iter'+str(i)+'.sam', shell=True)
        subprocess.call('grep ">" '+dir+'/iter'+str(i-1)+'_cap3_pass.fasta | sed \'s/^>//\' | while read -r line; do egrep "$line\s" '+dir+'/iter'+str(i)+'.sam | bam2fastx -a -A -s -Q -o '+dir+'/temp.fa - ; if [[ -s '+dir+'/temp.fa ]]; then timeout 20m cap3 '+dir+'/temp.fa; cat '+dir+'/temp.fa.cap.contigs >> ' + soapout + '.scafSeq; fi; done', shell=True)
        #subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta '+fa+'.cap.contigs '+fa+'.cap.singlets >> ' + soapout + '.scafSeq', shell=True)
        #subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta '+fa+'.cap.contigs >> ' + soapout + '.scafSeq', shell=True)
        subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta >> ' + soapout + '.scafSeq', shell=True)


    subprocess.call('cap3 ' + soapout + '.scafSeq -k 0 -p 75 -o 30 -h 80 -f 200 -g 4', shell=True)

    cap3 = dir + "/iter" + str(i) + "_cap3.fasta"
    subprocess.call('cat ' + soapout + '.scafSeq.cap.contigs ' + soapout + '.scafSeq.cap.singlets > ' + cap3, shell=True)

    passfile = dir + "/iter" + str(i) + "_cap3_pass.fasta"
    subprocess.call('bwa mem ' + args.cDNA + ' ' + cap3 + ' | grep "'+id+'"| bam2fastx -s -M -Q -a -o ' + passfile + ' - ', shell=True)

    #count = subprocess.check_output('ls -l ' + passfile + " | awk '{print $5}' ", shell=True)

    return id

def final_process (args, i, ID):
    dir = ID + "_files"
    print dir

    infile = dir +"/iter" + str(i) + "_cap3_pass.fasta"
    passfile = dir +"/iter" + str(i) + "_cap3_pass.fasta.renamed"

    sc = 0
    with open(passfile, 'w') as ins:
        for record in SeqIO.parse(infile, "fasta"):
            sc += 1
            ins.write(">Contig"+str(sc)+"\n"+str(record.seq)+"\n")
    ins.close()

    midfile = dir + "/final_all_seq.fasta"

    finalfile = dir +"/final_one_seq.fasta"

    seqhash = dict()
    for record in SeqIO.parse(passfile, "fasta"):
        seqhash[record.name] = record.seq

    subprocess.call("makeblastdb -dbtype nucl -in "+passfile, shell=True)

    revcom = dict()
    p1 = subprocess.Popen('blastn -db '+passfile+' -query '+args.cDNA+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        print l
        data = l.split("\t")
        if data[0] == ID:
            if (int(data[8]) > int(data[9])):
                #print "Reverse"
                revcom[data[1]] = 1

    for seqid in revcom:
        #print "Reversing "+seqid
        tmpseq = seqhash[seqid]
        seqhash[seqid] = tmpseq.reverse_complement()

    if len(seqhash) == 1:
        print "Only one sequence"
        with open(finalfile, 'w') as ins:
            for seqid in seqhash:
                ins.write(">"+ID+"\n")
                ins.write(str(seqhash[seqid])+"\n")
        ins.close()
        return

    p1 = subprocess.Popen("blastn -query "+passfile+" -db "+passfile+" -perc_identity 95 -outfmt '6 std qlen slen' ", shell=True,universal_newlines = True, stdout=subprocess.PIPE)

    gr = 0
    groups = dict()
    infodict = dict()
    done = []
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        print l
        data = l.split("\t")
        if data[0] != data[1]:
            if data[0] not in infodict:
                infodict[data[0]] = dict()
            if data[1] not in infodict[data[0]]:
                infodict[data[0]][data[1]] = dict()
            infodict[data[0]][data[1]]['qlen'] = int(data[12])
            infodict[data[0]][data[1]]['slen'] = int(data[13])
            if 'alen' not in infodict[data[0]][data[1]]:
                infodict[data[0]][data[1]]['alen'] = 0
            infodict[data[0]][data[1]]['alen'] += int(data[3])
            if 'qmin' not in infodict[data[0]][data[1]] or int(data[6]) < infodict[data[0]][data[1]]['qmin']:
                infodict[data[0]][data[1]]['qmin'] = int(data[6])
            if 'qmax' not in infodict[data[0]][data[1]] or int(data[7]) > infodict[data[0]][data[1]]['qmax']:
                infodict[data[0]][data[1]]['qmax'] = int(data[7])
            if int(data[8]) < int(data[9]):
                if 'smin' not in infodict[data[0]][data[1]] or int(data[8]) < infodict[data[0]][data[1]]['smin']:
                    infodict[data[0]][data[1]]['smin'] = int(data[8])
                if 'smax' not in infodict[data[0]][data[1]] or int(data[9]) > infodict[data[0]][data[1]]['smax']:
                    infodict[data[0]][data[1]]['smax'] = int(data[9])
            else:
                if 'smin' not in infodict[data[0]][data[1]] or int(data[9]) < infodict[data[0]][data[1]]['smin']:
                    infodict[data[0]][data[1]]['smin'] = int(data[9])
                if 'smax' not in infodict[data[0]][data[1]] or int(data[8]) > infodict[data[0]][data[1]]['smax']:
                    infodict[data[0]][data[1]]['smax'] = int(data[8])
    for i in infodict:
        for i2 in infodict[i]:
            print i+"\t"+i2+"\t"+str(infodict[i][i2]['alen'])
            # if infodict[i][i2]['qmin'] < 500 and infodict[i][i2]['qmax'] < infodict[i][i2]['qlen']-500 and (infodict[i][i2]['smin'] < 500 or infodict[i][i2]['smax'] > infodict[i][i2]['slen']-500):
            #     print "At start of query, and end/start of subject"
            #     continue
            # if infodict[i][i2]['qmax'] > infodict[i][i2]['qlen']-500 and infodict[i][i2]['qmin'] > 500 and (infodict[i][i2]['smin'] < 500 or infodict[i][i2]['smax'] > infodict[i][i2]['slen']-500):
            #     print "At end of query, and end/start of subject"
            #     continue
            if infodict[i][i2]['alen'] >= infodict[i][i2]['qlen']/2 or infodict[i][i2]['alen'] >= infodict[i][i2]['slen']/2:
                gotit = 0
                for g in groups:
                    if [(k, v) for (k, v) in groups[g].iteritems() if i in k or i2 in k]:
                        gotit += 1
                        groups[g][i+"+"+i2] = infodict[i][i2]['alen']
                if gotit == 0:
                    gr += 1
                    groups[gr] = dict()
                    groups[gr][i+"+"+i2] = infodict[i][i2]['alen']
                print "Good to align"

    out = open(midfile, 'w')
    tmpfile = dir+"/alntmp.fa"
    keep = dict()

    for g in groups:
        print "Group"+str(g)

        ids = dict()
        for i in groups[g]:
            ids[i.split("+")[0]] = 1
            ids[i.split("+")[1]] = 1

        with open(tmpfile, 'w') as ins:
            for i in ids:
                print "\t"+i
                ins.write(">"+i+"\n"+str(seqhash[i])+"\n")
        ins.close()

        p1 = subprocess.Popen("mafft --quiet --auto "+tmpfile, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

        # muscle_cline = MuscleCommandline(input=tmpfile)
        stdout, stderr = p1.communicate()
        align = AlignIO.read(StringIO(stdout), "fasta")
        print(align)

        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.dumb_consensus(ambiguous='N')
        print str(consensus).upper()

        n = 0
        for a in str(consensus).upper():
            if a == "N":
                n += 1
        pern = (float(n)/float(len(str(consensus))))*100
        print "Percent N: "+str(pern)+"%"
        if pern < 5:
            for i in ids:
                done.append(i)
            out.write(">Group"+str(g)+"\n"+str(consensus).upper()+"\n")
            keep["Group"+str(g)] = str(consensus).upper()
            continue

        if len(ids) == 2:
            continue

        print "Couldn't align all seq together, so doing it one at a time"
        alnorder = []
        aord = dict()
        conhash = dict()
        concount = 0
        for i in groups[g]:
            print "\t"+i+"\t"+str(groups[g][i])
            aord[groups[g][i]] = i
        for al in sorted(aord, reverse=True):
            print str(al)+"\t"+aord[al]
            alnorder.append(aord[al])

        newseq  = dict()

        for p in range(0,len(alnorder)):
            print "Aligning "+alnorder[p]
            i,i2 = alnorder[p].split("+")

            if i in conhash and i2 not in conhash:
                print "Adding to "+conhash[i]
                i = conhash[i]
            elif i2 in conhash and i not in conhash:
                print "Adding to "+conhash[i2]
                i2 = conhash[i2]
            elif i in conhash and i2 in conhash:
                print "Repeat!"
                continue
            with open(tmpfile, 'w') as ins:
                print "\t"+i+"\t"+i2
                ins.write(">"+i+"\n"+str(seqhash[i])+"\n")
                ins.write(">"+i2+"\n"+str(seqhash[i2])+"\n")
            ins.close()

            p1 = subprocess.Popen("mafft --quiet --auto "+tmpfile, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

            #muscle_cline = MuscleCommandline(input=tmpfile)
            stdout, stderr = p1.communicate()
            align = AlignIO.read(StringIO(stdout), "fasta")
            print(align)

            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.dumb_consensus(ambiguous='N')
            print str(consensus).upper()

            n = 0
            for a in str(consensus).upper():
                if a == "N":
                    n += 1
            print n
            pern = (float(n)/float(len(str(consensus))))*100
            print pern
            print "Percent N: "+str(pern)+"%"
            if pern < 5:
                if "Consensus" not in i and "Consensus" not in i2:
                    concount += 1
                if "Consensus" in i and "Consensus" in i2:
                    print "Aligned two consensus seq"
                    c1 = [int(s) for s in i.split() if s.isdigit()]
                    c2 = [int(s) for s in i2.split() if s.isdigit()]
                    print str(c1)+"\t"+str(c2)
                    concount += 1
                    newseq.pop(c1, None)
                    newseq.pop(c2, None)
                conhash[i] = "Consensus"+str(concount)
                conhash[i2] = "Consensus"+str(concount)
                seqhash["Consensus"+str(concount)] = consensus
                if concount not in newseq:
                    newseq[concount] = dict()
                    newseq[concount]['seq'] = consensus
                    newseq[concount]['ids'] = []
                newseq[concount]['seq'] = consensus
                newseq[concount]['ids'].append(i)
                newseq[concount]['ids'].append(i2)

        part = 0
        last = ""
        for c in newseq:
            for i in newseq[c]['ids']:
                done.append(i)
            out.write(">Group"+str(g)+"_"+str(c)+"\n"+str(newseq[c]['seq']).upper()+"\n")
            keep["Group"+str(g)+"_"+str(c)] = str(newseq[c]['seq']).upper()

    for i in seqhash:
        if i in done or "Consensus" in i:
            continue
        out.write(">"+i+"\n"+str(seqhash[i])+"\n")
        keep[i] = str(seqhash[i])

    out.close()

    if (len(list(SeqIO.parse(midfile, "fasta"))) == 1):
        print "Only one sequence left"
        subprocess.call("cat "+midfile+" | sed 's/^>.*$/>"+ID+"/' > "+finalfile, shell=True)
        return


    tmpfile = dir+"/blasttmp.fa"
    tmpfile2 = dir+"/blasttmp2.fa"
    subprocess.call("makeblastdb -dbtype nucl -in "+midfile, shell=True)

    orderdict = dict()
    p1 = subprocess.Popen('blastn -db '+midfile+' -query '+args.cDNA+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        print l
        data = l.split("\t")
        if data[0] == ID:
            if int(data[6]) not in orderdict:
                orderdict[int(data[6])] = dict()
            orderdict[int(data[6])][int(data[8])] = data[1]
    order = []
    for start in sorted(orderdict):
        for qstart in sorted(orderdict[start], reverse=True):
            print orderdict[start][qstart]+"\t"+str(start)+"\t"+str(qstart)
            if orderdict[start][qstart] not in order:
                order.append(orderdict[start][qstart])

    print order


    finalseq = []
    skip = 0

    for a in range(0,len(order)):
        if (skip == 1):
            skip = 0
            continue
        if (a == len(order) - 1):
            print "Last one"
            finalseq.append(str(keep[order[a]]))
            break

        b = a+1

        print "Current: "+order[a]+"\tNext: "+order[b]

        if len(keep[order[a]]) == 0:
            print "No sequence, moving onto next!"
            # a -= 1
            # print "Now: "+order[a]+"\tNext: "+order[b]
            continue

        with open(tmpfile, 'w') as ins:
            ins.write(">"+order[a]+"\n")
            ins.write(str(keep[order[a]])+"\n")
        ins.close()
        with open(tmpfile2, 'w') as ins:
            ins.write(">"+order[b]+"\n")
            ins.write(str(keep[order[b]])+"\n")
        ins.close()

        isgood = 0
        smin = 9999999
        smax = 0

        p1 = subprocess.Popen('blastn -subject '+tmpfile2+' -query '+tmpfile+' -outfmt 6 -evalue 1e-10',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
        for l in iter(p1.stdout.readline,''):
            l = l.rstrip()
            print l
            data = l.split("\t")
            if (int(data[7]) >= len(str(keep[order[a]])) - 500):
                isgood += 1
            if (int(data[8]) <= 500):
                isgood += 1
            if (int(data[6]) < smin):
                smin = int(data[6])-1
            if (int(data[9]) > smax):
                smax = int(data[9])
        if (isgood >= 2):
            overlap = []
            print "this is a good overlap, min: "+str(smin)+" max: "+str(smax)
            tmpseq = keep[order[a]]
            # print "First:"
            # print str(tmpseq)
            # print str(tmpseq[0:min])
            # print str(tmpseq[min:len(tmpseq)])
            overlap.append(tmpseq[smin:len(tmpseq)])
            #keep[order[a]] = tmpseq[0:min]

            # print "Second:"
            tmpseq = keep[order[b]]
            # print str(tmpseq)
            # print str(tmpseq[0:max])
            # print str(tmpseq[max:len(tmpseq)])
            overlap.append(tmpseq[0:smax])
            #keep[order[a+1]] = tmpseq[max:len(tmpseq)]

            with open(tmpfile, 'w') as ins:
                ins.write(">"+order[a]+"\n")
                ins.write(str(overlap[0])+"\n")
                ins.write(">"+order[b]+"\n")
                ins.write(str(overlap[1])+"\n")
            ins.close()

            p1 = subprocess.Popen("mafft --quiet --auto "+tmpfile, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

            #muscle_cline = MuscleCommandline(input=tmpfile)
            stdout, stderr = p1.communicate()
            align = AlignIO.read(StringIO(stdout), "fasta")
            print(align)
            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.dumb_consensus(ambiguous='N')
            print str(consensus).upper()

            n = 0
            for c in str(consensus).upper():
                if c == "N":
                    n += 1
            pern = (float(n)/float(len(str(consensus))))*100
            print "Percent N: "+str(pern)+"%"

            if pern < 15:
                tmpseq = keep[order[a]]
                keep[order[a]] = tmpseq[0:smin] + consensus
                tmpseq = keep[order[b]]
                keep[order[b]] = tmpseq[smax:len(tmpseq)]
                finalseq.append(str(keep[order[a]]).upper())
                # finalseq.append(str(consensus))
            else:
                print "Overlap not good enough"
                finalseq.append(str(keep[order[a]]))
                finalseq.append("N"*500)

        else:
            print "No overlap!"
            finalseq.append(str(keep[order[a]]))
            finalseq.append("N"*500)

    with open(finalfile, 'w') as ins:
        ins.write(">"+ID+"\n")
        ins.write("".join(finalseq)+"\n")
    ins.close()




def split_index (args):
    if os.path.exists(args.d):
        subprocess.call('rm -rf '+args.d, shell=True)
    subprocess.call("mkdir "+args.d, shell=True)
    if re.search("\.gz$", args.read1):
        subprocess.call("zcat "+args.read1+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l 4000000 -a 3 --additional-suffix=_R1.fastq - "+args.d+"/seq", shell=True)
        subprocess.call("zcat "+args.read2+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l 4000000 -a 3 --additional-suffix=_R2.fastq - "+args.d+"/seq", shell=True)
    else:
        subprocess.call("cat "+args.read1+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l 4000000 -a 3 --additional-suffix=_R1.fastq - "+args.d+"/seq", shell=True)
        subprocess.call("cat "+args.read2+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l 4000000 -a 3 --additional-suffix=_R2.fastq - "+args.d+"/seq", shell=True)

    subprocess.call('ls '+args.d+'/seq*.fastq | parallel -j '+str(args.t)+' bwa index {}', shell=True)

    with open(args.d+"/info.txt", 'w') as ins:
        ins.write(args.read1 + "\t" + args.read2)
    ins.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run an iterative "genome walking" assembly', usage='%(prog)s [options] cDNA.fa in1.fq in2.fq')
    parser.add_argument('cDNA', metavar = 'cDNA.fa', help='cDNA/mRNA sequence file')
    parser.add_argument('read1', metavar = 'in1.fq', help='FASTQ file of first paired-end reads')
    parser.add_argument('read2', metavar = 'in2.fq', help='FASTQ file of second paired-end reads')
    parser.add_argument('-t', nargs='?', metavar='INT', default=10, type=int, help='Number of processors to use (default: %(default)s)')
    parser.add_argument('-m', nargs='?', metavar='INT', default=50, type=int, help='Maximum number of iterations to make (default: %(default)s)')
    parser.add_argument('-d', nargs='?', metavar='dir/', default='Split_files/', help='Directory for fastq indexes, will be overwritten (default: %(default)s)')
    parser.add_argument('-i','--insert', nargs='?', metavar='INT', default=200, type=int, help='Average insert size (default: %(default)s)')
    parser.add_argument('-e','--endsize', nargs='?', metavar='INT', default=600, type=int, help='Number of bases from each end of the contigs to map (default: %(default)s)')
    parser.add_argument('-r','--remove', nargs='?', metavar='INT', default=200, type=int, help='After 10 iterations, remove contigs shorter than INT (default: %(default)s)')
    parser.add_argument('-f','--fastmap', nargs='?', metavar='INT', default=60, type=int, help='Minimum SMEM length permited in bwa fastmap (default: %(default)s)')
    parser.add_argument('--end_process_only', action='store_true', help='No iterative assembly will be performed, just the end process based on existing files (default: %(default)s)')

    args = parser.parse_args()

    pool = mp.Pool(processes=args.t)

    ref = args.cDNA

    subprocess.call('bwa index '+args.cDNA, shell=True)

    ids = load_ids(ref);
    print ids

    if os.path.exists(args.d) and os.path.exists(args.d+"/info.txt"):
        with open(args.d+"/info.txt", 'r') as ins:
            for l in ins:
                if not l == args.read1 + "\t" + args.read2:
                    split_index(args)
    else:
        split_index(args)

    last = dict()
    final = dict()

    if args.end_process_only:
        print "Will look at pre-existing files"
        for ID in ids:
            if os.path.exists(ID + "_files/iter" + str(args.m) + "_cap3_pass.fasta"):
                final[ID] = args.m
            elif not os.path.exists(ID+"_files/iter1_cap3_pass.fasta") or os.path.getsize(ID+"_files/iter1_cap3_pass.fasta") == 0:
                final[ID] = 0
            else:
                f = glob.glob(ID+"_files/*_cap3_pass.fasta")
                print f
                penultimate = sorted(f)[-2]
                print penultimate
                print penultimate[len(ID)+11:-16]
    else:
        for i in range(1,args.m+1):

            seqhash = dict()
            refseq = ''

            p1 = subprocess.Popen('ls '+args.d+'/seq*.fastq | parallel -k -j '+str(args.t)+' bwa fastmap -l '+ ("40" if i == 1 else str(args.fastmap)) +' {} '+ref,shell=True,universal_newlines = True, stdout=subprocess.PIPE)

            for l in iter(p1.stdout.readline,''):
                l = l.rstrip()
                data = l.split("\t")
                if re.match("SQ", l):
                    refseq = data[1]
                    refseq = re.sub("_contig.*$","",refseq)
                elif re.match("EM", l):
                    for a in range(4,len(data)):
                        id = data[a]
                        if (id == '*'):
                            continue
                        id = re.sub(":.*$","",id)
                        id = re.sub("/\d$","",id)
                        if refseq not in seqhash:
                            seqhash[refseq] = []
                        seqhash[refseq].append(id)



            new = dict()

            if i == 1:
                for ID in ids:
                    if ID not in seqhash:
                        final[ID] = 0

            idres = [pool.apply_async(assemble, args=(i,ID,seqhash[ID],args)) for ID in ids if ID not in final]
            idoutput = [p.get() for p in idres]
            #print idoutput

            for ID in idoutput:
                seqsum = 0
                maxseq = 0
                seqcount = 0
                for record in SeqIO.parse(ID + "_files/iter" + str(i) + "_cap3_pass.fasta", "fasta"):
                    seqcount += 1
                    seqsum += len(str(record.seq))
                    if len(str(record.seq)) > maxseq:
                        maxseq = len(str(record.seq))
                    if ID not in new:
                        new[ID] = dict()
                    new[ID][seqcount] = str(record.seq)
                print ID + "\t" + str(seqcount) + "\t" + str(maxseq) + "\t" + str(seqsum)
                if ID not in last:
                    if seqsum == 0:
                        print "No bases for "+ID+", exiting"
                        final[ID] = i-1
                        continue
                    last[ID] = dict()
                    last[ID]['sum'] = seqsum
                    last[ID]['max'] = maxseq
                    last[ID]['count'] = seqcount
                elif (last[ID]['sum'] >= seqsum and last[ID]['max'] >= maxseq) or (seqcount >= last[ID]['count']*3 and i > 2):
                    print "Haven't increased the total or max bp, or tripled the number of contigs for "+ID+", exiting"
                    final[ID] = i-1
                    continue
                last[ID]['sum'] = seqsum
                last[ID]['max'] = maxseq
                last[ID]['count'] = seqcount


            ref = "iter" + str(i+1) + "_ref.fasta"
            with open(ref, 'w') as ins:
                for id in new:
                    if id not in final:
                        for c in new[id]:
                            seq = new[id][c]
                            if i == 1 or len(seq) < args.endsize*2:
                                ins.write(">" + id + "_contig" + str(c) + "\n")
                                ins.write(seq + "\n")
                            elif i >= 10 and len(seq) <= args.remove:
                                pass
                            else:
                                ins.write(">" + id + "_contig" + str(c) + "_start\n")
                                ins.write(seq[:args.endsize] + "\n")
                                ins.write(">" + id + "_contig" + str(c) + "_end\n")
                                ins.write(seq[-args.endsize:] + "\n")
            ins.close()

    for ID in ids:
        if ID not in final:
            final[ID] = args.m
        if final[ID] > 0:
            pass
            #final_process(args, final[ID], ID)

    # finalres = [pool.apply_async(final_process), args=(args, final[ID], ID)) for ID in ids]
    # finaloutput = [p.get() for p in finalres]


    finalfa = "Final_sequences.fasta"
    # if os.path.exists(finalfa):
    #     subprocess.call("rm "+finalfa, shell=True)
    # for ID in ids:
        #if ID in final:
    subprocess.call("cat *_files/final_one_seq.fasta > "+finalfa, shell=True)
        #else:
        #    subprocess.call("cat "+ID+"_files/iter"+str(args.m)+"_cap3_pass.fasta | sed 's/^>/>"+ID+"_/' >> "+finalfa, shell=True)

    # fids = "Final_ids.txt"
    # f1 = "Final_R1.fastq"
    # f2 = "Final_R2.fastq"
    #
    # p1 = subprocess.Popen('ls '+args.d+'/seq*.fastq | parallel -k -j '+str(args.t)+' bwa fastmap -w 1 {} '+finalfa,shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    #
    # with open(fids, 'w') as ins:
    #     for l in iter(p1.stdout.readline,''):
    #         l = l.rstrip()
    #         data = l.split("\t")
    #     if re.match("EM", l):
    #         for a in range(4,len(data)):
    #             id = data[a]
    #             if (id == '*'):
    #                 next
    #             id = re.sub(":.*$","",id)
    #             id = re.sub("/\d$","",id)
    #             ins.write(id+'\n')
    # ins.close()
    #
    # subprocess.call("ls "+args.d+"/seq*_R1.fastq | parallel -j 2 -k 'cat {} | fqextract "+fids+"' > "+f1, shell=True)
    # subprocess.call("ls "+args.d+"/seq*_R2.fastq | parallel -j 2 -k 'cat {} | fqextract "+fids+"' > "+f2, shell=True)
    #
    # lib = "lib.txt"
    # with open(lib, 'w') as ins:
    #     ins.write("lib1 bowtie "+f1+' '+f2+' '+str(args.insert)+' 0.25 FR\n')
    # ins.close()
    #
    # subprocess.call('SSPACE_Standard_v3.0.pl -l '+lib+' -s '+finalfa+' -b Final_SSPACE', shell=True)
