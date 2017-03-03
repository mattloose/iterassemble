#! /usr/bin/env python
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
import bisect


def load_ids(file):

    idlist = []
    for record in SeqIO.parse(file, 'fasta'):
        idlist.append(record.name)
        dir = record.name + "_files"
        if not os.path.exists(dir):
            subprocess.call('mkdir '+dir, shell=True)
        with open(dir+"/transcript.fasta", 'w') as ins:
            ins.write(">"+record.name+"\n")
            ins.write(str(record.seq)+"\n")
        ins.close()

    return idlist


def assemble (i, id, arr1, args, upf1, upf2):
    dir = id + "_files"
    if not os.path.exists(dir):
        subprocess.call('mkdir '+dir, shell=True)

    f1 = dir + "/iter" + str(i) + "_R1.fastq"
    f2 = dir + "/iter" + str(i) + "_R2.fastq"
    fids = dir + "/iter" + str(i) + "_ids.txt"

    transcript = dir + "/transcript.fasta"

    with open(fids, 'w') as ins:
        ins.write("\n".join(arr1))
    ins.close()

    # if i > 1:
    #     subprocess.call("bwa index "+dir+"/iter"+str(i-1)+"_cap3_pass.fasta; bwa mem "+dir+"/iter"+str(i-1)+"_cap3_pass.fasta "+dir+"/iter"+str(i-1)+"_R1.fastq "+dir+"/iter"+str(i-1)+"_R2.fastq | samtools view -F 4 -f 8 - | awk '{print $1}' >> "+fids, shell=True)
    #
    # filelist = []
    # for fid in arr1:
    #     x = bisect.bisect(fidx[0],fid) - 1
    #     print fid + " is bigger than " + fidx[0][x]
    #     if fidx[1][x] not in filelist:
    #         filelist.append(fidx[1][x])
    # print filelist

    subprocess.call("cat "+upf1+" | fqextract "+fids+" > "+f1, shell=True)
    subprocess.call("cat "+upf2+" | fqextract "+fids+" > "+f2, shell=True)

    conf = dir + "/conf.txt"
    with open(conf, 'w') as ins:
        ins.write("max_rd_len="+str(args.length)+"\n[LIB]\navg_ins="+str(args.insert)+"\nreverse_seq=0\nasm_flags=1\nrd_len_cutoff="+str(args.length)+"\nrank=1\npair_num_cutoff=3\n")
        ins.write("q1="+f1+"\n")
        ins.write("q2="+f2+"\n")
    ins.close()

    soapout1 = dir + "/iter" + str(i) + "_63soap";
    soapout2 = dir + "/iter" + str(i) + "_fml";
    if os.path.exists(soapout1+".scafSeq"):
        subprocess.call("rm "+soapout1+"* "+soapout2+"* ", shell=True)
    soapout = dir + "/iter" + str(i) + "_soap";

    subprocess.call('timeout.sh 3600 '+args.soap+' all -s '+conf+' -K 63 -p 2 -o '+soapout1, shell=True)

    #subprocess.call('SOAPdenovo-63mer all -s '+conf+' -R -F -p 2 -o '+soapout2, shell=True)
    subprocess.call('cat ' + f1 + ' ' + f2 + " | fml-asm - | awk 'BEGIN{P=1;c=0}{if(P==1){c++;print \">FML_c\" c;}if (P==2){print}; if(P==4)P=0; P++}' > " + soapout2 + '.scafSeq', shell=True)

    subprocess.call('cat ' + soapout1 + '.scafSeq ' + soapout2 + '.scafSeq > ' +  soapout + '.scafSeq', shell=True)

    if i > 1:
        #fa = dir + '/iter' + str(i) + '.fasta'
        #subprocess.call('cat '+f1+' '+f2+' | awk \'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}\' > '+fa+ ' ; cap3 '+fa, shell=True)
        subprocess.call('bwa index '+dir+'/iter'+str(i-1)+'_cap3_pass.fasta; bwa mem '+dir+'/iter'+str(i-1)+'_cap3_pass.fasta '+f1+' '+f2+' > '+dir+'/iter'+str(i)+'.sam', shell=True)
        subprocess.call('grep ">" '+dir+'/iter'+str(i-1)+'_cap3_pass.fasta | sed \'s/^>//\' | while read -r line; do egrep "$line\s" '+dir+'/iter'+str(i)+'.sam | bam2fastx -a -A -s -Q -o '+dir+'/temp.fa - ; if [[ -s '+dir+'/temp.fa ]]; then timeout.sh 1200 cap3 '+dir+'/temp.fa; cat '+dir+'/temp.fa.cap.contigs >> ' + soapout + '.scafSeq; fi; done', shell=True)
        #subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta '+fa+'.cap.contigs '+fa+'.cap.singlets >> ' + soapout + '.scafSeq', shell=True)
        #subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta '+fa+'.cap.contigs >> ' + soapout + '.scafSeq', shell=True)
        subprocess.call('cat ' + dir + '/iter' + str(i-1) + '_cap3_pass.fasta >> ' + soapout + '.scafSeq', shell=True)


    subprocess.call('cap3 ' + soapout + '.scafSeq -k 0 -p 75 -o 30 -h 80 -f 200 -g 4', shell=True)

    cap3 = dir + "/iter" + str(i) + "_cap3.fasta"
    subprocess.call('cat ' + soapout + '.scafSeq.cap.contigs ' + soapout + ".scafSeq.cap.singlets | awk 'BEGIN{C=0} {if (/^>/){C+=1; print \">Contig_\" C;}else{print;}}' > " + cap3, shell=True)

    passfile = dir + "/iter" + str(i) + "_cap3_pass.fasta"
    if os.path.exists(passfile):
        subprocess.call("rm "+passfile,shell=True)

    keepseq = []

    subprocess.call('makeblastdb -in '+cap3+' -dbtype nucl -parse_seqids', shell=True)
    p1 = subprocess.Popen('blastn -db '+cap3+' -query '+transcript+' -outfmt 6 -culling_limit '+str(args.culling),shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        data = l.split("\t")
        if data[0] == id and data[1] not in keepseq:
            keepseq.append(data[1])
            subprocess.call('blastdbcmd -db '+cap3+' -entry '+data[1]+' -outfmt "%s" | awk \'BEGIN{print ">'+data[1]+'"}{print}\' >> '+passfile, shell=True)

    addseq = []

    p1 = subprocess.Popen('bwa index '+cap3+' ; bwa mem '+cap3+' '+f1+' '+f2+' | samtools view -F 2316 - | awk \'{if ($7 != "="){ print;}}\' | gawk \'{a[$3][$7]++} END {for (b in a){ for (c in a[b]){ print b "\t" c "\t" a[b][c];}}}\' ',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        print l
        data = l.split("\t")
        if data[0] in keepseq and data[2] >= args.nreads:
            if data[1] not in addseq and data[1] not in keepseq:
                addseq.append(data[1])
                subprocess.call('blastdbcmd -db '+cap3+' -entry '+data[1]+' -outfmt "%s" | awk \'BEGIN{print ">'+data[1]+'"}{print}\' >> '+passfile, shell=True)
        if data[1] in keepseq and data[2] >= args.nreads:
            if data[0] not in addseq and data[0] not in keepseq:
                addseq.append(data[0])
                subprocess.call('blastdbcmd -db '+cap3+' -entry '+data[0]+' -outfmt "%s" | awk \'BEGIN{print ">'+data[0]+'"}{print}\' >> '+passfile, shell=True)


    # passfile = dir + "/iter" + str(i) + "_cap3_pass.fasta"
    # subprocess.call('bwa mem ' + args.cDNA + ' ' + cap3 + ' | grep "'+id+'"| bam2fastx -s -M -Q -a -o ' + passfile + ' - ', shell=True)

    #count = subprocess.check_output('ls -l ' + passfile + " | awk '{print $5}' ", shell=True)

    return id

def final_process (args, i, ID):
    dir = ID + "_files"
    print dir

    finallog = dir + "/final.log"
    finallogout = open(finallog, "w")

    transcript = dir + "/transcript.fasta"

    infile = dir +"/iter" + str(i) + "_cap3_pass.fasta"
#    scaffile = dir +"/iter" + str(i) + "_cap3_pass.scaffolds.fasta"
    allr1 = dir + "/allR1.fastq"
    allr2 = dir + "/allR2.fastq"
    subprocess.call("cat "+dir+"/*_R1.fastq > "+allr1, shell=True)
    subprocess.call("cat "+dir+"/*_R2.fastq > "+allr2, shell=True)
    subprocess.call("rmPairedDuplicates.py "+allr1+" "+allr2, shell=True)

    passfile1 = dir +"/iter" + str(i) + "_cap3_pass.fasta.renamed"

    sc = 0
    startseqhash = dict()
    with open(passfile1, 'w') as ins:
        for record in SeqIO.parse(infile, "fasta"):
            if i >= 10 and len(str(record.seq)) <= args.remove:
                continue
            sc += 1
            ins.write(">Contig"+str(sc)+"\n"+str(record.seq)+"\n")
            startseqhash["Contig"+str(sc)] = record.seq

    ins.close()

    ## Add section that checks for chimeric contigs

    passfile = dir +"/iter" + str(i) + "_cap3_pass.fasta.renamed.nonchimeric"
    nonchiout = open(passfile, 'w')

    chires = dict()
    problems = []
    subprocess.call("makeblastdb -dbtype nucl -in "+passfile1, shell=True)
    p1 = subprocess.Popen('blastn -db '+passfile1+' -query '+transcript+' -outfmt 6 ',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        finallogout.write(l+"\n")
        data = l.split("\t")
        if data[1] not in chires:
            chires[data[1]] = dict()
            chires[data[1]]['st'] = []
            chires[data[1]]['en'] = []
        chires[data[1]]['st'].append(data[6])
        chires[data[1]]['en'].append(data[7])
    for c in chires:
        if len(chires[c]['st']) > 1:
            finallogout.write(c+"\n")
            for a in range(len(chires[c]['st'])):
                finallogout.write(str(a)+"\t"+chires[c]['st'][a] + "\t" + chires[c]['en'][a]+"\n")
                good = 0
                for b in range(len(chires[c]['st'])):
                    if b == a:
                        continue
                    if int(chires[c]['st'][a]) <= int(chires[c]['en'][b]) + 20 and int(chires[c]['st'][a]) >= int(chires[c]['en'][b]) - 20:
                        good += 1
                    if int(chires[c]['en'][a]) <= int(chires[c]['st'][b]) + 20 and int(chires[c]['en'][a]) >= int(chires[c]['st'][b]) - 20:
                        good += 1
                if good == 0:
                    finallogout.write("PROBLEMATIC\n")
                    if c not in problems:
                        problems.append(c)
    for c in problems:
        contigtempfile = dir + "/contig.temp.fasta"
        contigtempsam = dir +  "/contig.temp.sam"
        with open(contigtempfile, 'w') as ins:
            ins.write(">"+c+"\n")
            ins.write(str(startseqhash[c])+"\n")
        ins.close()
        subprocess.call('bwa index '+contigtempfile+'; bwa mem '+contigtempfile+' '+allr1+"_rmDup "+allr2+"_rmDup > "+contigtempsam, shell=True)
        contigtempperfect = dir + "/contig.temp.perfect.sam"
        subprocess.call('samtools view -H '+contigtempsam+' > '+contigtempperfect, shell=True)
        subprocess.call('samtools view -F 12 '+contigtempsam+' | grep "[[:space:]]100M[[:space:]]" | grep "MD:Z:100" >> '+contigtempperfect, shell=True)
        contigtempbam = dir + "/contig.temp.perfect.bam"
        subprocess.call("samtools sort -n -o "+contigtempbam+" "+contigtempperfect, shell=True)
        contigtempdepth = dir + "/contig.temp.perfect.depth.txt"
        with open(dir + "/chrom.sizes", 'w') as ins:
            ins.write(c+"\t"+str(len(str(startseqhash[c])))+"\n")
        ins.close()
        subprocess.call('bedtools bamtobed -i '+contigtempbam+' -bedpe | cut -f 1,2,6 | sort -k1,1 | bedtools genomecov -i - -g '+dir+'/chrom.sizes -d > '+contigtempdepth, shell=True)
        regions = dict()
        with open(contigtempdepth, 'r') as ins:
            p = -5
            b = 0
            for l in ins:
                l = l.rstrip()
                data = l.split("\t")
                if int(data[2]) == 0:
                    finallogout.write(l+'\n')
                    if int(data[1]) == p + 1:
                        p += 1
                    else:
                        if b in regions:
                            regions[b]['end'] = p
                        b += 1
                        if b not in regions:
                            regions[b] = dict()
                            regions[b]['start'] = int(data[1])
                        p = int(data[1])
            if b in regions:
                regions[b]['end'] = p
        finallogout.write(str(regions)+"\n")
        partcount = 0
        lastb = 0
        for b in regions:
            finallogout.write(str(b)+"\n")
            if regions[b]['start'] == 1 or regions[b]['end'] == len(str(startseqhash[c])):
                finallogout.write("At start/end\n")
            else:
                partcount += 1
                prev = 1
                if b > 1 and regions[b-1]['start'] != 1:
                    prev = regions[b-1]['end']
                newchunk = str(startseqhash[c])[prev:regions[b]['start']]
                finallogout.write("New chunk from "+ str(prev)+ " - "+str(regions[b]['start'])+"\n")
                nonchiout.write(">"+c+"_"+str(partcount)+"\n"+newchunk+"\n")
                lastb = b
        if partcount > 0:
            partcount += 1
            newchunk = str(startseqhash[c])[regions[lastb]['end']:len(startseqhash[c])]
            finallogout.write("New chunk from "+ str(regions[lastb]['end'])+ " - "+str(len(startseqhash[c]))+"\n")
            nonchiout.write(">"+c+"_"+str(partcount)+"\n"+newchunk+"\n")
        else:
            nonchiout.write(">"+c+"\n"+str(startseqhash[c])+"\n")

    for c in startseqhash:
        if c in problems:
            continue
        nonchiout.write(">"+c+"\n"+str(startseqhash[c])+"\n")

    nonchiout.close()

    # bam = dir + "/iter"+str(i)+"_cap3_pass.bam"
    # subprocess.call("bwa index "+infile+"; bwa mem "+infile+" "+allr1+"_rmDup "+allr2+"_rmDup | samtools view -b -F 2048 - > "+bam, shell=True)
    # subprocess.call("sga-bam2de.pl -n 5 --prefix "+dir+"/sgascaf "+bam, shell=True)
    # #subprocess.call("sga-astat.py "+bam+" > "+dir+"/sgascaf.astat", shell=True)
    # #subprocess.call("sga scaffold -m 200 --pe "+dir+"/sgascaf.de -a "+dir+"/sgascaf.astat -o "+dir+"/sgascaf.scaf "+infile, shell=True)
    # subprocess.call("sga scaffold -m 200 --pe "+dir+"/sgascaf.de -o "+dir+"/sgascaf.scaf "+infile, shell=True)
    # subprocess.call("sga scaffold2fasta -o "+scaffile+" -f "+infile+" "+dir+"/sgascaf.scaf", shell=True)

    # if not os.path.exists(scaffile):
    #     subprocess.call("cp "+infile+" "+scaffile, shell=True)


    midfile = dir + "/final_all_seq.fasta"

    finalfile = dir +"/final_one_seq.fasta"

    seqhash = dict()
    for record in SeqIO.parse(passfile, "fasta"):
        seqhash[record.name] = record.seq

    subprocess.call("makeblastdb -dbtype nucl -in "+passfile, shell=True)

    revcom = dict()
    p1 = subprocess.Popen('blastn -db '+passfile+' -query '+transcript+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        finallogout.write(l+"\n")
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
        finallogout.write("Only one sequence\n")
        with open(finalfile, 'w') as ins:
            for seqid in seqhash:
                ins.write(">"+ID+"\n")
                ins.write(str(seqhash[seqid])+"\n")
        ins.close()
        return ID

    p1 = subprocess.Popen("blastn -query "+passfile+" -db "+passfile+" -perc_identity 95 -outfmt '6 std qlen slen' ", shell=True,universal_newlines = True, stdout=subprocess.PIPE)

    gr = 0
    groups = dict()
    infodict = dict()
    done = []
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        finallogout.write(l+"\n")
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
            finallogout.write(i+"\t"+i2+"\t"+str(infodict[i][i2]['alen'])+"\n")
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
                finallogout.write("Good to align\n")

    out = open(midfile, 'w')
    tmpfile = dir+"/alntmp.fa"
    keep = dict()

    for g in groups:
        finallogout.write("Group"+str(g)+"\n")

        ids = dict()
        for i in groups[g]:
            ids[i.split("+")[0]] = 1
            ids[i.split("+")[1]] = 1

        with open(tmpfile, 'w') as ins:
            for i in ids:
                finallogout.write("\t"+i+"\n")
                ins.write(">"+i+"\n"+str(seqhash[i])+"\n")
        ins.close()

        p1 = subprocess.Popen("mafft --quiet --auto "+tmpfile, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

        # muscle_cline = MuscleCommandline(input=tmpfile)
        stdout, stderr = p1.communicate()
        align = AlignIO.read(StringIO(stdout), "fasta")
        finallogout.write(str(align)+"\n")

        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.dumb_consensus(ambiguous='N')
        finallogout.write(str(consensus).upper()+"\n")

        n = 0
        for a in str(consensus).upper():
            if a == "N":
                n += 1
        pern = (float(n)/float(len(str(consensus))))*100
        finallogout.write("Percent N: "+str(pern)+"%\n")
        if pern < 5:
            for i in ids:
                done.append(i)
            out.write(">Group"+str(g)+"\n"+str(consensus).upper()+"\n")
            keep["Group"+str(g)] = str(consensus).upper()
            continue

        if len(ids) == 2:
            continue

        finallogout.write("Couldn't align all seq together, so doing it one at a time\n")
        alnorder = []
        aord = dict()
        conhash = dict()
        concount = 0
        for i in groups[g]:
            finallogout.write("\t"+i+"\t"+str(groups[g][i])+"\n")
            aord[groups[g][i]] = i
        for al in sorted(aord, reverse=True):
            finallogout.write(str(al)+"\t"+aord[al]+"\n")
            alnorder.append(aord[al])

        newseq  = dict()

        for p in range(0,len(alnorder)):
            finallogout.write("Aligning "+alnorder[p]+"\n")
            i,i2 = alnorder[p].split("+")

            if i in conhash and i2 not in conhash:
                finallogout.write("Adding to "+conhash[i]+"\n")
                i = conhash[i]
            elif i2 in conhash and i not in conhash:
                finallogout.write("Adding to "+conhash[i2]+"\n")
                i2 = conhash[i2]
            elif i in conhash and i2 in conhash:
                finallogout.write("Repeat!\n")
                continue
            with open(tmpfile, 'w') as ins:
                finallogout.write("\t"+i+"\t"+i2+"\n")
                ins.write(">"+i+"\n"+str(seqhash[i])+"\n")
                ins.write(">"+i2+"\n"+str(seqhash[i2])+"\n")
            ins.close()

            p1 = subprocess.Popen("mafft --quiet --auto "+tmpfile, shell=True, universal_newlines = True, stdout=subprocess.PIPE)

            #muscle_cline = MuscleCommandline(input=tmpfile)
            stdout, stderr = p1.communicate()
            align = AlignIO.read(StringIO(stdout), "fasta")
            finallogout.write(str(align)+"\n")

            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.dumb_consensus(ambiguous='N')
            finallogout.write(str(consensus).upper()+'\n')

            n = 0
            for a in str(consensus).upper():
                if a == "N":
                    n += 1
            pern = (float(n)/float(len(str(consensus))))*100
            finallogout.write("Percent N: "+str(pern)+"%\n")
            if pern < 5:
                if "Consensus" not in i and "Consensus" not in i2:
                    concount += 1
                if "Consensus" in i and "Consensus" in i2:
                    finallogout.write("Aligned two consensus seq\n")
                    c1 = [int(s) for s in i.split() if s.isdigit()]
                    c2 = [int(s) for s in i2.split() if s.isdigit()]
                    finallogout.write(str(c1)+"\t"+str(c2)+"\n")
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
        finallogout.write("Only one sequence left\n")
        subprocess.call("cat "+midfile+" | sed 's/^>.*$/>"+ID+"/' > "+finalfile, shell=True)
        return ID


    tmpfile = dir+"/blasttmp.fa"
    tmpfile2 = dir+"/blasttmp2.fa"
    subprocess.call("makeblastdb -dbtype nucl -in "+midfile, shell=True)


    ### Modify this section to spot contigs which map multiple 'exons' out of sync. Select the longest exon?

    orderdict = dict()
    p1 = subprocess.Popen('blastn -db '+midfile+' -query '+transcript+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l = l.rstrip()
        finallogout.write(l+"\n")
        data = l.split("\t")
        if data[0] == ID:
            if int(data[6]) not in orderdict:
                orderdict[int(data[6])] = dict()
            orderdict[int(data[6])][int(data[8])] = data[1]
    order = []
    for start in sorted(orderdict):
        for qstart in sorted(orderdict[start], reverse=True):
            finallogout.write(orderdict[start][qstart]+"\t"+str(start)+"\t"+str(qstart)+"\n")
            if orderdict[start][qstart] not in order:
                order.append(orderdict[start][qstart])

    finallogout.write(str(order)+"\n")


    finalseq = []
    skip = 0

    for a in range(0,len(order)):
        if (skip == 1):
            skip = 0
            continue
        if (a == len(order) - 1):
            finallogout.write("Last one\n")
            finalseq.append(str(keep[order[a]]))
            break

        b = a+1

        finallogout.write("Current: "+order[a]+"\tNext: "+order[b]+"\n")

        if len(keep[order[a]]) == 0:
            finallogout.write("No sequence, moving onto next!\n")
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
            finallogout.write(l+"\n")
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
            finallogout.write("this is a good overlap, min: "+str(smin)+" max: "+str(smax)+"\n")
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
            finallogout.write(str(align)+"\n")
            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.dumb_consensus(ambiguous='N')
            finallogout.write(str(consensus).upper()+"\n")

            n = 0
            for c in str(consensus).upper():
                if c == "N":
                    n += 1
            pern = (float(n)/float(len(str(consensus))))*100
            finallogout.write("Percent N: "+str(pern)+"%\n")

            if pern < 15:
                tmpseq = keep[order[a]]
                keep[order[a]] = tmpseq[0:smin] + consensus
                tmpseq = keep[order[b]]
                keep[order[b]] = tmpseq[smax:len(tmpseq)]
                finalseq.append(str(keep[order[a]]).upper())
                # finalseq.append(str(consensus))
            else:
                finallogout.write("Overlap not good enough\n")
                finalseq.append(str(keep[order[a]]))
                finalseq.append("N"*500)

        else:
            finallogout.write("No overlap!\n")
            finalseq.append(str(keep[order[a]]))
            finalseq.append("N"*500)

    with open(finalfile, 'w') as ins:
        ins.write(">"+ID+"\n")
        ins.write("".join(finalseq)+"\n")
    ins.close()

    finallogout.close()

    return ID




def split_index (args):
    if os.path.exists(args.d):
        subprocess.call('rm -rf '+args.d, shell=True)
    subprocess.call("mkdir "+args.d, shell=True)
    if re.search("\.gz$", args.read1):
        subprocess.call("gzip -dc "+args.read1+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l "+str(args.split)+" -a 3 - "+args.d+"/seq", shell=True)
        subprocess.call("ls "+args.d+"/seq??? | awk '{system(\"mv \" $0 \" \" $0 \"_R1.fastq\")}'", shell=True)
        subprocess.call("gzip -dc "+args.read2+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l "+str(args.split)+" -a 3 - "+args.d+"/seq", shell=True)
        subprocess.call("ls "+args.d+"/seq??? | awk '{system(\"mv \" $0 \" \" $0 \"_R2.fastq\")}'", shell=True)
    else:
        subprocess.call("cat "+args.read1+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l "+str(args.split)+" -a 3 - "+args.d+"/seq", shell=True)
        subprocess.call("ls "+args.d+"/seq??? | awk '{system(\"mv \" $0 \" \" $0 \"_R1.fastq\")}'", shell=True)
        subprocess.call("cat "+args.read2+" | awk 'BEGIN{P=1}{if(P==1){gsub(/\s+.*$/,\"\"); gsub(/\/[1,2]$/, \"\")}; print; if(P==4)P=0; P++}' - | split -l "+str(args.split)+" -a 3 - "+args.d+"/seq", shell=True)
        subprocess.call("ls "+args.d+"/seq??? | awk '{system(\"mv \" $0 \" \" $0 \"_R2.fastq\")}'", shell=True)

    subprocess.call('ls '+args.d+'/seq*.fastq | parallel -j '+str(args.t)+' bwa index {}', shell=True)
    subprocess.call('ls '+args.d+'/seq*.fastq | parallel -j '+str(args.t)+' gzip {}', shell=True)
    #subprocess.call("ls "+args.d+"/seq*_R1.fastq | parallel -k --tag head -n 1 {} | awk '{print $2 \"\t\" $1}' | sed 's/^@//' | sed 's/_R1.fastq$//' > "+args.d+"/fq_to_file.txt", shell=True)
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
    parser.add_argument('-l','--length', nargs='?', metavar='INT', default=100, type=int, help='Maximum read length (default: %(default)s)')
    parser.add_argument('-e','--endsize', nargs='?', metavar='INT', default=600, type=int, help='Number of bases from each end of the contigs to map (default: %(default)s)')
    parser.add_argument('-r','--remove', nargs='?', metavar='INT', default=200, type=int, help='After 10 iterations, remove contigs shorter than INT (default: %(default)s)')
    parser.add_argument('-F','--fastmap1', nargs='?', metavar='INT', default=40, type=int, help='Minimum SMEM length permited in bwa fastmap for first iteration (default: %(default)s)')
    parser.add_argument('-f','--fastmap', nargs='?', metavar='INT', default=60, type=int, help='Minimum SMEM length permited in bwa fastmap for later iterations (default: %(default)s)')
    parser.add_argument('-c','--culling', nargs='?', metavar='INT', default=2, type=int, help="After each iteration blast step uses this culling_limit (default: %(default)s)")
    parser.add_argument('-n','--nreads', nargs='?', metavar='INT',default=5,type=int, help = "Contigs with n reads mapped across are kept after each iteration (default: %(default)s)")
    parser.add_argument('-M','--maxcontigs', nargs='?', metavar= 'INT', default=500, type=int, help="Maximum number of contigs per gene after each iteration (default: %(default)s)")
    parser.add_argument('-s','--split', nargs = '?', metavar='INT',default=4000000, type=int, help="FASTQ files will be split on this line number, must be divisible by 4 (default: %(default)s)")
    parser.add_argument('--continue_from', nargs='?', metavar='INT', default=0, type=int, help="Continue iterating from this number (default: %(default)s)")
    parser.add_argument('--soap', nargs = '?', default = "SOAPdenovo-63mer", help = "SOAP denovo program to use, must be in path (default: %(default)s)")
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

    # fidx1 = []
    # fidx2 = []
    # with open(args.d+"/fq_to_file.txt", 'r') as ins:
    #     for l in ins:
    #         s = l.rstrip().split("\t")
    #         fidx1.append(s[0])
    #         fidx2.append(s[1])
    # fidx = [fidx1,fidx2]
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
                if len(f) > 1:
                    penultimate = f[-2]
                else:
                    penultimate = f[0]
                final[ID] = int(penultimate[len(ID)+11:-16])
    else:

        logfile = "iterassemble.log"
        logout = open(logfile, 'w')
        logout.write("Iter\tID\tContig no.\tMax contig length\tTotal length\tStatus\n")

        if args.continue_from > 0:
            ref = "iter" + str(args.continue_from+1) + "_ref.fasta"

        for i in range(args.continue_from + 1, args.m + 1):

            seqhash = dict()
            refseq = ''
            idsfile = "iter"+str(i)+"_ids.txt"
            iout = open(idsfile, 'w')

            p1 = subprocess.Popen('ls '+args.d+'/seq*.fastq.gz | sed \'s/\.gz$//\' | parallel -k -j '+str(args.t)+' bwa fastmap -l '+ (str(args.fastmap1) if i == 1 else str(args.fastmap)) +' {} '+ref,shell=True,universal_newlines = True, stdout=subprocess.PIPE)

            for l in iter(p1.stdout.readline,''):
                l = l.rstrip()
                data = l.split("\t")
                if re.match("/*/*SQ", l):
                    refseq = data[1]
                    refseq = re.sub("_contig.*$","",refseq)
                elif re.match("EM", l):
                    for a in range(4,len(data)):
                        id = data[a]
                        if (id == '*'):
                            continue
                        id = re.sub(":[-+]\d+$","",id)
                        id = re.sub("/\d$","",id)
                        iout.write(id+"\n")
                        if refseq not in seqhash:
                            seqhash[refseq] = []
                        seqhash[refseq].append(id)

            iout.close()

            f1 = "iter"+str(i)+"_R1.fastq"
            f2 = "iter"+str(i)+"_R2.fastq"
            ex1 = subprocess.Popen("ls "+args.d+"/*_R1.fastq.gz | parallel -k -j 1 'gzip -dc {} | fqextract "+idsfile+"' > "+f1, shell=True)
            subprocess.call("ls "+args.d+"/*_R2.fastq.gz | parallel -k -j 1 'gzip -dc {} | fqextract "+idsfile+"' > "+f2, shell=True)
            ex1.wait()

            new = dict()

            for ID in ids:
                if ID not in final and ID not in seqhash:
                    final[ID] = i-1

            idres = [pool.apply_async(assemble, args=(i,ID,seqhash[ID],args,f1,f2)) for ID in ids if ID not in final]
            idoutput = [p.get() for p in idres]
            #print idoutput

            for ID in idoutput:
                seqsum = 0
                maxseq = 0
                seqcount = 0
                if os.path.exists(ID + "_files/iter" + str(i) + "_cap3_pass.fasta"):
                    for record in SeqIO.parse(ID + "_files/iter" + str(i) + "_cap3_pass.fasta", "fasta"):
                        seqcount += 1
                        seqsum += len(str(record.seq))
                        if len(str(record.seq)) > maxseq:
                            maxseq = len(str(record.seq))
                        if ID not in new:
                            new[ID] = dict()
                        new[ID][seqcount] = str(record.seq)
                print ID + "\t" + str(seqcount) + "\t" + str(maxseq) + "\t" + str(seqsum)
                logout.write(str(i)+"\t"+ID + "\t" + str(seqcount) + "\t" + str(maxseq) + "\t" + str(seqsum))
                if ID not in last:
                    if seqsum == 0:
                        print "No bases for "+ID+", exiting"
                        final[ID] = i-1
                        continue
                    last[ID] = dict()
                    last[ID]['sum'] = seqsum
                    last[ID]['max'] = maxseq
                    last[ID]['count'] = seqcount
                elif (last[ID]['sum'] >= seqsum and last[ID]['max'] >= maxseq) or (seqcount >= last[ID]['count']*3 and i > 2) or (seqcount >= args.maxcontigs):
                    print "Haven't increased the total or max bp, or tripled the number of contigs for "+ID+", or found too many contigs, exiting"
                    logout.write("\tENDED\n")
                    final[ID] = i-1
                    continue
                logout.write("\n")
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

            if os.path.getsize(ref) == 0:
                break


        logout.close()

    for ID in ids:
        if ID not in final:
            final[ID] = args.m
        # if final[ID] > 0:
        #     final_process(args, final[ID], ID)

    finalres = [pool.apply_async(final_process, args=(args, final[ID], ID)) for ID in ids if final[ID] > 0]
    finaloutput = [p.get() for p in finalres]


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
