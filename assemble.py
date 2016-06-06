#! /usr/bin/python
import subprocess
import sys, getopt, argparse
import re
import os.path
import multiprocessing as mp
from Bio import SeqIO
import atexit


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


    subprocess.call('cap3 ' + soapout + '.scafSeq -k 0 -p 75 -o 30', shell=True)

    cap3 = dir + "/iter" + str(i) + "_cap3.fasta"
    subprocess.call('cat ' + soapout + '.scafSeq.cap.contigs ' + soapout + '.scafSeq.cap.singlets > ' + cap3, shell=True)

    passfile = dir + "/iter" + str(i) + "_cap3_pass.fasta"
    subprocess.call('bwa mem ' + args.cDNA + ' ' + cap3 + ' | grep "'+id+'"| bam2fastx -s -M -Q -a -o ' + passfile + ' - ', shell=True)

    #count = subprocess.check_output('ls -l ' + passfile + " | awk '{print $5}' ", shell=True)

    return id

def final_process (args, i, ID):
    dir = ID + "_files"
    print dir

    passfile = dir +"/iter" + str(i) + "_cap3_pass.fasta"
    finalfile = dir +"/final_seq.fasta"
    keep = dict()

    tmpfile = dir+"/blasttmp.fa"

    subprocess.call("makeblastdb -dbtype nucl -in "+passfile, shell=True)

    for record in SeqIO.parse(passfile, "fasta"):
        with open(tmpfile, 'w') as ins:
            ins.write(">"+record.id+"\n")
            ins.write(str(record.seq)+"\n")
        ins.close()

        fail = False
        min = 9999999
        max = 0
        last = ""
        p1 = subprocess.Popen('blastn -db '+passfile+' -query '+tmpfile+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)

        for l in iter(p1.stdout.readline,''):
            l.rstrip()
            print l
            data = l.split("\t")
            if (data[0] != data[1]):
                if (last == ""):
                    last = data[1]
                elif(last != data[1]):
                    min = 9999999
                    max = 0
                if (data[6] < min):
                    min = data[6]

                if (data[7] > max):
                    max = data[7]

                print "Min: "+min+"\tMax: "+max
                if (min == 1 and max == len(str(record.seq))):
                    print "Matches a single hit from start to stop"
                    fail = True
        if (fail == False):
            print "keeping "+record.id
            keep[record.id] = str(record.seq)

    if (len(keep) == 1):
        print "Only 1 sequence is left!"
        with open(finalfile, 'w') as ins:
            for seqid in keep:
                ins.write(">"+ID+"\n")
                ins.write(keep[seqid]+"\n")
        ins.close()
        return

    orderdict = dict()
    p1 = subprocess.Popen('blastn -db '+passfile+' -query '+args.cDNA+' -outfmt 6',shell=True,universal_newlines = True, stdout=subprocess.PIPE)
    for l in iter(p1.stdout.readline,''):
        l.rstrip()
        print l
        data = l.split("\t")
        if (data[0] == ID and data[1] in keep):
            orderdict[data[6]] = data[1]
    order = []
    for start in sorted(orderdict):
        order.append(orderdict[start])

    print order


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

    for i in range(1,args.m+1):

        seqhash = dict()
        refseq = ''

        p1 = subprocess.Popen('ls '+args.d+'/seq*.fastq | parallel -k -j '+str(args.t)+' bwa fastmap -w 1 {} '+ref,shell=True,universal_newlines = True, stdout=subprocess.PIPE)

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
                        next
                    id = re.sub(":.*$","",id)
                    id = re.sub("/\d$","",id)
                    if refseq not in seqhash:
                        seqhash[refseq] = []
                    seqhash[refseq].append(id)



        new = dict()

        idres = [pool.apply_async(assemble, args=(i,ID,seqhash[ID],args)) for ID in ids if ID not in final]
        idoutput = [p.get() for p in idres]
        #print idoutput

        for ID in idoutput:
            seqsum = 0
            seqcount = 0
            for record in SeqIO.parse(ID + "_files/iter" + str(i) + "_cap3_pass.fasta", "fasta"):
                seqcount += 1
                seqsum += len(str(record.seq))
                if ID not in new:
                    new[ID] = dict()
                new[ID][seqcount] = str(record.seq)
            print ID + "\t" + seqsum
            if ID not in last:
                last[ID] = seqsum
            elif last[ID] >= seqsum:
                print "Haven't increased the total bp for "+ID+", exiting"
                final[ID] = i-1
                next
            last[ID] = seqsum


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

        final_process(args, final[ID], ID)

    # finalres = [pool.apply_async(final_process), args=(args, final[ID], ID)) for ID in ids]
    # finaloutput = [p.get() for p in finalres]


    finalfa = "Final_sequences.fasta"
    if os.path.exists(finalfa):
        subprocess.call("rm "+finalfa, shell=True)
    for ID in ids:
        #if ID in final:
            subprocess.call("cat "+ID+"_files/iter"+str(final[ID])+"_cap3_pass.fasta | sed 's/^>/>"+ID+"_/' >> "+finalfa, shell=True)
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
