#!/usr/bin/env python

import os, sys, argparse
import itertools

def ParseArg():
    p=argparse.ArgumentParser( description = 'Remove duplicated fastq reads which have same sequences for both forward and reverse reads. Choose the one appears first.', epilog = 'Library dependency: itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq file')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


if __name__ == '__main__':

    Unique_seqs=set()
    args=ParseArg()

    outfile1 = open(args.input1+"_rmDup","w")
    outfile2 = open(args.input2+"_rmDup","w")

    linecount = 0

    with open(args.input1, 'rt') as read1, open(args.input2, 'rt') as read2:
        for r1, r2 in itertools.izip(read1,read2):
            keepr1 = r1
            keepr2 = r2
            nextr1 = read1.next()
            nextr2 = read2.next()
            combi = str(nextr1.rstrip()+nextr2.rstrip())
            if combi not in Unique_seqs:
                Unique_seqs.add(combi)
                outfile1.write(keepr1+nextr1+read1.next()+read1.next())
                outfile2.write(keepr2+nextr2+read2.next()+read2.next())
            else:
                read1.next()
                read1.next()
                read2.next()
                read2.next()
