#!/usr/bin/env python

import os, sys, argparse
import itertools

def ParseArg():
    p=argparse.ArgumentParser( description = 'Remove paired fastq reads where at least one sequence contains a series of Ns', epilog = 'Library dependency: itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq file')
    p.add_argument('ntimes',type=int,metavar='Nx', help='Number of Ns to look for')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


if __name__ == '__main__':

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
            if N*args.ntimes in str(nextr1.rstrip()) or N*args.ntimes in str(nextr2.rstrip()):
                read1.next()
                read1.next()
                read2.next()
                read2.next()
            else:
                outfile1.write(keepr1+nextr1+read1.next()+read1.next())
                outfile2.write(keepr2+nextr2+read2.next()+read2.next())
