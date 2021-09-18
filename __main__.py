# python2.7
# -*- coding: utf-8 -*-
'''
Log:
released on Sept, 2021
'''

import os, sys
import argparse
import time

script_dir = os.path.dirname( __file__ )
sys.path.append(os.path.join(script_dir, 'src'))

try:
    from MainForCountCoverage import *
    import toolpath
    import PlotMatrix
    import PrepareAnnotation
except:
    print "Error: transcript_body_coverage scripts are not found"


def Arg():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')

    # Count the transcript body coverage
    parser_a = subparsers.add_parser('count', help='Count the transcript body coverage')
    parser_a.add_argument('--bam', help='RNAseq or 5end seq bam file', dest='bam')
    parser_a.add_argument('--number', help='number of bins', type=int, dest='number', default=101)
    parser_a.add_argument('--size_cutoff', help='transcript size cutoff', type=int, dest='size_cutoff', default=100)
    parser_a.add_argument('--type', help='type of input file: RNAseq of TSS', dest='type')
    parser_a.add_argument('--output', help='output matrix for plotting', dest='output_file', default='transcript_coverage.matrix')
    
    parser_a.add_argument('--samtools', help='path for samtools, e.g. ~/miniconda2/envs/my-python2-env/bin/', dest='samtoolspath', default='') 
    parser_a.add_argument('--bedtools', help='path for bedtools, e.g. ~/miniconda2/envs/my-python2-env/bin/', dest='bedtoolspath', default='') 
    
    parser_a.add_argument('--gtf', help='exon annotation generated using annotation function', dest='exon_gtf', default=None)
    parser_a.add_argument('--bin_gtf', help='encode annotation file with split bins for exons', dest='bin_gtf', default=None) # e.g. gencode.v24.exonForcoverage.100bin.gtf

    # conditional arguments for RNA-seq bam
    parser_a.add_argument('--FPKM_file', help='FPKM file from RSeQC', dest='FPKM_file', default=None)
    parser_a.add_argument('--FPKM_cutoff', help='FPKM cutoff for RNA-seq coverage', dest='FPKM_cutoff', type=float, default=10)
    
    # Plot the matrix
    parser_b = subparsers.add_parser('plot', help='Plot the transcript body coverage matrix')
    parser_b.add_argument('--input', nargs = '+', help = 'coverage matrix', dest='report_list')
    parser_b.add_argument('--name', nargs = '+', help = 'name label for each matrix file', dest='name_list', required=False)
    parser_b.add_argument('--png', help='output plot name', dest='png', default='transcript.coverage.png')
    parser_b.add_argument('--count_cutoff', help='count cutoff for TSS coverage', dest='Endseq_cutoff', type=int, default=-1)

    # Prepare the exon anntation file used for counting
    parser_c = subparsers.add_parser('annotation', help='prepare the exon anntation file used for counting')
    parser_c.add_argument('--input', help = 'gencode exon annotation file', dest='gtf')
    parser_c.add_argument('--output', help = 'exon annotation used by count function', dest='output_file')

    args = parser.parse_args()

    return args


#-------mainbody
if __name__ == '__main__':
    
    args = Arg()

    # Find the full path of all the required files
    localtime = time.asctime(time.localtime()) # return e.g. Mon Feb  4 15:33:23 2019
    print "Start at:", localtime

   
    if args.mode =='count':
        
        # Test bedtools and samtools, change the path based on where tools are installed
        toolpath.init(args.samtoolspath, args.bedtoolspath)
        toolresult = toolpath.Tools() # bedtools return, samtools return
        print toolresult[0]
        print toolresult[1]

        type = args.type
        bam = os.path.abspath(args.bam)
        if not os.path.exists(''.join([bam, '.bai'])):
            print "Error: The index file bam.bai is not found."
            quit()
        output_file = os.path.abspath(args.output_file)
        size_cutoff = args.size_cutoff
        number = args.number
        
        # neee at least one file provided
        if args.exon_gtf: 
            exon_gtf = os.path.abspath(args.exon_gtf)
        else:
            exon_gtf = None
        if args.bin_gtf:
            bin_gtf = os.path.abspath(args.bin_gtf)
        else:
            bin_gtf = None
        if exon_gtf or bin_gtf:
            pass
        else:
            print "Error: An Exon gtf file or A bin file needs to be provided."
            quit()        

        tempdir = CreateTemp()
        os.mkdir(tempdir)

        if type == 'RNAseq':
            if args.FPKM_file:
                FPKMfile = os.path.abspath(args.FPKM_file)
            else:
                FPKMfile = None
            FPKM_cutoff = args.FPKM_cutoff
            os.chdir(tempdir)
            
            try:
                CountCoverageRNAseq(bam, FPKMfile, FPKM_cutoff, exon_gtf, number, size_cutoff, output_file, bin_gtf)
            except Exception as e:
                print e
                os.chdir('..')
                check_call('rm -r {}'.format(tempdir), shell=True)
                quit()
        else:
            os.chdir(tempdir)
            try:
                CountCoverageEndseq(bam, exon_gtf, number, size_cutoff, output_file, bin_gtf)
            except Exception as e:
                print e
                os.chdir('..')
                check_call('rm -r {}'.format(tempdir), shell=True)
                quit()
        
        os.chdir('..')
        check_call('rm -r {}'.format(tempdir), shell=True)
    
    elif args.mode =='plot':
        # if name not provide, use 'seq1, seq2...'
        if not args.name_list:
            ls = ['seq']*3
            name_list = []
            for idx, item in enumerate(ls, 1):
                name_list.append('{}{}'.format(item, idx))
        else:
            name_ls = args.name_list
        PlotMatrix.Plot(args.report_list, name_list, args.png, args.Endseq_cutoff)
    elif args.mode =='annotation':
        PrepareAnnotation.PrepareGencode(args.gtf, args.output_file)
    else:
        print 'Error: Choose either annotation, count or plot as subparser.'

    localtime = time.asctime(time.localtime())
    print "End at :", localtime  
