# python2.7
# -*- coding: utf-8 -*-
'''
main for counting the transcript body covarge
'''

from CountCoverage import *
from subprocess import check_call

def CountCoverageRNAseq(bam, FPKMfile, FPKM_cutoff, exon_gtf, number, size_cutoff, output_file, bin_gtf):
    '''
    For RNAseq data
    
    bam: RNA-seq bam file containing primary reads
    FPKMfile: 
    A file generated by RSeQC FPKM_count.py function, e.g.
    #chrom  st      end     accession       mRNA_size       gene_strand     Frag_count      FPM     FPKM
    chr1    8352396 8817465 ENST00000337907.7       8026.0  -       489.0   121.103992915   15.0889599944

    FPKM_cutoff: Only transcripts with FPKM ≥ cutoff are used for RNA-seq coverage calculation
    exon_gtf: An annotation gtf file containing exon information for each transcript.
    number: number of bins used for each transcript
    size_cutoff: 
    Only transcripts with length ≥ max(size_cutoff, number-1) will be used for RNA-seq coverage calculation.
    output_file: a matrix file can be used for plotting
    transcriptID, count in Bin1 ..., Count in Bin number-1
    '''
    print "Calculate the transcript body coverage for RNAseq data {}.".format(bam)
    print "RNAseq bam file: {}".format(bam)
    print "The number of bins used: {}".format(number-1)
    print "Transcript size cutoff: {}".format(size_cutoff)
    print "The output matrix: {}".format(output_file)

    
    if FPKMfile:
        if not exon_gtf:
            print "Error: An Exon gtf file needs to be provided to create bins."
            quit()
        test_gtf = BinGtf(exon_gtf, number)
        test_gtf.ExtractGeneRNAseq(FPKMfile, FPKM_cutoff)
        test_gtf.generateGtfBin('Transcript.coverage.bin', size_cutoff)
        bin_gtf_file = 'Transcript.coverage.bin'
        print "RNAseq FPKM file: {}".format(FPKMfile)
        print "Transcripts with FPKM>={} are used for analysis.".format(FPKM_cutoff)
    else:
        print "Not find FPKM file. So Use all the transripts in the exon gtf for analysis."
        if bin_gtf:
            bin_gtf_file = bin_gtf
            print "Use the bins of exons provided in {} for analysis.".format(bin_gtf)
        else:
            if not exon_gtf:
                print "Error: An Exon gtf file needs to be provided to create bins."
                quit()
            print "Create bins for exons in the exon gtf {} for analysis.".format(exon_gtf)
            test_gtf = BinGtf(exon_gtf, number)
            test_gtf.ExtractGeneEndseq()
            test_gtf.generateGtfBin('Transcript.coverage.bin', size_cutoff)
            bin_gtf_file = 'Transcript.coverage.bin'
            command = 'cp {} ..'.format(bin_gtf_file)
            check_call(command, shell=True)
            print "The created bins are saved in {}.".format('Transcript.coverage.bin')
   
    multicov_ls = BedtoolsCoveragePathos(bin_gtf_file, bam, 'RNAseq') # multiple process on cluster or a single computer with multiple CPUS.
    test_count = CountCoverage(multicov_ls, number, output_file)
    test_count.CaculateCountInTranscript()
               
def CountCoverageEndseq(bam, exon_gtf, number, size_cutoff, output_file, bin_gtf):
    '''
    For Endseq data

    bam: 1bp bam file containing the primary reads.
    exon_gtf: An annotation gtf file containing exon information for each transcript.
    number: number of bins used for each transcript
    size_cutoff: 
    Only transcripts with length ≥ max(size_cutoff, number-1) will be used for RNA-seq coverage calculation.
    output_file: a matrix file can be used for plotting
    transcriptID, count in Bin1 ..., Count in Bin number-1
    
    bin_gtf: A gtf file containing N-1 bins for each transcript. Not required.
    e.g.
    chr1    .       .       11869   11884   .       +       .       ENST00000456328.2;bin1;
    chr1    .       .       11885   11901   .       +       .       ENST00000456328.2;bin2;
    chr1    .       .       11902   11917   .       +       .       ENST00000456328.2;bin3;
    chr1    .       .       11918   11934   .       +       .       ENST00000456328.2;bin4;
    chr1    .       .       11935   11950   .       +       .       ENST00000456328.2;bin5;
    '''
    print "Calculate the transcript body coverage for TSS data {}.".format(bam)
    print "TSS onebase bam file: {}".format(bam)
    print "The number of bins used: {}".format(number-1)
    print "Transcript size cutoff: {}".format(size_cutoff)
    print "The output matrix: {}".format(output_file)

    if bin_gtf: # provided gtf files with bins, save running time
        print "Use the bins of exons provided in {} for analysis.".format(bin_gtf)
        multicov_ls = BedtoolsCoveragePathos(bin_gtf, bam, 'TSS')
    else:
        if not exon_gtf:
            print "Error: An Exon gtf file needs to be provided to create bins."
            quit()
        print "Create bins for exons in the exon gtf {} for analysis.".format(exon_gtf)
        test_gtf = BinGtf(exon_gtf, number)
        test_gtf.ExtractGeneEndseq()
        test_gtf.generateGtfBin('Transcript.coverage.bin', size_cutoff) # Transcript.coverage.bin -> gencode.v24.exonForcoverage.100bin.gtf for 5'end seq
        multicov_ls = BedtoolsCoveragePathos('Transcript.coverage.bin', bam, 'TSS')
        command = 'cp {} ..'.format('Transcript.coverage.bin')
        check_call(command, shell=True)
        print "The created bins are saved in {}.".format('Transcript.coverage.bin')
    
    test_count = CountCoverage(multicov_ls, number, output_file)
    test_count.CaculateCountInTranscript()