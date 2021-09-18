# python2.7
# -*- coding: utf-8 -*-
'''
Based on python 2.7

Use to prepare the exon annotation used for Transcript body coverage counting.

gtf: gencode.v24.exon.gtf

chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; tag "basic"; transcript_support_level "1"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA  exon    12613   12721   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; tag "basic"; transcript_support_level "1"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
..
chr1    HAVANA  exon    29534   29570   .       -       .       gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P-001"; exon_number 1; exon_id "ENSE00001890219.1"; level 2; ont "PGO:0000005"; tag "basic"; transcript_support_level "NA"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1    HAVANA  exon    24738   24891   .       -       .       gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P-001"; exon_number 2; exon_id "ENSE00003507205.1"; level 2; ont "PGO:0000005"; tag "basic"; transcript_support_level "NA"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1    HAVANA  exon    18268   18366   .       -       .       gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P-001"; exon_number 3; exon_id "ENSE00003477500.1"; level 2; ont "PGO:0000005"; tag "basic"; transcript_support_level "NA"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
Note:
For both positive and negative strand, exon number 1 is always the most 5'end exon;
exons are ordered based on exon number.
See CoverageSortForCount.pptx to explain how to sort the coverage file
'''

from collections import OrderedDict
import re

def transcriptID(x):
    return re.findall('transcript_id \"(.*?)\";',x)[0]


def PrepareGencode(gtf, output_file):
    '''
    gtf: gencode.v24.exon.gtf
    output_file: gencode.v24.exonForcoverage.gtf

    Use to extract exons of transcript from gencode exon gtf annoation and calculate the length of transcript
    
    output gencode.v24.exonForcoverage.gtf e.g.:
    <chr><.><.><transcript start><transcript end><id; transcript size><strand><.><exon_number 1 start; exon_number 1 end; exon_number 2 start; exon_number 2 end ..>
    chr1    .       .       11869   14409   ENST00000456328.2;1657; +       .       11869;12227;12613;12721;13221;14409;
    chr1    .       .       12010   13670   ENST00000450305.2;632;  +       .       12010;12057;12179;12227;12613;12697;12975;13052;13221;13374;13453;13670;
    chr1    .       .       14404   29570   ENST00000488147.1;1351; -       .       29534;29570;24738;24891;18268;18366;17915;18061;17606;17742;17233;17368;16858;17055;16607;16765;15796;15947;15005;15038;14404;14501;
    '''
    print "Prepare the exon annotation used by count function."
    dic_transcript = OrderedDict() # {(chr, strand, id) : [(exon_number 1 start, exon_number 1 end) ...]}
    with open(gtf) as f:
        for line in f:
            id, strand, chr = transcriptID(line.strip().split('\t')[8]), line.strip().split('\t')[6], line.strip().split('\t')[0]
            if dic_transcript.has_key((chr, strand, id)):
                dic_transcript[(chr, strand, id)].append((line.strip().split('\t')[3], line.strip().split('\t')[4]))
            else:
                dic_transcript[(chr, strand, id)] = [(line.strip().split('\t')[3], line.strip().split('\t')[4])]

    output = open(output_file, 'w')
    for item in dic_transcript:
        exon_boundary = [int(x) for sublist in dic_transcript[item] for x in sublist] # convert [('11869', '12227'), ('12613', '12721'), ('13221', '14409')] to flat list ['11869', '12227', '12613', '12721', '13221', '14409']
        # <chr><.><.><transcript_start><transcript_end><transcript_id, transcript_size><strand><.><(exon1 start, exon1 end), (exon2 start, exon2 end)...>
        exon_size = [int(exon[1])-int(exon[0])+1 for exon in dic_transcript[item]]
        infor = [item[2], str(sum(exon_size)), '']
        exon_ls = map(str, exon_boundary)
        exon_ls.append('')
        temp = [item[0], '.', '.', str(min(exon_boundary)), str(max(exon_boundary)), ';'.join(infor), item[1], '.', ';'.join(exon_ls)]
        print>>output, '\t'.join(temp)
    output.close()
    return 1
