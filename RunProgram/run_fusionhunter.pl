#! /usr/bin/env perl

use IO::File;

$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/fusionhunter-1.4/FusionHunter-v1.4/bin";


if (!defined($ENV{PERL5LIB}) || $ENV{PERL5LIB} !~ /\S/) {
    $ENV{PERL5LIB} = "/gne/research/workspace/twu/lib/perl5/site_perl/5.10.0";
} else {
    $ENV{PERL5LIB} = "/gne/research/workspace/twu/lib/perl5/site_perl/5.10.0:" . $ENV{PERL5LIB};
}


# Must be named FusionHunter.cfg
$FP = new IO::File(">FusionHunter.cfg") or die "Cannot write to file FusionHunter.cfg";
write_conf($FP,$fastq1,$fastq2);
close($FP);

system("$BINDIR/FusionHunter.pl FusionHunter.cfg > LOG");
system("rm -f FusionHunter.cfg");

$FP = new IO::File("FusionHunter.fusion");
while (defined($line = <$FP>)) {
    if ($line =~ /^\# Fusion/) {
	($clusterid) = $line =~ /\/(\S+)/;

	($strand1,$strand2) = $line =~ /\[(.)(.)\]/;
	$align = <$FP>;		# accession
	$align = <$FP>;		# alignment
	($gene1,$gene2) = $align =~ /(\S+) x (\S+)/;

	($chr1,$start1,$end1,$chr2,$start2,$end2) = $align =~ /(\S+):(\d+)-(\d+) (\S+):(\d+)-(\d+)/;
	if ($strand1 eq "-") {
	    $chrpos1 = $start1;
	} elsif ($strand1 eq "+") {
	    $chrpos1 = $end1;
	} else {
	    die "strand1 is $strand1.  line is $line";
	}
	if ($strand2 eq "-") {
	    $chrpos2 = $end2;
	} elsif ($strand2 eq "+") {
	    $chrpos2 = $start2;
	} else {
	    die "strand2 is $strand2.  line is $line";
	}

	print ">$clusterid $gene1..$gene2\n";
	print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
    }
}
close($FP);


exit;




sub write_conf {
    my ($FP, $fastq1, $fastq2) = @_;

    print $FP <<CONF;

########################
# version1.4
# modification 2011/06/09
##########################
# configuration file for FusionHunter
# '#' indicates comments
# '=' is in need for each configuration line
# edit this file in your working directory 
##########################


# hg19 annotation packages are added, replace all hg18 links with hg19 if you are using hg19



##########################
# Basic options, you must get them changed for your data!! 
##########################

#set 1 if running on hg18 or hg19; otherwise 0
IF_HUMAN = 1

# left part of read pairs, should be in fastq format
L = $fastq1

# right part of read pairs should be in fastq format
R = $fastq2

# reference dir/name, reference should be in fasta format
Reference =/gne/research/workspace/twu/fusion-bakeoff/hg19/hg19.fa

# the directory containing Bowtie index/basename of Bowtie index
# NO '/' in the end
#

BowtieIdx=/gne/research/workspace/twu/fusion-bakeoff/hg19/hg19

# the directory and name of gene annotation list, we suggest UCSC annotation, these AnnotationFiles are included in FusionHunter package
Gene_annotation = /gne/research/workspace/twu/fusion-bakeoff/fusionhunter-1.4/AnnotationFiles_hg19/hg19.ucscKnownGene

# directory and file name repeats region annotation, these AnnotationFiles are included in FusionHunter package
Repeats = /gne/research/workspace/twu/fusion-bakeoff/fusionhunter-1.4/AnnotationFiles_hg19/hg19.repeats

# the directory and file name of self alignment regions, these AnnotationFiles are included in FusionHunter package
SelfAlign = /gne/research/workspace/twu/fusion-bakeoff/fusionhunter-1.4/AnnotationFiles_hg19/hg19.chain.pairs

# the directory and file name of human EST database, these AnnotationFiles are included in FusionHunter package
EST = /gne/research/workspace/twu/fusion-bakeoff/fusionhunter-1.4/AnnotationFiles_hg19/hg19.SpliceEST

# size of the segmented reads, we strongly suggest it should not be longer than half of full read length e.g. <=25 if your RNA-seq read is 50bp
segment_size = 25

# number of cores for bowtie and other parallel processes, for sake of speed, we suggest you use as many cores as possible
CORE = 1

#  min number of paired-end reads that support a fusion (used in regionPairsList)
#  the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get, and shall take longer time
PAIRNUM = 2

# min number of junction spanning reads to support a fusion
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get
#
# TO BE NOTED: if you set the MINSPAN = 1, in order to reduce false positives, any candidate junction supported 
# by only 1 spanning read would be discarded unless the fusion junction point is exactly on annotated exon boundary. This process is embeded in FusionHunter.
MINSPAN = 1

# min size of the maximum base coverage on either side of the junction
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get
MINOVLP = 8

# Size of exact match for each junction flanking tile, should not be larger than MINOVLP
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get, and shall take longer time
# We strongly suggest it should be >=3
TILE = 4

# total mismatch
# the smaller this number, the more specificity you get; the larger this number, the more sensitivity you get, and shall take longer time
# We strongly suggest it should be <=4
MISMATCH = 2



#########################
# Advanced options, you may change them for your preference
# changes are optional
########################

# number of multi hits for partial reads
M1 = 20

# number of reads to keep for partial reads
K1 = 8

# number of multi hits for full reads
M2 = 1

# number of reads to keep for full reads
K2 = 1

# max allowed repeat proportion of a read (used in reduceBwt)
REAPTOVLP = 0.6

# number of chains to overlap with a read (used in reduceBwt)
CHAINNUM = 20

# max allowed repeat proportion of a read (more stringent, used in leftRightOvlp)
RPTOVLP = 0.2

# max allowed alignment proportion between a pair of reads (used in leftRightOvlp)
CHAINOVP = 0.2

# distance to self-chain boundary (used in postLeftRightOvlp)
CHAINDIS = 200000

# proportion of a overlaps with a region (used in regionPairs)
READOVLP = 0.8

# screen out repetitive regions of reference when doing gapped alignment
MASK = 1

# output of fusion by FusionHunter
fusion_output = FusionHunter.fusion

# output of readthrough by FusionHunter

readthrough_output = FusionHunter.readthrough

CONF

  return;
}

