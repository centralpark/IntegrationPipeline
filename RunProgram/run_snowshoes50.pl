#! /usr/bin/env perl

use IO::File;
use Getopt::Std;
$opt_l = 50;			# read length (50 or 75)
getopts("l:");

if (!defined($readlength = $opt_l)) {
    die "Need to specify readlength with -l flag";
}

# For megablast and samtools
$ENV{'PATH'} = $ENV{'PATH'} . ":/gne/research/workspace/twu/bin";
$ENV{'PATH'} = $ENV{'PATH'} . ":/gne/research/apps/gne-sequence-analysis/bin/x86_64-linux-2.6-sles11";


$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/PERL_scripts";

# Needs to be version 0.5.10 or earlier.  Versions 0.6.1 and 0.6.2 fail with given indices
$BWA = "/gne/research/workspace/twu/bin/bwa";
$HG19 = "/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/ref_fusion/BWA_indexed_genome_build37/allchr_hg19.fa";

if ($readlength == 75) {
    $JUNCTIONS = "/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/ref_fusion/BWA_indexed_Junction_build37/hg19_refFlat_filtered.75.junction";
} elsif ($readlength = 50) {
    $JUNCTIONS = "/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/ref_fusion/BWA_indexed_Junction_build37/hg19_refFlat_filtered.50.junction";
} else {
    die "readlength must be 50 or 75";
}

# Must be named sample_infor.txt
$FP = new IO::File(">sample_infor.txt") or die "Cannot write to file sample_infor.txt";
write_sample_info($FP);
close($FP);


# Must be named configure_file.txt
$FP = new IO::File(">configure_file.txt") or die "Cannot write to file configure_file.txt";
write_conf($FP);
close($FP);


# Must be named s_1_1_sequence and s_1_2_sequence
$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

$command = "ln -s $fastq1 s_1_1_sequence";
system($command);
$command = "ln -s $fastq2 s_1_2_sequence";
system($command);
$fastq1 = "s_1_1_sequence";
$fastq2 = "s_1_2_sequence";


$command = "$BINDIR/0_make_directory_structures -p .";
system($command);

# Align reads
$command = "$BWA aln -l 32 -t 4 $HG19 $fastq1 > BWA_aln_genome/s_1_aln_1.sai";
system($command);
$command = "$BWA aln -l 32 -t 4 $HG19 $fastq2 > BWA_aln_genome/s_1_aln_2.sai";
system($command);
$command = "$BWA aln -l 32 -t 4 $JUNCTIONS $fastq1 > BWA_aln_Junction/s_1_aln_1.sai";
system($command);
$command = "$BWA aln -l 32 -t 4 $JUNCTIONS $fastq2 > BWA_aln_Junction/s_1_aln_2.sai";
system($command);

# Make SAM files
$command = "$BWA samse $HG19 BWA_aln_genome/s_1_aln_1.sai $fastq1 > BWA_aln_genome/s_1_1.sam";
system($command);
$command = "$BWA samse $HG19 BWA_aln_genome/s_1_aln_2.sai $fastq2 > BWA_aln_genome/s_1_2.sam";
system($command);
$command = "$BWA samse $HG19 BWA_aln_Junction/s_1_aln_1.sai $fastq1 > BWA_aln_Junction/s_1_1.sam";
system($command);
$command = "$BWA samse $HG19 BWA_aln_Junction/s_1_aln_2.sai $fastq2 > BWA_aln_Junction/s_1_2.sam";
system($command);

# Make sorted SAM files
$command = "samtools view -Sb BWA_aln_genome/s_1_1.sam -o BWA_aln_genome/s_1_1.bam";
system($command);
$command = "samtools view -Sb BWA_aln_genome/s_1_2.sam -o BWA_aln_genome/s_1_2.bam";
system($command);
$command = "samtools sort -n -m 4000000000 BWA_aln_genome/s_1_1.bam BWA_aln_genome/s_1_1_ID_sorted";
system($command);
$command = "samtools sort -n -m 4000000000 BWA_aln_genome/s_1_2.bam BWA_aln_genome/s_1_2_ID_sorted";
system($command);
$command = "samtools view -X BWA_aln_genome/s_1_1_ID_sorted.bam -o BWA_aln_genome/sorted_sam/s_1_1_ID_sorted.sam";
system($command);
$command = "samtools view -X BWA_aln_genome/s_1_2_ID_sorted.bam -o BWA_aln_genome/sorted_sam/s_1_2_ID_sorted.sam";
system($command);

$command = "samtools view -Sb BWA_aln_Junction/s_1_1.sam -o BWA_aln_Junction/s_1_1.bam";
system($command);
$command = "samtools view -Sb BWA_aln_Junction/s_1_2.sam -o BWA_aln_Junction/s_1_2.bam";
system($command);
$command = "samtools sort -n -m 4000000000 BWA_aln_Junction/s_1_1.bam BWA_aln_Junction/s_1_1_ID_sorted";
system($command);
$command = "samtools sort -n -m 4000000000 BWA_aln_Junction/s_1_2.bam BWA_aln_Junction/s_1_2_ID_sorted";
system($command);
$command = "samtools view -X BWA_aln_Junction/s_1_1_ID_sorted.bam -o BWA_aln_Junction/sorted_sam/s_1_1_ID_sorted.sam";
system($command);
$command = "samtools view -X BWA_aln_Junction/s_1_2_ID_sorted.bam -o BWA_aln_Junction/sorted_sam/s_1_2_ID_sorted.sam";
system($command);

# Analyze alignments for fusions
$command = "$BINDIR/1_process_sorted_SAM_to_mapped_reads -p . > LOG";
system($command);
$command = "$BINDIR/2_Mapped_Reads_processed -p . >> LOG";
system($command);
$command = "$BINDIR/3_get_fusion_transcript_and_protein_results -p . >> LOG";
system($command);


# Look at results in results/final_fusion_report_RNA.txt
$FP = new IO::File("results/final_fusion_report_RNA.txt");
if (defined($line = <$FP>)) {
    # Skip header
}
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    ($gene1,$gene2) = $fields[2] =~ /(\S+)->(\S+)/;
    @infoa = split ":",$fields[10];
    @infob = split ":",$fields[11];
    $genea = $infoa[2];
    $geneb = $infob[2];
    if ($genea eq $gene1 && $geneb eq $gene2) {
	@info1 = split ":",$fields[10];
	@info2 = split ":",$fields[11];
    } elsif ($geneb eq $gene1 && $genea eq $gene2) {
	@info1 = split ":",$fields[11];
	@info2 = split ":",$fields[10];
    } else {
	die "$genea/$geneb do not match $gene1/$gene2.  info is $fields[10]/$fields[11] line is $line";
    }

    $chr1 = $info1[1];
    $strand1 = $info1[6];
    if ($strand1 eq "-") {
	$chrpos1 = $info1[4] + 1;
    } elsif ($strand1 eq "+") {
	$chrpos1 = $info1[5];
    } else {
	die "strand1 is $strand1.  line is $line";
    }

    $chr2 = $info2[1];
    $strand2 = $info2[6];
    if ($strand2 eq "-") {
	$chrpos2 = $info2[5];
    } elsif ($strand2 eq "+") {
	$chrpos2 = $info2[4] + 1;
    } else {
	die "strand2 is $strand2.  line is $line";
    }

    print ">$gene1..$gene2\n";
    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";

}
close($FP);


exit;


sub write_sample_info {
    my ($FP) = @_;

    print $FP <<SAMPLE;
Sample ID	Lab ID	Sample ID
1	sim	simulated-data
SAMPLE

    return;
}


sub write_conf {
    my ($FP) = @_;

    print $FP <<CONF;
################################################################################################################
## NOTE: please don't change the format of this file. The only things to change 
## are the values after the "=" sign. please also don't insert any space before or after "=" in the line;
################################################################################################################
## NOTE: \$snowshoes_home is where the perl scripts are
\$snowshoes_home=/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/PERL_scripts
## NOTE: \$refdata is where all the reference data are. All content downoaded from the ref_fusion folder should be here.
\$refdata=/gne/research/workspace/twu/fusion-bakeoff/snowshoes/SnowShoes-FTD_2.0_Build37/ref_fusion
## NOTE: \$read_length is the read length of reads
\$read_length=$readlength
## NOTE: \$distance is the the minimum distance between two genes for fusion candidates if they are on the same chromosome? (e.g: type in 5000 for 5kb). We recommend >5000
\$distance=50000
## NOTE: \$lib_size is the size of the RNA fragments during mRNA-Seq library construction. The default is 300 bp
\$lib_size=200
## NOTE: is the minimal number of read pairs that supports each fusion gene candidate. We recommend a minimal value of 10
\$minimal=5
## NOTE: \$max_fusion_isoform is the maximum number of fusion isoforms between two fusion parnters. The recommanded value is 2
\$max_fusion_isoform=5
CONF

    return;
}
