#! /usr/bin/env perl

use IO::Dir;
use IO::File;


$BOWTIE = "/gne/research/workspace/twu/bin/bowtie";
$BOWTIE_BUILD = "/gne/research/workspace/twu/bin/bowtie-build";
$BINDIR = "/gne/research/workspace/twu/bin";
$ANNOTDIR = "/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/FusionSeq_Data_hg19";
$BOWTIE_INDICES = "/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/Bowtie_Indices";
$RSEQTOOLS = "/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/tusharb/RSEQtools-0.6";

# Need to have blat and bowtie be in current path
$ENV{PATH} = "/gne/research/workspace/twu/bin:" . $ENV{PATH};


$FP = new IO::File(">fusionseqrc") or die "Cannot write to file fusionseqrc";
write_conf($FP);
close($FP);
$ENV{FUSIONSEQ_CONFPATH} = "fusionseqrc";


$ENV{LD_LIBRARY_PATH} = "/gne/research/apps/gsl/1.15/x86_64-linux-2.6-sles11/lib:" . $ENV{LD_LIBRARY_PATH};


$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];


# Format reads
system("mkdir -p data");
$FP1 = new IO::File($fastq1) or die;
$FP2 = new IO::File($fastq2) or die;
$OUT = new IO::File(">data/interleaved.fq") or die;
$READS = new IO::File(">data_allReads.txt") or die;
interleave($OUT,$READS,$FP1,$FP2);
close($READS);
close($OUT);
close($FP2);
close($FP1);

system("gzip --force data_allReads.txt");


# Align reads with Bowtie
# -p is nthreads
# Junctions do not align paired
#system("$BOWTIE -q $BOWTIE_INDICES/hg19/hg19 -1 $fastq1 -2 $fastq2 > sim1_hg19.bowtie.paired");
system("$BOWTIE -q $BOWTIE_INDICES/hg19/hg19 data/interleaved.fq > sim1_hg19.bowtie.single");
#system("$BOWTIE -q $BOWTIE_INDICES/knownGene_2x60_spliceJunctions data/interleaved.fq > sim1_jcns.bowtie.single");

# Make MRF files (but only unpaired singles appear to contribute to data.1.gfr)
#system("$RSEQTOOLS/bowtie2mrf paired -sequence -IDs < sim1_hg19.bowtie.paired > paired.mrf");
system("$RSEQTOOLS/bowtie2mrf genomic -sequence -IDs < sim1_hg19.bowtie.single > single.mrf");
#system("$RSEQTOOLS/bowtie2mrf junctions -sequence -IDs < sim1_jcns.bowtie.single > junctions.mrf");

$IN = new IO::File("single.mrf") or die;
$OUT = new IO::File(">unpaired.mrf") or die;
merge_single_mrf($OUT,$IN);
close($OUT);
close($IN);


########################################################################

# Step 1
system("$BINDIR/geneFusions data 4 < unpaired.mrf > data.1.gfr");


# Step 2
# Have to create an empty blackList.txt file
# Need to have blat and bowtie be in the current path, so 
system("rm -f blackList.txt");
system("touch blackList.txt");
system("$BINDIR/gfrMitochondrialFilter < data.1.gfr | $BINDIR/gfrRepeatMaskerFilter $ANNOTDIR/hg19_repeatMasker.interval 5 | $BINDIR/gfrCountPairTypes | $BINDIR/gfrExpressionConsistencyFilter | $BINDIR/gfrAbnormalInsertSizeFilter 0.01 | $BINDIR/gfrPCRFilter 4 4 | $BINDIR/gfrProximityFilter 1000 | $BINDIR/gfrAddInfo | $BINDIR/gfrAnnotationConsistencyFilter ribosomal | $BINDIR/gfrAnnotationConsistencyFilter pseudogenes | $BINDIR/gfrBlackListFilter blackList.txt | $BINDIR/gfrLargeScaleHomologyFilter | $BINDIR/gfrRibosomalFilter | $BINDIR/gfrSmallScaleHomologyFilter > data.gfr");

# Step 3
# Have to create data.meta manually
$nreads = `grep -v AlignmentBlocks unpaired.mrf | grep -v "#" | wc -l`;
chop $nreads;
$FP = new IO::File(">data.meta");
print $FP "Mapped_reads\t$nreads\n";

system("$BINDIR/gfrConfidenceValues data < data.gfr > data.confidence.gfr");


# Step 4
# For 50-mers, 40 means 10 nt mapped to either of the tiles
# For 75-mers, use 65

system("$BINDIR/gfr2bpJunctions data.confidence.gfr 65 200 1.0");


# Step 5: Process joblist 1
system("mkdir -p Bowtie_Indices");

$job_name = "job.$$";
$FP = new IO::File("data_jobList1.txt");
while (defined($line = <$FP>)) {
    chop $line;
    $line =~ s/bowtie-build/$BOWTIE_BUILD/;
    $line =~ s/bowtie /$BOWTIE /;
    $line =~ s/$BOWTIE_INDICES/Bowtie_Indices/g;    # Put Bowtie indices locally
    $command = "bsub -o /dev/null -J $job_name \"$line\"";
    print STDERR $command . "\n";
    system($command);
}
close($FP);


# Wait for all jobs to finish
$command = "bsub -o /dev/null -K -w 'done($job_name)' sleep 1";
system($command);
print STDERR "Waiting for all jobs to finish...";


# Step 6: Process joblist 2
$FP = new IO::File("data_jobList2.txt");
while (defined($line = <$FP>)) {
    chop $line;
    system($line);
}
close($FP);

$DP = new IO::Dir(".");
while (defined($file = $DP->read)) {
    if ($file =~ /(\S+)\.bp$/) {
	$root = $1;
	if ($root =~ /filtered/) {
	    # Skip
	} else {
	    $command = "$BINDIR/validateBpJunctions < $file | $BINDIR/bpFilter 4 4 100 0.01 30 > $root.filtered.bp";
	    print STDERR $command . "\n";
	    system($command);

	    if (-z "$root.filtered.bp") {
		# Skip empty files, which cause bp2alignment to crash
	    } else {
		$command = "$BINDIR/bp2alignment < $root.filtered.bp > $root.alignments.txt";
		print STDERR $command . "\n";
		system($command);
	    
		parse_alignments($root,"$root.alignments.txt");

#	    $command = "$BINDIR/bp2wig $root.filtered.bp";
#	    print STDERR $command . "\n";
#	    system($command);
	    }
	}
    }
}
close($DP);


exit;


sub interleave {
    my ($OUT, $READS, $FP1, $FP2) = @_;
    undef $line1;
    undef $line2;

    while (defined($entry1 = read_entry_1($FP1)) && defined($entry2 = read_entry_2($FP2))) {
	print $OUT join("\n",@ {$entry1}) . "\n";
	print $OUT join("\n",@ {$entry2}) . "\n";
	print $READS $ {$entry1}[1] . "\n";
	print $READS $ {$entry2}[1] . "\n";
    }
    return;
}

sub read_entry_1 {
    my ($FP1) = @_;
    @entry1 = ();

    if (!defined($line1)) {
	if (!defined($line1 = <$FP1>)) {
	    return;
	} else {
	    chop $line1;
	}
    }
    push @entry1,$line1;

    while (defined($line1 = <$FP1>) && $line1 !~ /^@/) {
	chop $line1;
	push @entry1,$line1;
    }
    if (defined($line1)) {
	chop $line1;
    }

    return \@entry1;
}

sub read_entry_2 {
    my ($FP2) = @_;
    @entry2 = ();

    if (!defined($line2)) {
	if (!defined($line2 = <$FP2>)) {
	    return;
	} else {
	    chop $line2;
	}
    }
    push @entry2,$line2;

    while (defined($line2 = <$FP2>) && $line2 !~ /^@/) {
	chop $line2;
	push @entry2,$line2;
    }
    if (defined($line2)) {
	chop $line2;
    }

    return \@entry2;
}

sub merge_single_mrf {
    my ($OUT, $IN) = @_;

    if (defined($header = <$IN>)) {
	print $OUT $header;
    }

    while (defined($line = <$IN>)) {
	chop $line;
	@fields = split /\t/,$line;
	$block = $fields[0];
	$sequence = $fields[1];
	$acc = $fields[2];
	$acc =~ s/\/\d$//;
	if (defined($lastacc) && $acc eq $lastacc) {
	    print $OUT $lastblock . "|" . $block;
	    print $OUT "\t";
	    print $OUT $lastsequence . "|" . $sequence;
	    print $OUT "\t";
	    print $OUT $lastacc . "|" . $acc;
	    print $OUT "\n";
	    undef $lastacc;
	}
	$lastblock = $block;
	$lastsequence = $sequence;
	$lastacc = $acc;
    }

    if (defined($lastacc) && $acc eq $lastacc) {
	print $OUT $lastblock . "|" . $block;
	print $OUT "\t";
	print $OUT $lastsequence . "|" . $sequence;
	print $OUT "\t";
	print $OUT $lastacc . "|" . $acc;
	print $OUT "\n";
    }

    return;
}


sub parse_alignments {
    my ($root, $alignment_file) = @_;
    my $FP;
    my $line;
    my $printp = 0;

    $FP = new IO::File($alignment_file);
    while (defined($line = <$FP>)) {
	if ($line =~ /Tile 1: (\S+):\d+-(\d+)/) {
	    $chr1 = $1;
	    $chrpos1 = $2;
	    $line = <$FP>;
	    if ($line =~ /Tile 2: (\S+):(\d+)-\d+/) {
		$chr2 = $1;
		$chrpos2 = $2;
		if ($printp == 0) {
		    print ">$root\n";
		    $printp = 1;
		}
		print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
	    }
	}
    }
    close($FP);
    return;
}


sub write_conf {
    my ($FP) = @_;

    print $FP <<CONF;
# --------------------------------- This section is required ---------------------------------
# Location of the bowtie indexes of the human genome and the composite model
BOWTIE_INDEXES="/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/Bowtie_Indices" 
# the subdirectory of BOWTIE_INDEXES where the human genome is indexed by bowtie-build
BOWTIE_GENOME="hg19"
# the subdirectory of BOWTIE_INDEXES where the composite model is indexed by bowtie-build
BOWTIE_COMPOSITE="hg19_knownGeneAnnotationTranscriptCompositeModel"

# Pointer to the program twoBitToFa part of the blat suite
BLAT_TWO_BIT_TO_FA="/gne/research/workspace/twu/bin/twoBitToFa" 
# Location and filename of the reference genome in 2bit format (to be used by blat)
BLAT_DATA_DIR="/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/hg19"
BLAT_TWO_BIT_DATA_FILENAME="hg19.2bit"

# Location and name of the transcript composite model sequence and interval files
TRANSCRIPT_COMPOSITE_MODEL_DIR="/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/FusionSeq_Data_hg19"
TRANSCRIPT_COMPOSITE_MODEL_FA_FILENAME="knownGeneAnnotationTranscriptCompositeModel.fa"
TRANSCRIPT_COMPOSITE_MODEL_FILENAME="knownGeneAnnotationTranscriptCompositeModel.txt"

# location of the annotation files
ANNOTATION_DIR="/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/FusionSeq_Data_hg19"
# conversion of knownGenes to gene symbols, description, etc. 
KNOWN_GENE_XREF_FILENAME="kgXref.txt" 
# conversion of knownGenes to TreeFam
KNOWN_GENE_TREE_FAM_FILENAME="knownToTreefam.txt" 

# Location and filename of the ribosomal library
RIBOSOMAL_DIR="/gne/research/workspace/twu/fusion-bakeoff/fusionseq-0.7/FusionSeq_Data_hg19"
RIBOSOMAL_FILENAME="ribosomal.2bit"

# ----------------------- This section is optional: visualization tools -------------------------
# URL of the cgi directory on the web server
WEB_URL_CGI="http://cgiURL" 
# location of the data directory on the web server, as seen from the web server
WEB_DATA_DIR="/path/to/data" 
# URL of the data directory on the web server
WEB_DATA_LINK="http://dataURL" 
# Number of nucleotides flanking the region (for UCSC Genome Browser)
UCSC_GENOME_BROWSER_FLANKING_REGION=500 
# URL of the public website (non cgi)
WEB_PUB_DIR="http://publicURL" 
# Location of the structural data for Circos
WEB_SDATA_DIR="/path/to/structural/Data/Circos" 
# Location of Circos installation
WEB_CIRCOS_DIR="/path/to/circos" 
CONF

    return;
}

