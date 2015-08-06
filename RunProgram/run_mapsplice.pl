#! /usr/bin/env perl

use IO::File;
use Getopt::Std;
undef $opt_G;			# do not use gene annotation
getopts("G");


$ENV{PATH} = "/gne/research/workspace/twu/fusion-bakeoff/mapsplice/MapSplice_multi_threads_2.0.1.6/bin:" . $ENV{PATH};

$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/mapsplice/MapSplice_multi_threads_2.0.1.6/bin";
$ANNOTDIR = "/gne/research/workspace/twu/fusion-bakeoff/mapsplice/annot";
$GENES = "$ANNOTDIR/Homo_sapiens.GRCh37.66.gtf";


$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

$command = "python $BINDIR/mapsplice_multi_thread.py -c $ANNOTDIR/hg19 -1 $fastq1 -2 $fastq2 -x $ANNOTDIR/hg19 --fusion";
if (defined($opt_G)) {
    # Do not use known genes
} else {
    $command .= " --gene-gtf $GENES";
}
$command .= " 2>/dev/null";
print STDERR $command . "\n";
system($command);


$FP = new IO::File("mapsplice_out/fusions_annotated.txt");
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $clusterid = $fields[3];

    ($chr1,$chr2) = $fields[0] =~ /(\S+)~(\S+)/;
    $chrpos1 = $fields[1];
    $chrpos2 = $fields[2];
    
    $genes1 = $fields[48];
    $genes1 =~ s/,$//;
    $genes1 =~ s/,/\|/g;
    $genes2 = $fields[49];
    $genes2 =~ s/,$//;
    $genes2 =~ s/,/\|/g;

    $encompass_count = $fields[27];

    print ">$clusterid $genes1..$genes2 $encompass_count\n";
    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
    $seenp{"$chr1:$chrpos1 $chr2:$chrpos2"} = 1;
}
close($FP);

# Handle novel splice sites
$FP = new IO::File("mapsplice_out/fusions_candidates.txt");
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $clusterid = $fields[3];

    ($chr1,$chr2) = $fields[0] =~ /(\S+)~(\S+)/;
    $chrpos1 = $fields[1];
    $chrpos2 = $fields[2];
    
    if (!defined($seenp{"$chr1:$chrpos1 $chr2:$chrpos2"})) {
	# $genes1 = $fields[48];
	# $genes1 =~ s/,$//;
	# $genes1 =~ s/,/\|/g;
	# $genes2 = $fields[49];
	# $genes2 =~ s/,$//;
	# $genes2 =~ s/,/\|/g;

	$encompass_count = $fields[27];

	print ">$clusterid Novel..Novel $encompass_count\n";
	print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
    }
}
close($FP);

exit;


