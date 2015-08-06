#! /usr/bin/env perl

use IO::File;
use Getopt::Std;
undef $opt_S;			# summarize only
getopts("S");


$ENV{PYTHONPATH} = "/gne/research/workspace/twu/lib64/python2.6/site-packages";
$ENV{PATH} = "/gne/research/workspace/twu/bin:" . $ENV{PATH};

$BINDIR = "/gne/research/workspace/twu/bin";
$INDEXDIR = "/gne/research/workspace/twu/fusion-bakeoff/chimerascan/hg19";
$OUTPUTDIR = ".";

$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

if (!defined($opt_S)) {
    $command = "python $BINDIR/chimerascan_run.py -v --quals solexa $INDEXDIR $fastq1 $fastq2 $OUTPUTDIR";
    system($command);
}


$FP = new IO::File("$OUTPUTDIR/chimeras.bedpe");
if (defined($header = <$FP>)) {
    # Skip
}
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $clusterid = $fields[6];
    $gene1 = $fields[12];
    $gene2 = $fields[13];


    $chr1 = $fields[0];
    $chr2 = $fields[3];
    if ($fields[8] eq "+") {
	$chrpos1 = $fields[2] + 1;
    } elsif ($fields[8] eq "-") {
	$chrpos1 = $fields[1] + 1;
    } else {
	die "Cannot parse +/- in field 7 of $OUTPUTDIR/chimeras.bedpe";
    }
    if ($fields[9] eq "+") {
	$chrpos2 = $fields[4] + 1;
    } elsif ($fields[9] eq "-") {
	$chrpos2 = $fields[5] + 1;
    } else {
	die "Cannot parse +/- in field 7 of $OUTPUTDIR/chimeras.bedpe";
    }

    $nfrags = $fields[16];
    $n_spanning_frags = $fields[17];
    @breakpoint_spanning_frags = split ",",$fields[21];

    printf ">$clusterid $gene1..$gene2 $nfrags $n_spanning_frags %d\n",$#breakpoint_spanning_frags + 1;
    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
}
close($FP);


exit;


