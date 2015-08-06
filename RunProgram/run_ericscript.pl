#! /usr/bin/env perl

use IO::File;

$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/ericscript-0.4.0";

$ENV{'PATH'} = $ENV{'PATH'} . ":/gne/research/workspace/twu/bin";

$FASTQ1 = $ARGV[0];
if (!defined($FASTQ2 = $ARGV[1])) {
    $FASTQ2 = "";
}

if (-e output) {
    system("rm -rf output");
}

# Took out --remove, so we can debug
$command = "$BINDIR/ericscript.pl --genomeref=/gne/research/workspace/twu/fusion-bakeoff/hg19/hg19.fa --refid=homo_sapiens --outputfolder=output --samplename=sample $FASTQ1 $FASTQ2";
system($command);


# Process results

$FP = new IO::File("output/sample.results.filtered.tsv") or die "Cannot find results file";

if (defined($header = <$FP>)) {
    undef $header;
}

while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $gene1 = $fields[0];
    $gene2 = $fields[1];
    $chr1 = "chr" . $fields[2];
    $chrpos1 = $fields[3];
    $chr2 = "chr" . $fields[5];
    $chrpos2 = $fields[6];

    $fusion = "$gene1..$gene2";
    $crossing_reads{$fusion} = $fields[10];
    $spanning_reads{$fusion} = $fields[11];

    push @ {$coords{$fusion}},"$chr1\t$chrpos1\t$chr2\t$chrpos2";
}

close($FP);


# Print results

foreach $fusion (keys %coords) {
    print ">$fusion";
    print " $crossing_reads{$fusion}";
    print " $spanning_reads{$fusion}\n";
    foreach $coord (@ {$coords{$fusion}}) {
	print $coord . "\n";
    }
}
    
exit;

