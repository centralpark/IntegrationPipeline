#! /usr/bin/env perl

use IO::File;
use Getopt::Std;
$opt_p = 0.50;			# Cutoff for probability
getopts("p:");


$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/defuse-0.5/defuse-0.5.0/scripts";
$CONFIG = "$BINDIR/config.txt";


# Needs to have .1.fastq and .2.fastq
$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

if (0) {
    ($datadir,$filename1) = $fastq1 =~ /(\S+)\/(\S+)/;
    ($datadir,$filename2) = $fastq2 =~ /(\S+)\/(\S+)/;
    if ($filename1 =~ /(\S+)\.1\.fastq/) {
	$root = $1;
    } else {
	die "Cannot determine fileroot from $filename1";
    }

} else {
    $command = "ln -s $fastq1 raw.1.fastq";
    system($command);
    $command = "ln -s $fastq2 raw.2.fastq";
    system($command);
    $datadir = `pwd`;
    chop $datadir;
    $root = "reads";
    $filename1 = "reads.1.fastq";
    $filename2 = "reads.2.fastq";
}



my $reads_index_filename = $root.".fqi";
my $reads_names_filename = $root.".names";
my $reads_sources_filename = $root.".sources";

# Requires that retreive_fastq.pl be modified to rely upon this file
$FP = new IO::File(">$reads_sources_filename");
print $FP "raw.1.fastq\n";
print $FP "raw.2.fastq\n";
close($FP);


$command = "$BINDIR/retreive_fastq.pl -c $CONFIG -d $datadir -1 $filename1 -2 $filename2 -i $reads_index_filename -n $reads_names_filename -s $reads_sources_filename -f > LOG";
print STDERR $command . "\n";
system($command);

# Do not call with -d flag, because already called retreive_fastq.pl directly above
$command = "$BINDIR/defuse.pl -c $CONFIG -o . >> LOG";
print STDERR $command . "\n";
system($command);


$FP = new IO::File("results.classify.tsv");
if (defined($line = <$FP>)) {
    # Skip header
}
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $clusterid = $fields[0];
    $gene1 = $fields[31];
    $gene2 = $fields[32];

    $chr1 = "chr" . $fields[25];
    $chrpos1 = $fields[38];
    $chr2 = "chr" . $fields[26];
    $chrpos2 = $fields[39];
    $prob = $fields[65];
    if ($prob < $opt_p) {
	print "#>$clusterid $gene1..$gene2 $prob\n";
	print "#$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
    } else {
	print ">$clusterid $gene1..$gene2 $prob\n";
	print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
    }
}
close($FP);

exit;



