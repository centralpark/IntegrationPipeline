#! /usr/bin/env perl

use warnings;
use IO::File;
use Getopt::Std;
undef $opt_i;			# input sample name
$opt_o = "sample";		# output sample name
undef $opt_d;			# genome
undef $opt_m;			# known genes
undef $opt_s;			# known splicesites
undef $opt_B;
undef $opt_g;			# gunzip
undef $opt_P;			# chr add prefix
getopts("i:o:d:m:s:B:gP:");

$GSTRUCTDIR = "/gne/research/apps/gstruct/gstruct-2013-05-22/x86_64-linux-2.6-sles11/bin";
$SAMPLE = $opt_o;
if (!defined($genome = $opt_d)) {
    $genome = "hg19_gamma";
}
if (!defined($genes = $opt_m)) {
    $genes = "hg19.known.genes";
}
if (!defined($splicesites = $opt_s)) {
    $splicesites = "hg19.known.splicesites";
}


$fastqs = join(" ",@ARGV);

# Run GSTRUCT-fusion
$command = "$GSTRUCTDIR/gstruct-fusion";
if (defined($opt_g)) {
    $command .= " -g";
}
if (defined($opt_B)) {
    $command .= " -B $opt_B";
}
## $command .= " -B sample.align";
if (defined($opt_P)) {
    $command .= " -P $opt_P";
}
if (defined($opt_i)) {
    $command .= " -i $opt_i";
}
$command .= " -d $genome -m $genes -o $SAMPLE -n 50 -s hg19.known.splicesites $fastqs";
print STDERR $command . "\n";
system($command);


if (0) {
# Process output (gstruct-fusion does this already)
    $FP = new IO::File("$SAMPLE.results.distinct");

    $prev_gene1 = "";
    $prev_gene2 = "";
    while (defined($line = <$FP>)) {
	if ($line =~ /^>/) {
	    chop $line;
	    ($nreads) = $line =~ /\((\d+)\)/;

	    @fields = split "-",$line;
	    $gene1 = $fields[2];
	    $gene2 = $fields[3];

	    ($info1,$info2) = $fields[1] =~ /(\S+)\.\.(\S+)/;
	    ($strand1,$chr1,$chrpos1) = $info1 =~ /([_+])(\S+)\@(\d+)/;
	    ($strand2,$chr2,$chrpos2) = $info2 =~ /([_+])(\S+)\@(\d+)/;
	    undef $strand1;
	    undef $strand2;

	    if ($gene1 ne $prev_gene1 || $gene2 ne $prev_gene2) {
		print ">$gene1..$gene2 $nreads\n";
		$prev_gene1 = $gene1;
		$prev_gene2 = $gene2;
	    }
	    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
	}
    }
}

exit;

