#! /usr/bin/env perl

use IO::File;

$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

# For bowtie
$ENV{PATH} = "/gne/research/workspace/twu/bin:" . $ENV{PATH};

$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/tophat2/tophat-2.0.0.Linux_x86_64";
$INDEXDIR = "/gne/research/workspace/twu/fusion-bakeoff/hg19/hg19";

system("ln -s /gne/research/workspace/twu/fusion-bakeoff/tophat2/annot/refGene.txt .");
system("ln -s /gne/research/workspace/twu/fusion-bakeoff/tophat2/annot/ensGene.txt .");

# Output dir must contain the string "tophat_", because tophat-fusion-post has 'string.find(dir,"tophat_")'
$command1 = "$BINDIR/tophat -o tophat_output -p 8 --fusion-search --keep-fasta-order --bowtie1";
$command1 .= " --no-coverage-search -r 50 --mate-std-dev 80 --fusion-min-dist 100000";
$command1 .= " --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM";
$command1 .= " $INDEXDIR $fastq1 $fastq2";

$command2 = "$BINDIR/tophat-fusion-post -p 8 --num-fusion-reads 1 --num-fusion-pairs 2";
$command2 .= " --num-fusion-both 5 --skip-blast $INDEXDIR";

print STDERR $command1 . "\n";
#system($command1);

print STDERR $command2 . "\n";
#system($command2);


$FP = new IO::File("tophatfusion_out/result.txt");
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $gene1 = $fields[1];
    $gene2 = $fields[4];

    $chr1 = $fields[2];
    $chrpos1 = $fields[3];
    $chr2 = $fields[5];
    $chrpos2 = $fields[6];

    print ">$gene1..$gene2\n";
    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
}
close($FP);

exit;

