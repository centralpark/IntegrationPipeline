#! /usr/bin/env perl


use IO::File;

$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

$ENV{PYTHONPATH} = "/gne/research/workspace/twu/lib64/python2.6/site-packages";

$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/shortfuse/ShortFuse-0.2";


$working_dir = "working_dir";

if (1) {
system("mkdir -p $working_dir");

if ($fastq1 =~ /\//) {
    ($root) = $fastq1 =~ /.*\/([^\/]+)$/;
} else {
    $root = $fastq1;
}
system("touch $working_dir/$root.discord.bust");

$FP = new IO::File(">shortfuse.conf") or die "Cannot write to file shortfuse.conf";
write_conf($FP,$working_dir);
close($FP);

# Error: Previously had --conf $$.conf
system("python $BINDIR/run_pipeline.py --conf shortfuse.conf $fastq1 $fastq2 > LOG");
#system("rm -f shortfuse.conf");
}


$FP = new IO::File("$working_dir/fusion_counts.bedpe");
if (defined($header = <$FP>)) {
    # Skip
}


while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $gene1 = $fields[10];
    $gene2 = $fields[11];

    $chr1 = $fields[0];
    $chr2 = $fields[3];
    $chrpos1_low = $fields[1] + 1;
    $chrpos1_high = $fields[2];

    $chrpos2_low = $fields[4] + 1;
    $chrpos2_high = $fields[5] + 1;

    push @ {$coords{"$gene1..$gene2"}},"$chr1\t$chrpos1_low\t$chr2\t$chrpos2_low";
    push @ {$coords{"$gene1..$gene2"}},"$chr1\t$chrpos1_low\t$chr2\t$chrpos2_high";
    push @ {$coords{"$gene1..$gene2"}},"$chr1\t$chrpos1_high\t$chr2\t$chrpos2_low";
    push @ {$coords{"$gene1..$gene2"}},"$chr1\t$chrpos1_high\t$chr2\t$chrpos2_high";
}
close($FP);

foreach $genepair (keys %coords) {
    print ">$genepair\n";
    print join("\n",@ {$coords{$genepair}}) . "\n";
}

exit;




sub write_conf {
    my ($FP, $working_dir) = @_;

    print $FP <<CONF;

#The directory where Bowtie is installed. It should have bowtie and bowtie-build binaries in it
bowtie_dir                  /gne/research/workspace/twu/bin

#The directory where ShortFuse is installed. It should have the EM and sift binaries in it
#along with the Python files included in the download.
shortfuse_dir               /gne/research/workspace/twu/fusion-bakeoff/shortfuse/ShortFuse-0.2

#The directory containing the ShortFuse reference files. It should have the
#refseq_transcripts and transcripts_plus_genome Bowtie index files in it.
ref_dir                     /gne/research/workspace/twu/fusion-bakeoff/shortfuse/ShortFuse_ref

#The directory where the intermediate files and results from ShortFuse will be written.
#If it does not exist, it will be created.
working_dir                 $working_dir

#The number of processes to use when calling Bowtie, the -p flag. It's a good idea for
#this to be one less than the total number of available processors.
processes                   1

#Whether or not to overwrite files if they already exist.
#If this is set to True, ShortFuse will perform steps even if the expected output
#files already exist. If False, it will skip those steps.
overwrite_existing          True

#How to call Python 2.7 . This probably doesn't need to be changed, but if you
#have a local installation of python or your python points to Python 3, you
#will need to change this.
python_bin                  python

CONF

  return;
}

