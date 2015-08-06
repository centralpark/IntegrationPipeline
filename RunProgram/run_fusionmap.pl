#! /usr/bin/env perl

use IO::File;

$MONO = "/gne/research/workspace/twu/bin/mono";
$EXE = "/gne/research/workspace/twu/fusion-bakeoff/fusionmap-2013-02-01/FusionMap/bin/FusionMap.exe";
$ANNOTDIR = "/gne/research/workspace/twu/fusion-bakeoff/fusionmap-2013-02-01/annot";

$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];


$FP = new IO::File(">secontrol.config") or die "Cannot write to file secontrol.config";
write_conf($FP);
close($FP);


$command = "$MONO $EXE --semap $ANNOTDIR Human.B37 RefGene secontrol.config 2>/dev/null";
print STDERR $command . "\n";
system($command);

$FP = new IO::File("output/TestDataset_InputPEReads.FusionReport.txt");
if (defined($line = <$FP>)) {
    # Skip header
}
while (defined($line = <$FP>)) {
    chop $line;
    @fields = split /\t/,$line;
    $chr1 = "chr" . $fields[5];
    $chrpos1 = $fields[6];
    $chr2 = "chr" . $fields[7];
    $chrpos2 = $fields[8];
    $gene1 = $fields[9];
    $gene2 = $fields[13];
    print ">$gene1..$gene2\n";
    print "$chr1\t$chrpos1\t$chr2\t$chrpos2\n";
}
close($FP);

exit;





# Note: Use32BitMode must be set to True, otherwise program crashes early
# But even this fix still gives a crash
sub write_conf {
    my ($FP) = @_;

    print $FP <<CONF;
<Files>
$fastq1
$fastq2

<Options>
QualityEncoding=Sanger
MonoPath=/gne/research/workspace/twu/bin/mono
PairedEnd=True
RnaMode=True
Use32BitMode=False
ThreadNumber=8 //Possible values: 1-100. Default value=1
FileFormat=FASTQ // Possible values: FASTQ, QSEQ, FASTA. Default value=FASTQ
MinimalFusionAlignmentLength=25 //Possible values: 15-50. Default value=25 (alpha)
FusionReportCutoff=1 // Possible values: 1-1000. Default value=1 (beta)
NonCanonicalSpliceJunctionPenalty=4 //Possible values: 01-. Default value = 2 (G)
MinimalHit=2 // Minimal distinct read; Possible values: 1-10000, Default value =2 
MinimalRescuedReadNumber=1 // Minimal rescued read number. Default value = 1
OutputFusionReads=True // Possible values: True, False. Default value = True
FilterBy=DefaultList //advanced filtering using default black list from FusionMap, set to None to avoid automatic downloading

<Output>
TempPath=/tmp
OutputPath=output
OutputName=TestDataset_InputPEReads
CONF

    return;
}

