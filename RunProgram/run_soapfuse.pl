#! /usr/bin/env perl


use IO::File;


$BINDIR = "/gne/research/workspace/twu/fusion-bakeoff/soapfuse/SOAPfuse-v1.25";
$CONFIG = "$BINDIR/config/config.txt";


$fastq1 = $ARGV[0];
$fastq2 = $ARGV[1];

if (1) {
system("mkdir -p A/Lib-a");
system("gzip -c $fastq1 > A/Lib-a/Run-a_1.fq.gz");
system("gzip -c $fastq2 > A/Lib-a/Run-a_2.fq.gz");

$FP = new IO::File(">input.txt") or die "Cannot write to file input.txt";
print $FP "A\tLib-a\tRun-a\t75\n";
close($FP);


$FP = new IO::File(">config.txt") or die "Cannot write to file config.txt";
write_conf($FP);
close($FP);

system("$BINDIR/SOAPfuse-RUN.pl -c config.txt -fd . -l input.txt -o .");
}

parse_results("./final_fusion_genes/A/A.final.Fusion.specific.for.genes");

system("rm -rf A");

exit;


sub parse_results {
    my ($filename) = @_;

    $FP = new IO::File($filename) or die "Cannot open file $filename";
    if (defined($header = <$FP>)) {
    }
    while (defined($line = <$FP>)) {
	chop $line;
	@fields = split /\t/,$line;
	$gene1 = $fields[0];
	$gene2 = $fields[5];
	$chr1 = $fields[1];
	$chrpos1 = $fields[3];
	$chr2 = $fields[6];
	$chrpos2 = $fields[8];

	$nspanning{"$gene1..$gene2"} += $fields[10];
	$njunction{"$gene1..$gene2"} += $fields[11];
	push @ {$coords{"$gene1..$gene2"}},"$chr1\t$chrpos1\t$chr2\t$chrpos2";
    }
    close($FP);

    foreach $genepair (keys %coords) {
	print ">$genepair $nspanning{$genepair} $njunction{$genepair}\n";
	foreach $coord (@ {$coords{$genepair}}) {
	    print $coord . "\n";
	}
    }

    return;
}



sub write_conf {
    my ($FP) = @_;

    print $FP <<CONF;

DB_db_dir  =   /gne/research/workspace/twu/fusion-bakeoff/soapfuse/data

#----------------------------------#
##   database files of SOAPfuse   ##
#----------------------------------#

DB_wg_soap_ref       =     \$(db_dir)/WG_index_soap/human.fa.index
DB_trans_soap_ref    =     \$(db_dir)/transcript_index_soap/transcript.fa.index
DB_trans_bwa_ref     =     \$(db_dir)/transcript_index_bwa/transcript.fa
DB_trans_psl         =     \$(db_dir)/transcript.psl
DB_trans_gtf         =     \$(db_dir)/Homo_sapiens.gtf
DB_gene_psl          =     \$(db_dir)/gene.psl
DB_gene_fa           =     \$(db_dir)/gene.fa
DB_genefamily        =     \$(db_dir)/gene_family_HGNC/gene_family_HGNC.brief.txt
DB_blast_homo_list   =     \$(db_dir)/blast_homo_gene.m8.gz

######################################
######     All programs used   #######
######################################

#----------------------------------#
##   directory of all programs    ##
#----------------------------------#

PG_pg_dir   =    /gne/research/workspace/twu/fusion-bakeoff/soapfuse/SOAPfuse-v1.25/source/bin

#----------------------------------#
##  programs in SOAPfuse pipeline ##
#----------------------------------#

PG_soap               =     \$(pg_dir)/aln_bin/soap2.21
PG_bwa                =     \$(pg_dir)/aln_bin/bwa
PG_blat               =     \$(pg_dir)/aln_bin/blat
PG_bwt                =     \$(pg_dir)/aln_bin/2bwt-builder2.20
PG_convert            =     \$(pg_dir)/aln_bin/convert

####################################################
#### All perl scripts used in SOAPfuse pipeline ####
####################################################

#----------------------------------#
##   directory of perl scripts    ##
#----------------------------------#

PS_ps_dir   =   /gne/research/workspace/twu/fusion-bakeoff/soapfuse/SOAPfuse-v1.25/source

#----------------------------------#
##   perl scripts of each step    ##
#----------------------------------#

PS_s01            =     \$(ps_dir)/SOAPfuse-01-alignWG.pl
PS_s02            =     \$(ps_dir)/SOAPfuse-02-align_unmap_transcript.pl
PS_s03            =     \$(ps_dir)/SOAPfuse-03-align_trim_unmap_transcript.pl
PS_s04            =     \$(ps_dir)/SOAPfuse-04-change_SE.pl
PS_s05            =     \$(ps_dir)/SOAPfuse-05-candidate.pl
PS_s06            =     \$(ps_dir)/SOAPfuse-06-divide_soap_denovo_unmap.pl
PS_s07            =     \$(ps_dir)/SOAPfuse-07-junction_seq_deal.pl
PS_s08            =     \$(ps_dir)/SOAPfuse-08-final_fusionGene.pl
PS_s09            =     \$(ps_dir)/SOAPfuse-09-deeper_analysis.pl

######################################
######   The output directory   ######
######################################

#----------------------------------------------------------------#
## the directory which will stores sub-directories of each step ##
#----------------------------------------------------------------#

# else you can set it when run SOAPfuse, in its command line options. (perl SOAPfuse-RUN.pl -o PD_all_out)
PD_all_out           =     .

#--------------------------------#
## sub-directories of each step ##
#--------------------------------#

# SOAPfuse will make (mkdir) them automatically

PD_alignWG                   =     \$(all_out)/alignWG
PD_align_unmap_Tran          =     \$(all_out)/align_unmap_Tran
PD_align_trim_unmap_Tran     =     \$(all_out)/align_trim_unmap_Tran
PD_change_SE                 =     \$(all_out)/change_SE
PD_candidate                 =     \$(all_out)/candidate
PD_denovo_unmap              =     \$(all_out)/denovo_unmap
PD_junction_seq              =     \$(all_out)/junction_seq
PD_final_fusion_genes        =     \$(all_out)/final_fusion_genes

##############################################
#### Parameters used in SOAPfuse pipeline ####
##############################################

#----------------------------------#
##    SOAPfuse operation mode     ##
#----------------------------------#

# operation fusion detection in the SOMATIC mode.
# default is 'no', set 'yes' to enable it.
# 'no' means only detect fusions for single sample, not in tumor-vs-control mode.
# SOAPfuse can seek the somatic fusions, you must write both Normal and Tumor info in a single sample-list.
PA_all_somatic_mode        =     no

# the postfixes of sample-ID to distinguish Normal and Tumor samples, effective when the SOMATIC mode is enabled.
# please prepare all strings that tag after the patient-ID.
# e.g., like 'K101-N' and 'K101-CA', '-N' is the postfix of Normal tissue, and '-CA' is the Tumor's.
# please following the default format, add 'N:' or 'T:' before the postfixes to specify the tissue.
# if you have multiple postfixes for one type tissue (Normal or Tumor), just write them with ';' as separator.
# we have prepared several common postfixes of sample-ID in cancer research.
# NOTE: SOAPfuse use perl regular module to distinguish the postfix in a case-insensitive manner.
PA_all_postfix_of_tissue   =     'N:-N;N:-Normal;N:-B;N:-Blood;T:-CA;T:-C;T:-T;T:-Tumor;T:-Cancer'

#----------------------------------#
##  Global parameters in pipeline ##
#----------------------------------#

# the postfix of fq/fa file name.
# default is 'fq.gz'. [ the fq/fa file must be named as {Lane_name}_[12].(postfix) ]
PA_all_fq_postfix     =     fq.gz

# the cpu process used when do alignment.
# default is 8.
PA_all_process_of_align_software     =     8

# the shortest length that unmapped reads will be trimmed to.
# default is 30.
# this length should be >= 30, or will be adjusted to 30.
# set 'skip' to avoid trim operation.
# SOAPfuse will do the trimming operation automatically when your RNA-Seq data has short insert size
# that leads to the overlapping of two ends from the same pair.
PA_all_shortest_length_trim_unmap_to     =     40

# maximum number of genes that the trimmed reads mapped.
# default is 1.
# SOAPfuse trims the unmapped reads, and realigns them to seek more alignments with reduced uniqueness.
# set this parameter to take the non-uniquely mapped trimmed reads into account as long as the number of
# their mapped genes is not larger than setting.
PA_all_maximum_genome_loc_trimmed_read_mapped     =     2

# maximum number of genes that intact reads mapped.
# default is 1.
PA_all_maximum_genome_loc_intact_read_mapped     =     1

# the length extended at the exon edge into intron for seek the junction-site in intron.
# default is 0, means not seek any junction-site in intron.
# if set as non-zero, more computational resource may need, esp. running-time of step s07.
# steps involved by this parameter are s07, s08 and s09.
# if you want seek junction-site in intron, set it as 100 (suggested), controlled by the maximum 200.
PA_all_intron_len_extend_from_exon_edge     =     default

#------------------------------------------------#
##  parameters about realigning unmapped reads  ##
#------------------------------------------------#

# realign unmapped reads to make them more clean.
# default is 'no', set 'yes' to enable it.
# this realignment uses bwa to filter unmap-read caused by samll indels.
PA_s02_realign     =     yes

#----------------------------------------------#
##  parameters about filtering candidate set  ##
#----------------------------------------------#

# sign to save genes whose name has the dot character ('.').
# default is 'no', set 'yes' to enable it.
# generally, we donot care these genes for their less research values as unspecific functions.
PA_s05_save_genes_name_with_dot     =     default

# sign to save genes from same gene family.
# default is 'no', set 'yes' to enable it.
# The gene family data is get from 'http://www.genenames.org/genefamily.html'.
# if you think this database is not ok for your research, set this parameter as 'yes'.
PA_s05_save_genes_from_same_family     =     default

# sign to enable the amass control of span-reads, set 'no' to disable it.
# It is suggestted that enable it ('yes') when you want to get a result with low FP.
# Once enabled, the memory and time the following steps consume will decrease apparently.
PA_s05_amass_control_of_span_reads     =     yes

# the minimum number of paired span-reads that support the fusion gene pairs.
# default is 1.
# used when filter the fusion candidate set.
# the larger this parameter is, the less time and memory the following steps will consume.
PA_s05_the_minimum_span_reads_for_candidate     =     2

#------------------------------------------------------------#
##  parameters about dealing devided-denovo unmapped reads  ##
#------------------------------------------------------------#

# sign to save the reads which have mismatch at the end.
# default is 'no', set 'yes' to enable it.
# if 'yes', some false junc-reads may remain, leading to some false positives.
PA_s06_save_reads_have_mismatch_around_fusepos     =     default

#%%%% followed two parameters are effective when 'PA_s06_save_reads_have_mismatch_around_fusepos' is 'no' %%%%#

# the flank region near the end of mapped half-reads for filtering mismatch.
# default is 5.
PA_s06_number_of_flank_bases_near_read_end_for_filter_mismatch     =     default

# the maximum number of mismatch allowed in flank region when filter.
# default is 0.
PA_s06_the_maximum_mismatch_in_flank_region     =     default

#-------------------------------------------------------------------------#
##  parameters about construct junction library for junc-reads alignment ##
#-------------------------------------------------------------------------#

# the minimum number of span-reads when picks gene pairs from candidate set for junction library construction.
# default is 2.
PA_s07_the_minimum_span_reads_for_junction_construction     =     2

# the number of extended bases around the ends of paired-end reads when calculate the region of exhaustion.
# default is 0.
# the larger this parameter is, the more computation amount of local_exhaustion is, but not obviously.
PA_s07_extended_bases_near_pe_read_end     =     default

# the minimum consistency used for define the credible regions to construct fusion junction sequence.
# default is 0 [Maximum is 0.8].
# e.g., 0.8, means that bases are considered as useful only when number of span-reads that cover them
# is >= 80% of all span-reads that support the same fusion gene pairs.
# if you want use all regions that span-reads support, set this parameter as 0.
# the larger this parameter is, the less time and memory SOAPfuse will use in the following operations.
# But large parameter may cause the losing of some low-frequency fusion events.
PA_s07_the_min_cons_for_credible_fuse_region     =     0.5

# maximum number of mismatch allowed in alignment to putative junction suquences library.
# default is 3.
PA_s07_maximum_mismatch_for_align_junction_reads     =     default

# the number of flank bases for checking mismatch around junction point.
# default is 5, which means [fusepos-5bp,fusepos+5bp].
PA_s07_flank_bases_around_fuse_point_for_check_mismatch     =     default

# maximum number of mismatch allowed in the flank region for mismatch checking.
# default is 0.
PA_s07_maximum_mismatch_in_flank_region     =     default

# the number of bases that a junc-read must map to both sides of the fusion at least.
# default is >=5.
# proper number will make the junc-reads more credible.
# Large number may lose some fusions supported by a small quantity of junc-reads.
PA_s07_junc_read_map_both_sides_at_least     =     7

#--------------------------------------------#
##  parameters about get final fusions      ##
#--------------------------------------------#

# the number of bases extended at the inner ends of span-reads.
# default is 0.
# It is suggested to set it as non-zero only when your RNA-Seq data has short insert size.
PA_s08_number_of_extend_bases     =     default

# use the insert size to select the most proper supporting.
# default is 'no', set 'yes' to enable it.
# Set it as 'yes' when your data has good library construction (thin distribution of insert size, small SD).
PA_s08_insert_control_sup     =   default

# the minimum sum number of junc_read and span-read of one fusion event.
# default is 2.
PA_s08_min_sum_reads     =     3

#--------------- e.g. 'PA_s08_min_support_reads_for_both_edge' is '1,2' -------------------#
#-  it means only when span_reads >= 1 and junc_read >=2, then remain edge-edge fusions   -#
#------   edge: edge of exon;    internal: internal part of exon or intron region    ------#
# all three parameters default as '1,1'.
PA_s08_min_support_reads_for_both_edge     =     1,1
PA_s08_min_support_reads_for_one_edge_one_internal     =     2,2
PA_s08_min_support_reads_for_both_internal     =     2,2

# the minimum distance between two intra-chromosome gene pairs.
# the gene pair are from same strand and fused in their original genomic orentation.
# obviously, this parameter is set for read-through type fusion.
# default is 1000.
PA_s08_min_intrachr_distance     =     default

# the minimum bases covered at both up (5' end) and down (3' end) sides around the junction point.
# default is 5.
PA_s08_min_bases_covered_both_sides_around_fuse_point     =     10

# only remain the fusion events happened at edge of exon.
# default is 'no', set 'yes' to enable it.
PA_s08_only_remain_edge_case     =     default

#------------------------------------------------------#
##  parameters about analysis based on final fusions  ##
#------------------------------------------------------#

# draw SOAPfuse fusion figures for all final fusion events.
# default is 'no', set 'yes' to enable it.
PA_s09_draw_fusion_expression_svg     =     yes

CONF

    return;
}

