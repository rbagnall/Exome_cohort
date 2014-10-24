#!/usr/bin/perl

#This script aligns a cohort of paired end NGS data to the illumina TruSeq or Nextera Exome or sureselect_EZ_Exome_v4 + UTR (+50 nucleotides of padding either side of exome intervals)
#written by richard bagnall (r.bagnall@centenary.org.au)
#Usage: Exome_cohort_v.pl -fastq [path/2/user] -cohort [cohort_name] -exome [truseq or sureselectv4]
#save raw reads as [name1].1.fastq.gz, [name1].2.fastq.gz , [name2].1.fastq.gz , [name2].2.fastq.gz, etc...
#save the raw reads in a directory called Rawdata e.g. /home/shared/NGS/human/[USER]/Rawdata

use strict; use warnings;
use File::Copy;
use Getopt::Long;
use Parallel::ForkManager;
use List::Util qw( min max );

####################################################
# Setting up command line options and usage errors #
####################################################

# set a command line option variable (-cohort); is required and used to assign a filename to the cohort's VCF file in unified_genotyper subroutine
my $cohort = ''; 
# set a command line option variable (-fastq); is required and is the path to the rawdata
my $input ='';
# set a command line option variable (-exome); is required and is the exome enrichment file used (either truseq for the illumina TruSeq exome kit, or sureselectv4 for the Sureselect_Exomev4+UTR kit)
my $exome ='';

GetOptions ("cohort=s" => \$cohort,
			"fastq=s" => \$input,
            "exome=s" => \$exome);

my $offset = length($input); # use this extensively in subroutines for getting names of samples, files and folders

# -cohort is required, else print usage and die
if ($cohort eq '') {
	print "\n\n\t*** ERROR: You need to define a name for this cohort with the -cohort option\n";
	usage();
    exit;
}
# -fastq is required, else print usage and die
if ($input eq '') {
	print "\n\n\t*** ERROR: You need to define the path to the raw fastq files with the -fastq option\n";
	usage();
    exit;
}
# pre-empt common error
elsif ($input =~ m/\/$/) {
	print "\n\n\t*** ERROR: Please remove the last / from the -fastq option\n";
	usage();
    exit;
}
# -exome is required, else print usage and die
if ($exome !~ m/^(truseq|sureselectV4)$/) {
	print "\n\n\t*** ERROR: You need to define the exome enrichment kit used with the -exome option\n";
	usage();
    exit;
}

print "\n\n\n\t\t---------EXOME SEQUENCING ANALYSIS PIPELINE v1.2.1-----------\n";
print "\t\tAnalysing $cohort cohort exome data \n";
print "\t\tAligning Paired-end sequence reads (Illumina 1.8+)\n";
print "\t\tUsing the $exome exome enrichment kit\n";
print "\t\tUsing BWA MEM for read alignment\n";
print "\t\tUsing Picard to remove duplicate reads\n";
print "\t\tUsing GATK v3 for realigning\n";
print "\t\tUsing GATK v3 for recalibration\n";
print "\t\tUsing GATK v3 for Unified Genotyping\n";
print "\t\tUsing Hard filters and VQSR\n\n";
print "\t\tUsing SeattleSeq Annotation website\n\n";
print "\t\t-------------------------------------------------------------\n\n";
print "\n";

# create temporary folders
mkdir ("$input/tempBAMs") or die "Unable to create tempBAMs directory: <$!>\n";
mkdir ("$input/tempVCFs") or die "Unable to create tempVCFs directory: <$!>\n";

# paths
my $path2rawdata = "$input/Rawdata";
my $path2bam = "$input/tempBAMs";
my $path2vcf = "$input/tempVCFs";
my $path2gatk = '/home/groups/cardio/Applications/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar';
my $path2indel1kg = '/home/groups/cardio/References/INDELS/1000G_phase1.indels.b37.vcf';
my $path2indelmills = '/home/groups/cardio/References/INDELS/Mills_and_1000G_gold_standard.indels.b37.sites.vcf';
my $path2mergedexome50 = '/home/groups/cardio/References/Exome/Exome_seq_pipeline/merged_plus50.intervals';
my $path2mergedexome10 = '/home/groups/cardio/References/Exome/Exome_seq_pipeline/merged_plus10.intervals';
my $path2targetregions = "/home/groups/cardio/References/Exome/Exome_seq_pipeline/$exome.targetregions.bed";
my $path2SNP =   '/home/groups/cardio/References/SNP/dbsnp135.b37.vcf';
my $path2omni = '/home/groups/cardio/References/SNP/1000G_omni2.5.b37.sites.vcf';
my $path2hapmap = '/home/groups/cardio/References/SNP/hapmap_3.3.b37.sites.vcf';
my $path2ref = '/home/groups/cardio/References/Bwa_b37/hs37d5.fa';

#####################
#   Start BWA mem   #
#####################

my @fq = glob("$path2rawdata/*.fastq.gz"); # make array of fastq.gz

# check that there are the same number of forward and reverse files, and they have correct extensions
my @f_fq = grep(/1.fastq.gz$/i, @fq);
my @r_fq = grep( /2.fastq.gz$/i, @fq );
if (scalar(@f_fq) != scalar(@r_fq)) {
    rmdir "$path2bam";
    rmdir "$path2vcf";
    die "*** error: There must be a 1.fastq.gz and 2.fastq.gz file for each sample\n\n";
}

my $eightfork_manager = Parallel::ForkManager->new(3);

for (my $i = 0; $i < @fq; $i = $i ++) {
    my @fq_pair = splice(@fq, $i, 2);
    $eightfork_manager->start and next;
    bwa_mem(join(" ", @fq_pair));
    $eightfork_manager->finish;
}

$eightfork_manager-> wait_all_children;

print "\n\n*** BWA mapping complete ***\n";
clock();

#######################
#  Realign Bamfile    #
#######################

my $sixfork_manager = Parallel::ForkManager->new(2);

my @sort_ddbamfiles = glob("$path2bam/*.sorted.bam"); #make array of dedupped.sorted.bam files
for (my $i = 0; $i < @sort_ddbamfiles; $i++) {
    $sixfork_manager->start and next;
	indexbam($sort_ddbamfiles[$i]); # loop through array of bamfiles and pass to indexbam subroutine
	realigner_target_creator($sort_ddbamfiles[$i]); # pass to realigner_target_creator subroutine
	indel_realigner($sort_ddbamfiles[$i]); # pass to indel_realigner subroutine
    unlink ($sort_ddbamfiles[$i]); # cleanup: delete bam files
    unlink ("$sort_ddbamfiles[$i].bai"); # cleanup: delete bam index files
    $sixfork_manager->finish;
}

$sixfork_manager-> wait_all_children;

my @indel_intervals = glob("$path2bam/*.indel.intervals"); # loop through array of .intervals files and delete them
for (my $i = 0; $i < @indel_intervals; $i++) { unlink ("$indel_intervals[$i]") }


print "\n\n*** GATK indel realigner complete ***\n";
clock();

#####################################
# base quality score recalibration  #
#####################################

my $fstsixfork_manager = Parallel::ForkManager->new(2);

my @realigned = glob("$path2bam/*realigned.bam"); #make array of .realigned.bam files
for (my $i = 0; $i < @realigned; $i++) { # loop through array of bamfiles
    $fstsixfork_manager->start and next;
	base_recalibration($realigned[$i]); # pass to base_recalibration subroutine
    $fstsixfork_manager->finish;
}

$fstsixfork_manager-> wait_all_children;

my @grp = glob("$path2bam/*.grp"); # loop through array of .grp files and delete them
for (my $i = 0; $i < @grp; $i++) { unlink ("$grp[$i]") }

print "\n\n*** GATK BQSR complete ***\n";
clock();

#################################
# Depth of coverage calculation #
#################################

my $sndsixfork_manager = Parallel::ForkManager->new(2);

my @recalibrated = glob("$path2bam/*recalibrated.bam"); #make array of .recalibrated.bam files
for (my $i = 0; $i < @recalibrated; $i++) { # loop through array of bamfiles
    $sndsixfork_manager->start and next;
	coverage($recalibrated[$i]); # pass to coverage sub
    $sndsixfork_manager->finish;
}

$sndsixfork_manager-> wait_all_children;

sleep (10);

print "\n\n*** Depth of coverage calculation complete ***\n";
clock();

sleep (10);

# prompt user for permission to continue with genotyping, or quit

yn();

#####################
# Unified Genotyper #
#####################

# make a string of the recalibrated.bam files, separated by -I
my @recalibrated_bam_files = glob("$path2bam/*recalibrated.bam"); # make array of recalibrated bam files
my $recalibrated_joined = join(" -I ",@recalibrated_bam_files); # $input/tempBAMs/sample1.reduced.bam -I $input/tempBAMs/sample2.reduced.bam
unified_genotyper($recalibrated_joined);

###################
# Select Variants #
###################

my @raw_vcf = glob("$path2vcf/*raw.vcf"); # make array of raw vcf files (even though there is only 1 file)

for (my $i = 0; $i < @raw_vcf; $i++) {
my $sample_name = substr $raw_vcf[$i], ($offset + 10), -8; # get [samplename]
print "\n\n*** Selecting $sample_name indel variations ***\n";
clock();
system("java -jar $path2gatk -T SelectVariants -R $path2ref -V $raw_vcf[$i] -o $path2vcf/$sample_name.indels.vcf -selectType INDEL");
print "\n\n*** Selecting $sample_name SNP variations ***\n";
clock();
system("java -jar $path2gatk -T SelectVariants -R $path2ref -V $raw_vcf[$i] -o $path2vcf/$sample_name.snps.vcf -selectType SNP");
}

######################
# Variant Filtration #
######################

# filter indels first: make array of indel.vcf
my @raw_indel_vcf = glob("$path2vcf/*indels.vcf");

for (my $i = 0; $i < @raw_indel_vcf; $i++) {
my $sample_name = substr $raw_indel_vcf[$i], ($offset + 10), -11; # get [samplename]
print "\n\n*** Filtering $sample_name indel variations for quality ***\n";
clock();
system("java -jar $path2gatk -T VariantFiltration -R $path2ref -V $raw_indel_vcf[$i] -o $path2vcf/$sample_name.filtered.indels.vcf --filterExpression 'FS > 200.0 ' --filterName 'FS'  --filterExpression 'QD < 1.8 ' --filterName 'QD' --filterExpression 'ReadPosRankSum < -20.0' --filterName 'RPRS' --filterExpression 'InbreedingCoeff < -0.8' --filterName 'IC' ");
# cleanup: delete unfiltered indel and indel index files
unlink ($raw_indel_vcf[$i]);
unlink ("$raw_indel_vcf[$i].idx");
}

# filter snps: make array of snp.vcf
my @raw_snps_vcf = glob("$path2vcf/*snps.vcf");
for (my $i = 0; $i < @raw_snps_vcf; $i++) {
my $sample_name = substr $raw_snps_vcf[$i], ($offset + 10), -9; # get [samplename]
print "\n\n*** Filtering $sample_name SNP variations for quality ***\n";
clock();
system("java -jar $path2gatk -T VariantFiltration -R $path2ref -V $raw_snps_vcf[$i] -o $path2vcf/$sample_name.filtered.snps.vcf --filterExpression 'FS > 60.0 ' --filterName 'FS'  --filterExpression 'QD < 2.0' --filterName 'QD' --filterExpression 'ReadPosRankSum < -8.0' --filterName 'RPRS' --filterExpression 'HaplotypeScore >13.0' --filterName 'HS' --filterExpression 'MQ < 25.0' --filterName 'MQ' --filterExpression 'MQRankSum <-15.0' --filterName 'MQRS' --mask $path2vcf/$sample_name.filtered.indels.vcf --maskName 'INDEL' ");

unlink ($raw_snps_vcf[$i]); # cleanup: Delete unfiltered snp and index files
unlink ("$raw_snps_vcf[$i].idx");
}

#############################
#   Move files and cleanup  #
#############################

print "\n\n*** Renaming files and folders ***\n";
clock();

cleanup();

#################################
#   Begin Annotation of snps    #
#################################

print "\n\n*** Starting snp annotation pipeline ***\n";
clock ();

mkdir ("/home/groups/cardio/.do_not_delete_this_folder") or die "Unable to create do_not_delete_this_folder directory <$!>\n";
# make test.vcf with headers for submission to seattleseq annotation
open(TEST, ">/home/groups/cardio/.do_not_delete_this_folder/test.vcf") or die "Could not open test.vcf: <$!>\n";
print TEST "\# autoFile testAuto.txt\n\# compressAuto true\n"; # need these two header lines for submission
# get headers and PASS sites from vcf file
system("grep '^#\\|PASS' < $input/$cohort/Variants/$cohort.cohort.filtered.snps.vcf >>/home/groups/cardio/.do_not_delete_this_folder/test.vcf");
# load test.vcf to seattleseq
chdir('/home/groups/cardio/Applications/perl_codes'); # need to run from this directory [bloody java!]
system("java SubmitSeattleSeqAnnotationAutoJob");
system("gunzip /home/groups/cardio/.do_not_delete_this_folder/test.vcf.txt.gz");

# now have an file called text.vcf.txt and need to get it back to sample name
my $current_file_name = "/home/groups/cardio/.do_not_delete_this_folder/test.vcf.txt";
my $new_file_extension = "$input/$cohort/Variants/$cohort.snp.annotated.txt";
move($current_file_name, $new_file_extension) or die "Unable to move $current_file_name to $new_file_extension: $!\n";
close TEST;
# cleanup
unlink("/home/groups/cardio/.do_not_delete_this_folder/test.vcf");
rmdir "/home/groups/cardio/.do_not_delete_this_folder" or die "Unable to delete do_not_delete_this_folder directory: <$!>\n";

###################################
#   Begin secondary annotations   #
###################################



print STDERR "\n\n***Starting secondary annotation of snp variants ***\n";
clock();

# make a hash of ESP6500 variants
print STDERR "Retrieving ESP6500 variation frequencies\n";
open(ESP6500, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/ESP6500_all_snps.anno") or die "Can't open ESP6500.anno: <$!>\n";
my %ESP6500_hash;
for (<ESP6500>) {
    chomp;
    my($esp_position,$esp_freq)=split(/@/); #splits the row into elements
    my $ESP6500_key="$esp_position"; #puts first three columns of anno file into one key
    my $ESP6500_value="$esp_freq"; #puts rs number EA and AA freq into value
    $ESP6500_hash{$ESP6500_key} = $ESP6500_value; #sets the key of the hash to the value
}
close ESP6500;
print STDERR "Done\n\n";

# make a hash of ExAC variants
print STDERR "Retrieving ExAC variation counts\n";
open(ExAC, "/home/shared/NGS/human/richardb/Annotations/ExAC.r0.1.sites.snp.anno") or die "Can't open ExAC.anno: <$!>\n";
my %ExAC_hash;
for (<ExAC>) {
    chomp;
    my($exac_chr,$exac_pos,$exac_ref,$exac_alt,$exac_count)=split(/\t/); #splits the row into elements
    my $ExAC_key="$exac_chr\t$exac_pos\t$exac_ref\t$exac_alt"; #puts chrom, position, ref and alt of anno file into one key
    my $ExAC_value="$exac_count"; #puts allele count into value
    $ExAC_hash{$ExAC_key} = $ExAC_value; #sets the key of the hash to the value
}
close ExAC;
print STDERR "Done\n\n";

# make a hash of 1KG variants
print STDERR "Retrieving 1000 genomes allele frequencies\n";
open(KG_2011, "/home/shared/NGS/human/richardb/Annotations/ALL.autosomes.phase3.snp.anno") or die "Can't open 1KG.anno: <$!>\n";
my %KG2011_hash;
for (<KG_2011>) {
    chomp;
    my($kg_chr,$kg_pos,$kg_ref,$kg_alt,$kg_AFfreq,$kg_AMR,$kg_EAS,$kg_SAS,$kg_AFR,$kg_EUR)=split(/\t/); #splits the row into elements
    my $kg_key="$kg_chr\t$kg_pos\t$kg_ref\t$kg_alt"; #puts chr position ref and alt into one key
    my $kg_value="$kg_AFfreq\t$kg_AMR\t$kg_EAS\t$kg_SAS\t$kg_AFR\t$kg_EUR"; #puts frequences into value
    $KG2011_hash{$kg_key} = $kg_value; #sets the key of the hash to the value
}
close KG_2011;
print STDERR "Done\n\n";

# make a hash of gene info
print STDERR "Retrieving gene names and descriptions\n";
open(GENE_INFO, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/gene_info.anno") or die "Can't open GENE_INFO.anno: <$!>\n";
my %geneinfo_hash;
for (<GENE_INFO>) {
    chomp;
    my($geneinfo_position,$gene_info)=split(/@/); #splits the row into elements (chr plus gene name, and fkpm)
    my $geneinfo_key="$geneinfo_position\t"; #puts chr and gene name and tab into one key
    my $geneinfo_value="$gene_info"; #puts gene info into value
    $geneinfo_hash{$geneinfo_key} = $geneinfo_value; #sets the key of the hash to the value
}
close GENE_INFO;
print STDERR "Done\n\n";

# make a hash of heart expresssion
print STDERR "Retrieving gene expression values for heart tissue\n";
open(HEART_EXP, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/humanheart_fkpm.anno") or die "Can't open HEART_EXP.anno: <$!>\n";
my %heartexp_hash;
for (<HEART_EXP>) {
    chomp;
    my($heartexp_position,$heartexp_info)=split(/@/); #splits the row into elements
    my $heartexp_key="$heartexp_position\t"; #puts chr and gene name and tab into one key
    my $heartexp_value="$heartexp_info"; #puts gene fkpm into value
    $heartexp_hash{$heartexp_key} = $heartexp_value; #sets the key of the hash to the value
}
close HEART_EXP;
print STDERR "Done\n\n";

# make a hash of aortic valve expresssion
print STDERR "Retrieving gene expression values for aortic valve tissue\n";
open(AOVALVE_EXP, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/humanaorticvalve_fkpm.anno") or die "Can't open AOVALVE.anno: <$!>\n";
my %aovalveexp_hash;
for (<AOVALVE_EXP>) {
    chomp;
    my($aovalveexp_position,$aovalveexp_info)=split(/@/); #splits the row into elements
    my $aovalveexp_key="$aovalveexp_position\t"; #puts chr and gene name and tab into one key
    my $aovalveexp_value="$aovalveexp_info"; #puts gene fkpm into value
    $aovalveexp_hash{$aovalveexp_key} = $aovalveexp_value; #sets the key of the hash to the value
}
close AOVALVE_EXP;
print STDERR "Done\n\n";

# make a hash of brain expresssion
print STDERR "Retrieving gene expression values for brain tissue\n";
open(BRAIN_EXP, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/humanbrain_fkpm.anno") or die "Can't open BRAIN_EXP.anno: <$!>\n";
my %brainexp_hash;
for (<BRAIN_EXP>) {
    chomp;
    my($brainexp_position,$brainexp_info)=split(/@/); #splits the row into elements
    my $brainexp_key="$brainexp_position\t"; #puts chr and gene name and tab into one key
    my $brainexp_value="$brainexp_info"; #puts gene fkpm into value
    $brainexp_hash{$brainexp_key} = $brainexp_value; #sets the key of the hash to the value
}
close BRAIN_EXP;
print STDERR "Done\n\n";

# make a hash of disease genes in the omim database
print STDERR "Retrieving disease genes from the omim database\n";
open(OMIM, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/omim.anno") or die "Can't open OMIM.anno: <$!>\n";
my %omim_hash;
for (<OMIM>) {
    chomp;
    my($omim_chr,$omim_gene,$omim_disease)=split(/\t/); #splits the row into elements
    my $omim_key="$omim_chr\t$omim_gene\t"; #puts elements into one key
    my $omim_value="$omim_disease"; #puts text into value
    $omim_hash{$omim_key} = $omim_value; #sets the key of the hash to the value
}
close OMIM;
print STDERR "Done\n\n";

# make a hash of clinvar snps
print STDERR "Retrieving SNP variations from the ClinVar database\n";
open(CLINVAR, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/clinvar_snp.anno") or die "Can't open clinvar_snp.anno: <$!>\n";
my %clinvar_hash;
for (<CLINVAR>) {
    chomp;
    my($clinvar_chr,$clinvar_pos,$clinvar_id,$clinvar_prediction)=split(/\t/); #splits the row into elements
    my $clinvar_key="$clinvar_chr\t$clinvar_pos\t"; #puts elements into one key
    my $clinvar_value="$clinvar_id\t$clinvar_prediction"; #puts text into value
    $clinvar_hash{$clinvar_key} = $clinvar_value; #sets the key of the hash to the value
}
close CLINVAR;
print STDERR "Done\n\n";

# make a hash of hgmd snps
print STDERR "Retrieving SNP variations from the Human Gene Mutation database\n";
open(HGMD, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/hgmd_snps.anno") or die "Can't open hgmd_snps.anno: <$!>\n";
my %hgmd_hash;
for (<HGMD>) {
    chomp;
    my($hgmd_chr,$hgmd_pos,$hgmd_id)=split(/\t/); #splits the row into elements
    my $hgmd_key="$hgmd_chr\t$hgmd_pos\t"; #puts elements into one key
    my $hgmd_value="$hgmd_id"; #puts text into value
    $hgmd_hash{$hgmd_key} = $hgmd_value; #sets the key of the hash to the value
}
close HGMD;
print STDERR "Done\n\n";

# make a hash of transcript strand
print STDERR "Retrieving transcript strand information\n";
open(STRAND, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/strand.anno") or die "Can't open strand.anno: <$!>\n";
my %strand_hash;
for (<STRAND>) {
    chomp;
    my($strand_chr,$strand_transcript,$strand_strand)=split(/\t/); #splits the row into elements
    my $strand_key="$strand_chr\t$strand_transcript\t"; #puts elements into one key
    my $strand_value="$strand_strand"; #puts text into value
    $strand_hash{$strand_key} = $strand_value; #sets the key of the hash to the value
}
close STRAND;
print STDERR "Done\n\n";

# make a hash of RVIS
print STDERR "Retrieving residual variation intolerance scores\n";
open(RVIS, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/RVIS.anno") or die "Can't open RVIS.anno: <$!>\n";
my %rvis_hash;
for (<RVIS>) {
    chomp;
    my($rvis_gene,$rvis_score,$rvis_percent)=split(/\t/); #splits the row into elements
    my $rvis_key="$rvis_gene\t"; #puts elements into one key
    my $rvis_value="$rvis_score\t$rvis_percent"; #puts text into value
    $rvis_hash{$rvis_key} = $rvis_value; #sets the key of the hash to the value
}
close RVIS;
print STDERR "Done\n\n";

# make vcf meta and full genotype hash
my @metadata_vcf = glob("$input/$cohort/Variants/*.cohort.filtered.snps.vcf");
for (my $i = 0; $i < @metadata_vcf; $i++) {
    vcf_meta($metadata_vcf[$i]);
}

my %dp_hash;
open(VCFMETA, "<$input/VCFMETA_anno") or die "error reading $input/VCFMETA_anno: <$!>\n";
for (<VCFMETA>) {
    chomp;
    my($dp_chr,$dp_pos,$dp_ref,$AC,$AF,$AN,$QD,$MQ,$MQ0,$dp_ave,$dp_range)=split(/\t/); #splits the row into elements
    my $dp_key="$dp_chr\t$dp_pos\t$dp_ref"; #puts elements into one key
    my $dp_value="$QD\t$MQ\t$MQ0\t$dp_ave\t$dp_range"; #puts text into value
    $dp_hash{$dp_key} = $dp_value; #sets the key of the hash to the value
}
close VCFMETA;
unlink ("$input/VCFMETA_anno");

my %genotypemeta_hash;
open(GENOTYPEMETA, "<$input/VCFFULL_anno") or die "error reading $input/VCFFULL_anno: <$!>\n";
for (<GENOTYPEMETA>) {
    chomp;
    my($genotypemeta_position,$genotypemeta_meta)=split(/@/); #splits the row into elements
    my $genotypemeta_key="$genotypemeta_position"; #puts elements into one key
    my $genotypemeta_value="$genotypemeta_meta"; #puts text into value
    $genotypemeta_hash{$genotypemeta_key} = $genotypemeta_value; #sets the key of the hash to the value
}
close GENOTYPEMETA;
unlink ("$input/VCFFULL_anno");

# create a final annotation file
open(FINAL_ANNO, ">>$input/final_anno") or die "error creating final_anno file: <$!>\n";
my @bed = ();
my @vcf_meta =();

# get sample names and add to a the header line of sample_genotypes
open (SAMPLES, "<$input/$cohort/Variants/$cohort.samples.txt") or die "error reading $input/$cohort/Variants/$cohort.samples.txt: <$!>\n";
my @samplenames = <SAMPLES>;
chomp(@samplenames);
close SAMPLES;

open(SAMPLE_GENOTYPES, ">>$input/sample_genotypes") or die "error creating sample_genotypes file: <$!>\n";
print SAMPLE_GENOTYPES join ("\t", @samplenames), "\n";

print STDERR "Adding annotations to Seattle Seq file\n";
# loop through seattle seq file
open(SEATTLE,"$input/$cohort/Variants/$cohort.snp.annotated.txt"); #the seattle seq annotation file
for (<SEATTLE>) {
    chomp;
	if  (/^\#/) { # skip all lines beginning with #
        next;
}
my @seattleline=split("\t"); # get array from each line of seattleseq annotation

# get data for hash table lookups
my $alt;
my @genotype = split("/", $seattleline[5]); # split sample alleles that were observed in the samples (NB sometimes they are triallelic in dbsnp)
if ($genotype[0] =~ m/$seattleline[3]/) {$alt = $genotype[1]} # if first observed allele matches the ref, the second observed allele is alt
else  {$alt = $genotype[0]} # else alt is the first observed allele
my $position = "$seattleline[1]\t$seattleline[2]\t$seattleline[3]"; # chromosome, position, reference base for ESP6500_snp
my $full_position = "$seattleline[1]\t$seattleline[2]\t$seattleline[3]\t$alt"; # get chromosome, position, refbase and alt base for 1KG_2011, ExAC and cadd lookup
my $chr_gene_tab = "$seattleline[1]\t$seattleline[21]\t"; # get chromosome and gene name and tab character for gene, heart, brain, aortic valve expression, omim lookup
my $chr_pos_tab = "$seattleline[1]\t$seattleline[2]\t"; # get chromosome and gene name and tab character for clinvar and hgmd lookup
my $chr_transcript_tab = "$seattleline[1]\t$seattleline[7]\t"; # get chromosome and transcript and tab character for strand
my $gene_tab = "$seattleline[21]\t"; # get gene and tab for RVIS lookup
push(@bed, $seattleline[1],"\t",($seattleline[2]-1),"\t",$seattleline[2],"\n"); # make bed file for protein annotation and targetscan lookup

# format final file for mysql database
# replace GERP 'NA' with 0
$seattleline[17] =~ s/^NA$/0/g;
# replace Grantham NA with 0
$seattleline[15] =~ s/^NA$/0/g;

# replace Grantham , with . and limit to two elements
my @grantham = split(",", $seattleline[15]);
if (scalar(@grantham)>1){$seattleline[15] = join(".", @grantham[0..1]);}
# replace none with NULL
foreach(@seattleline) {s/^none$/NULL/g;}
# replace unknown with NULL
foreach(@seattleline) {s/^unknown$/NULL/g;}
# replace NA with NULL
foreach(@seattleline) {s/^NA$/NULL/g;}

# annotate intron-near-splice
if (($seattleline[8] =~ m/intron-near-splice/)&&($seattleline[17] >=5.5)) {$seattleline[8] =~ s/intron-near-splce/intron-near-splice-hc/g};
if (($seattleline[8] =~ m/intron-near-splice/)&&($seattleline[17] <5.5)&&($seattleline[17] >4.0)) {$seattleline[8] =~ s/intron-near-splce/intron-near-splice-c/g};

# then annotate using ESP6500 hash
if (exists $ESP6500_hash{$position}) {push(@seattleline, "$ESP6500_hash{$position}")}
else                       {push(@seattleline, "0\t0\t0")}
# then annotate using 1KG_2011 hash
if (exists $KG2011_hash{$full_position}) {push(@seattleline, "$KG2011_hash{$full_position}")}
else                       {push(@seattleline, "0\t0\t0\t0\t0\t0")}
# then annotate using gene_info hash
if (exists $geneinfo_hash{$chr_gene_tab}) {push(@seattleline, "$geneinfo_hash{$chr_gene_tab}")}
else                       {push(@seattleline, "NULL\tNULL")}
# then annotate using heart_exp hash
if (exists $heartexp_hash{$chr_gene_tab}) {push(@seattleline, "$heartexp_hash{$chr_gene_tab}")}
else                       {push(@seattleline, "NULL")}
# then annotate using aovalve_exp hash
if (exists $aovalveexp_hash{$chr_gene_tab}) {push(@seattleline, "$aovalveexp_hash{$chr_gene_tab}")}
else                       {push(@seattleline, "NULL")}
# then annotate using brain_exp hash
if (exists $brainexp_hash{$chr_gene_tab}) {push(@seattleline, "$brainexp_hash{$chr_gene_tab}")}
else                       {push(@seattleline, "NULL")}
# then annotate using omim hash
if (exists $omim_hash{$chr_gene_tab}) {push(@seattleline, "$omim_hash{$chr_gene_tab}")}
else                       {push(@seattleline, "NULL")}
# then annotate using clinvar hash
if (exists $clinvar_hash{$chr_pos_tab}) {push(@seattleline, "$clinvar_hash{$chr_pos_tab}")}
else                       {push(@seattleline, "NULL\tNULL")}
# then annotate using hgmd hash
if (exists $hgmd_hash{$chr_pos_tab}) {push(@seattleline, "$hgmd_hash{$chr_pos_tab}")}
else                       {push(@seattleline, "NULL")}
# then annotate using ExAC hash
if (exists $ExAC_hash{$full_position}) {push(@seattleline, "$ExAC_hash{$full_position}")}
else                       {push(@seattleline, "0")}
# then annotate using rvis hash
if (exists $rvis_hash{$gene_tab}) {push(@seattleline, "$rvis_hash{$gene_tab}")}
else                       {push(@seattleline, "NULL\tNULL")}
# then annotate cDNA numbering (missense, silent and nonsense only)
my $strand; # value of + or -

if ($seattleline[13] =~ m/NULL/) {
    push(@seattleline, "NULL");
} elsif (($seattleline[7] =~ m/_/)&&($seattleline[13] !~ m/NULL/)) {
    my @transcript = split (/\./, $seattleline[7]); # split NM_12345.6 at the full stop (hash table does not have values after the full stop)
    # see if the transcript number or gene name is in the strand hash
    my $new_transcript_value = "$seattleline[1]\t$transcript[0]\t"; # 1 NM_123
    my $new_gene_value = "$seattleline[1]\t$seattleline[21]\t"; # 1 MYBPC3
    if (exists $strand_hash{$new_transcript_value}) {
        $strand = "$strand_hash{$new_transcript_value}";
        if    ($strand =~ m/\+/) {
            push (@seattleline, "$seattleline[7]:c.$seattleline[13]$seattleline[3]>$alt")
        } elsif ($strand =~ m/\-/) {
            my $rev_ref = $seattleline[3]; $rev_ref =~ tr/ACGT/TGCA/;
            my $rev_alt = $alt; $rev_alt =~ tr/ACGT/TGCA/;
            push(@seattleline, "$seattleline[7]:c.$seattleline[13]$rev_ref>$rev_alt");
            }
        }
    # if not, see if the gene name is in the strand hash
    elsif (exists $strand_hash{$new_gene_value}) {
        $strand = "$strand_hash{$new_gene_value}";
        if ($strand =~ m/\+/) {
            push (@seattleline, "$seattleline[7]:c.$seattleline[13]$seattleline[3]>$alt");
        }
        elsif ($strand =~ m/\-/) {
            my $rev_ref = $seattleline[3]; $rev_ref =~ tr/ACGT/TGCA/;
            my $rev_alt = $alt; $rev_alt =~ tr/ACGT/TGCA/;
            push(@seattleline, "$seattleline[7]:c.$seattleline[13]$rev_ref>$rev_alt");
        }
    }

    else {
        push(@seattleline, "$seattleline[7]:c.$seattleline[13]");
    }
} else {
    push(@seattleline, "NULL");
}

# add strand as + or -
# see if the transcript or gene is in the strand hash
if  (/^\#/) { # skip all lines beginning with #
next;
}
elsif ($seattleline[7] =~ m/_/) {
    my @transcript =split (/\./, $seattleline[7]);
    my $new_transcript_value = "$seattleline[1]\t$transcript[0]\t";
    my $new_gene_value = "$seattleline[1]\t$seattleline[21]\t"; # 1 MYBPC3
    if (exists $strand_hash{$new_transcript_value}) {
        push(@seattleline, "$strand_hash{$new_transcript_value}");
    } elsif (exists $strand_hash{$new_gene_value}) {
        push(@seattleline, "$strand_hash{$new_gene_value}");
    } else {
        push(@seattleline, "NULL");
    }
}
else {
    push(@seattleline, "NULL");
}

# add alt allele
push(@seattleline, "$alt");
# add Alamut ID
my @alamut = ("Chr","$seattleline[1]","(GRCh37):g.","$seattleline[2]","$seattleline[3]",">","$alt");
push(@seattleline, (join("", @alamut)));

# add protein substitutions
if  (/^\#/) { # skip all lines beginning with #
next;
}
elsif ($seattleline[11] =~ m/NULL/) {push (@seattleline, "NULL")}

elsif ($seattleline[8] =~ m/synonymous/) {$seattleline[11] =~ tr/A-Z/a-z/;
    my $aa_ref = ucfirst($seattleline[11]);
    my @protein_residue = split ("/", $seattleline[12]); # 343/682
    my @synonymous = ("$aa_ref","$protein_residue[0]","$aa_ref"); # Ser343Ser
    push (@seattleline, (join("", @synonymous)))}


else {my @aa_substitution = split (",", $seattleline[11]); # SER,ARG
    $aa_substitution[0] =~ tr/A-Z/a-z/;
    my $aa_ref = ucfirst($aa_substitution[0]);
    $aa_substitution[1] =~ tr/A-Z/a-z/;
    my $aa_alt = ucfirst($aa_substitution[1]);
    my @protein_residue = split ("/", $seattleline[12]); # 343/682
    my @protein_substitution = ("$aa_ref","$protein_residue[0]","$aa_alt"); # Ser343Arg
    push (@seattleline, (join("", @protein_substitution)))}

# add transcript notation

# then annotate depth(meta) hash
if (exists $dp_hash{$position}) {push(@seattleline, "$dp_hash{$position}")}
else                       {push(@seattleline, "NULL\tNULL\tNULL\tNULL\tNULL")}
# then annotate vcf_meta_data array
if (exists $genotypemeta_hash{$position}) {push(@vcf_meta, "$genotypemeta_hash{$position}")}
else                       {push(@vcf_meta, "METAFAIL")}

# determine zygosity
my @zygosity = ('NULL', 'NULL', 'NULL', 'NULL'); # empty elements for the sanger sequenced, by, date, comments columns
my @genotypes = split(",", $seattleline[4]);
for (my $i = 0; $i <@genotypes; $i++) {
    if ($genotypes[$i] =~ m/$seattleline[3]/) {push (@zygosity, "0/0")}
    elsif ($genotypes[$i] =~ /^(A|C|G|T)$/) {push (@zygosity, "1/1")}
    elsif ($genotypes[$i] =~ /^(R|Y|S|W|K|M)$/) {push (@zygosity, "0/1")}
    elsif ($genotypes[$i] =~ m/^N$/) {push (@zygosity, "./.")}
}
my $catgenotypes = join (" ", @zygosity);
my $altcount = ($catgenotypes =~ tr/1/1/);
my $refcount = ($catgenotypes =~ tr/0/0/);
my $allelecount = ($refcount + $altcount);
my $allelefreq = sprintf("%.4f", ($altcount / ($altcount + $refcount)));
unshift (@zygosity, ("$altcount", "$allelefreq", "$allelecount"));

print FINAL_ANNO join ("\t", @seattleline), "\n";

print SAMPLE_GENOTYPES join ("\t", @zygosity), "\n";
} # end of seattleseq <>

open(VCF_META_FILE,">$input/vcf_meta_file");

my @vcf_meta_header=();
for (my $i = 0; $i < @samplenames; $i++) {
    push(@vcf_meta_header, "$samplenames[$i]_meta");
}
print VCF_META_FILE join ("\t", @vcf_meta_header), "\n";
print VCF_META_FILE join ("\n", @vcf_meta);

close SEATTLE;
close FINAL_ANNO;
close SAMPLE_GENOTYPES;
close VCF_META_FILE;
print STDERR "Done\n\n";

########################## annotate with protein domains using bedtools ###########################

print STDERR "Creating bed file of variants\n";

open(BED, ">$input/positions.bed") or die "error reading position.bed for reading";
print BED @bed;

print STDERR "Finding varints residing in protein domains\n";

system("intersectBed -a $input/positions.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/protein_domains.bed -wb > $input/domains.bed");
system("intersectBed -a $input/positions.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/protein_domains.bed -v -wb | cat $input/domains.bed - | sort -k1,1V -k2,2g > $input/domains.sorted.bed");

print STDERR "Done\n\n";
print STDERR "Retrieving protein domains. This may take a minute...\n";


open(DOMAINS,"<$input/domains.sorted.bed");
while (<DOMAINS>) {
    chomp;
# fill in NULL values where there was no match
# get array from each line of domains
    my @domain_line=split("\t");
    if ($domain_line[6]) {} # if array [6] is not empty, do nothing
    else {(push(@domain_line, "NULL", "NULL", "NULL", "NULL"))}

    open(FOLD, ">>$input/folds.txt") or die "error reading $input/folds.txt for reading";
    print FOLD $domain_line[6],"\n";

}
close DOMAINS;
close FOLD;
print STDERR "Done\n\n";

############################ annotate with targetscan domains using bedtools ############################
print STDERR "Finding varints residing in microRNA binding sites\n";

system("intersectBed -a $input/positions.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/targetscan.bed -wb > $input/targetscan.bed");
system("intersectBed -a $input/positions.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/targetscan.bed -v -wb | cat $input/targetscan.bed - | sort -k1,1V -k2,2g > $input/targetscan.sorted.bed");

print STDERR "Done\n\n";
print STDERR "Retrieving target scan sites. This may take a minute...\n";

open(TARGETSCAN,"<$input/targetscan.sorted.bed");
while (<TARGETSCAN>) {
    chomp;
# fill in NULL values where there was no match
# get array from each line of domains
    my @targetscan_line=split("\t");
    if ($targetscan_line[6]) {} # if array [6] is not empty, do nothing
    else {(push(@targetscan_line, "NULL", "NULL", "NULL", "NULL", "NULL"))}

    open(MIRNA, ">>$input/mirna_sites.txt") or die "error reading $input/mirna_sites.txt for reading";
    print MIRNA "$targetscan_line[6]\t$targetscan_line[7]\n";

}
close TARGETSCAN;
close MIRNA;
print STDERR "Done\n\n";

########################## create final snp annotation file ################################
print STDERR "Printing final annotation file\n";
system("echo \"Chr\tNucleotide\tRef\tAlt\tAlamut\tSample_alleles\tGene\tGene_name\tGene_description\tHeart_fkpm\tAortic_valve_fkpm\tBrain_fkpm\tTranscript_ID\tStrand\tcDNA_variation\tDistance_to_splice\tFunctionGVS\tFunctionDBSNP\tAmino_acid\tGrantham\tProtein_residue\tProtein_domain\trs_ID\tESP6500_EA\tESP6500_AA\tESP6500_alleles\t1KG_AF\t1KG_AMR\t1KG_EAS\t1KG_SAS\t1KG_AFR\t1KG_EUR\tExAC_count\t46way_Phast\tGERP\tOMIM\tKEGG\tRVIS\tRVIS_percent\tCADD_score\tPolyPhen\tClinVar_ID\tClinVar_prediction\tHGMD_ID\tClinical_association\tMiRNA_site\tMiRNA_score\tQD\tMQ\tMQ0\tDepth_ave\tDepth_range\tAlt_count\tAlt_freq\tAllele_count\tSanger_sequenced\tSequenced_by\tDate\tComments\" > $input/$cohort/Variants/$cohort.annotated_snps.text");

system("paste $input/final_anno $input/folds.txt $input/mirna_sites.txt| awk -v OFS=\"\t\" -F\"\t\" '{print\$2,\$3,\$4,\$62,\$63,\$6,\$22,\$48,\$49,\$50,\$51,\$52,\$8,\$61,\$60,\$31,\$9,\$10,\$64,\$16,\$13,\$70,\$11,\$40,\$41,\$36,\$42,\$43,\$44,\$45,\$46,\$47,\$57,\$17,\$18,\$53,\$33,\$58,\$59,\$19,\$15,\$54,\$55,\$56,\$30,\$71,\$72,\$65,\$66,\$67,\$68,\$69}' >> $input/$cohort/Variants/$cohort.annotated_snps.text");

system("paste $input/$cohort/Variants/$cohort.annotated_snps.text $input/sample_genotypes $input/vcf_meta_file > $input/$cohort/Variants/$cohort.annotated_snps.txt");

unlink ("$input/domains.bed");
unlink ("$input/domains.sorted.bed");
unlink ("$input/final_anno");
unlink ("$input/folds.txt");
unlink ("$input/positions.bed");
#unlink ("$input/$cohort/Variants/$cohort.snp.annotated.txt");
unlink ("$input/sample_genotypes");
unlink ("$input/$cohort/Variants/$cohort.annotated_snps.text");
unlink ("$input/vcf_meta_file");


#################################### Begin annotation of INDELS ################################################

print "\n\n*** Starting indel annotation pipeline ***\n";
clock ();

mkdir ("/home/groups/cardio/.do_not_delete_this_indel_folder") or die "Unable to create do_not_delete_this_indel_folder directory <$!>\n";
# make test.vcf with headers for submission to seattleseq annotation
open(TEST, ">/home/groups/cardio/.do_not_delete_this_indel_folder/indel_test.vcf") or die "Could not open indel_test.vcf: <$!>\n";
print TEST "\# autoFile testAuto.txt\n\# compressAuto true\n"; # need these two header lines for submission
# get headers and PASS sites from vcf file
system("grep '^#\\|PASS' < $input/$cohort/Variants/$cohort.cohort.filtered.indels.vcf >>/home/groups/cardio/.do_not_delete_this_indel_folder/indel_test.vcf");
# load test.vcf to seattleseq
chdir('/home/groups/cardio/Applications/seattle_indel_java'); # need to run from this directory [bloody java!]
system("java SubmitSeattleSeqAnnotationAutoJob");
system("gunzip /home/groups/cardio/.do_not_delete_this_indel_folder/indel_test.vcf.txt.gz");

# now have an file called indel_text.vcf.txt and need to get it back to sample name
my $current_indel_file_name = "/home/groups/cardio/.do_not_delete_this_indel_folder/indel_test.vcf.txt";
my $new_indel_file_extension = "$input/$cohort/Variants/$cohort.indels.annotated.txt";
move($current_indel_file_name, $new_indel_file_extension) or die "Unable to move $current_indel_file_name to $new_indel_file_extension: $!\n";
close TEST;
# cleanup
unlink("/home/groups/cardio/.do_not_delete_this_indel_folder/indel_test.vcf");
rmdir "/home/groups/cardio/.do_not_delete_this_indel_folder" or die "Unable to delete do_not_delete_this_indel_folder directory: <$!>\n";

###################################
#   Begin secondary annotations   #
###################################

print STDERR "\n\n***Starting secondary annotation of indel variants ***\n";
clock ();

print STDERR "Collecting ESP6500 indel variation frequencies\n";
# make a hash of ESP6500 variants
open(ESP6500_indel, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/ESP6500_all_indels.anno") or die "Can't open ESP6500.indel.anno: <$!>\n";
my %ESP6500_indel_hash;
for (<ESP6500_indel>) {
    chomp;
    my($esp_indel_chr,$esp_indel_position,$esp_indel_ref,$esp_indel_alt,$esp_indel_EUfreq,$esp_indel_AAfreq,$esp_indel_ALLfreq)=split(/\t/); #splits the row into elements
    my $ESP6500_indel_key="$esp_indel_chr\t$esp_indel_position\t$esp_indel_ref\t$esp_indel_alt"; #puts first 4 columns of anno file into one key
    my $ESP6500_indel_value="$esp_indel_EUfreq\t$esp_indel_AAfreq\t$esp_indel_ALLfreq"; #puts gene name and description into value
    $ESP6500_indel_hash{$ESP6500_indel_key} = $ESP6500_indel_value; #sets the key of the hash to the value
}
close ESP6500_indel;
print STDERR "Done\n\n";

# make a hash of 1KG variants
print STDERR "Collecting 1000 genomes phase1.b37 indel variation frequencies\n";
open(KG_2011_indel, "/home/shared/NGS/human/richardb/Annotations/ALL.autosomes.phase3.indel.anno") or die "Can't open 1KG_2011_indel.anno: <$!>\n";
my %KG2011_indel_hash;
for (<KG_2011_indel>) {
    chomp;
    my($kgindel_chr,$kgindel_position,$kgindel_ref,$kgindel_alt,$kgindel_AFfreq,$kgindel_AMR,$kgindel_EAS,$kgindel_SAS,$kgindel_AFR,$kgindel_EUR)=split(/\t/); #splits the row into elements
    my $kgindel_key="$kgindel_chr\t$kgindel_position\t$kgindel_ref\t$kgindel_alt"; #puts chr position and ref into one key
    my $kgindel_value="$kgindel_AFfreq\t$kgindel_AMR\t$kgindel_EAS\t$kgindel_SAS\t$kgindel_AFR\t$kgindel_EUR"; #puts freq into value
    $KG2011_indel_hash{$kgindel_key} = $kgindel_value; #sets the key of the hash to the value
}
close KG_2011_indel;
print STDERR "Done\n\n";

# make a hash of clinvar indels
print STDERR "Retrieving INDEL variations from the ClinVar database\n";
open(CLINVAR_INDEL, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/clinvar_indels.anno") or die "Can't open clinvar_indels.anno: <$!>\n";
my %clinvar_indel_hash;
for (<CLINVAR_INDEL>) {
    chomp;
    my($clinvar_indel_chr,$clinvar_indel_pos,$clinvar_indel_id,$clinvar_indel_prediction)=split(/\t/); #splits the row into elements
    my $clinvar_indel_key="$clinvar_indel_chr\t$clinvar_indel_pos\t"; #puts elements into one key
    my $clinvar_indel_value="$clinvar_indel_id\t$clinvar_indel_prediction"; #puts text into value
    $clinvar_indel_hash{$clinvar_indel_key} = $clinvar_indel_value; #sets the key of the hash to the value
}
close CLINVAR_INDEL;
print STDERR "Done\n\n";

# make a hash of hgmd indels
print STDERR "Retrieving INDEL variations from the Human Gene Mutation database\n";
open(HGMD_INDEL, "/home/groups/cardio/References/Annotations/Exome_seq_pipeline/hgmd_indels.anno") or die "Can't open hgmd_indels.anno: <$!>\n";
my %hgmd_indel_hash;
for (<HGMD_INDEL>) {
    chomp;
    my($hgmd_indel_chr,$hgmd_indel_pos,$hgmd_indel_id)=split(/\t/); #splits the row into elements
    my $hgmd_indel_key="$hgmd_indel_chr\t$hgmd_indel_pos\t"; #puts elements into one key
    my $hgmd_indel_value="$hgmd_indel_id"; #puts text into value
    $hgmd_indel_hash{$hgmd_indel_key} = $hgmd_indel_value; #sets the key of the hash to the value
}
close HGMD_INDEL;
print STDERR "Done\n\n";

# make a hash of ExAC indel variants
print STDERR "Retrieving ExAC indel variation counts\n";
open(INDELExAC, "/home/shared/NGS/human/richardb/Annotations/ExAC.r0.1.sites.indel.anno") or die "Can't open ExAC.indel.anno: <$!>\n";
my %indelExAC_hash;
for (<INDELExAC>) {
    chomp;
    my($indelexac_chr,$indelexac_pos,$indelexac_ref,$indelexac_alt,$indelexac_count)=split(/\t/); #splits the row into elements
    my $indelExAC_key="$indelexac_chr\t$indelexac_pos\t$indelexac_ref\t$indelexac_alt"; #puts chrom, position, ref and alt of anno file into one key
    my $indelExAC_value="$indelexac_count"; #puts allele count into value
    $indelExAC_hash{$indelExAC_key} = $indelExAC_value; #sets the key of the hash to the value
}
close INDELExAC;
print STDERR "Done\n\n";

# make indel vcf meta and full genotype hashes
my @indelmetadata_vcf = glob("$input/$cohort/Variants/*.cohort.filtered.indels.vcf");
for (my $i = 0; $i < @indelmetadata_vcf; $i++) {
    vcf_meta($indelmetadata_vcf[$i]);
}

my %indelmeta_hash;
open(INDELVCFMETA, "<$input/VCFMETA_anno") or die "error reading $input/VCFMETA_anno: <$!>\n";
for (<INDELVCFMETA>) {
    chomp;
    my($ind_chr,$ind_pos,$ind_ref,$indAC,$indAF,$indAN,$indQD,$indMQ,$indMQ0,$ind_ave,$ind_range)=split(/\t/); #splits the row into elements
    my $ind_key="$ind_chr\t$ind_pos\t$ind_ref"; #puts elements into one key
    my $ind_value="$indQD\t$indMQ\t$indMQ0\t$ind_ave\t$ind_range"; #puts text into value
    $indelmeta_hash{$ind_key} = $ind_value; #sets the key of the hash to the value
}
close INDELVCFMETA;
unlink ("$input/VCFMETA_anno");

my %indelgenotypemeta_hash;
open(INDELGENOTYPEMETA, "<$input/VCFFULL_anno") or die "error reading $input/VCFFULL_anno: <$!>\n";
for (<INDELGENOTYPEMETA>) {
    chomp;
    my($indgenotypemeta_position,$indgenotypemeta_meta)=split(/@/); #splits the row into elements
    my $indgenotypemeta_key="$indgenotypemeta_position"; #puts elements into one key
    my $indgenotypemeta_value="$indgenotypemeta_meta"; #puts text into value
    $indelgenotypemeta_hash{$indgenotypemeta_key} = $indgenotypemeta_value; #sets the key of the hash to the value
}
close INDELGENOTYPEMETA;
unlink ("$input/VCFFULL_anno");

# create a final annotation file
open(FINAL_INDEL_ANNO, ">>$input/final_indel_anno") or die "error creating final_indel_anno file: <$!>\n";
my @indelvcf_meta =();

# get sample names and add to a the header line of sample_genotypes
open(SAMPLE_INDEL_GENOTYPES, ">>$input/sample_indel_genotypes") or die "error creating sample_indel_genotypes file: <$!>\n";
print SAMPLE_INDEL_GENOTYPES join ("\t", @samplenames), "\n";

my @bed_indel = ();

# loop through seattle seq file
open(SEATTLE_INDEL,"$input/$cohort/Variants/$cohort.indels.annotated.txt"); #the seattle seq annotation file
for (<SEATTLE_INDEL>) {
    chomp;
    # skip all lines beginning with #
	if  (/^\#/) {
        next;
}
# get array from each line of seattleseq annotation
my @seattleline_indel=split("\t");
# get chromosome, position, refbase for DP lookup
my $position_indel = "$seattleline_indel[1]\t$seattleline_indel[2]\t$seattleline_indel[3]";
# get chromosome, position, refbase and alt base for ESP6500 and 1KG_2011 lookup
my @alt_indels = split(",", $seattleline_indel[5]);
my $full_position_indel = "$seattleline_indel[1]\t$seattleline_indel[2]\t$seattleline_indel[3]\t$alt_indels[0]";
# get chromosome and gene name for gene info lookup
my $chr_gene_indel = "$seattleline_indel[1]\t$seattleline_indel[20]";
# get chromosome and gene and tab for omim info lookup
my $chr_gene_tab_indel = "$seattleline_indel[1]\t$seattleline_indel[20]\t";
# get chromosome position tab
my $chr_pos_tab_indel = "$seattleline_indel[1]\t$seattleline_indel[2]\t";
# get gene tab for rvis lookup
my $gene_tab_indel = "$seattleline_indel[20]\t";
push(@bed_indel, $seattleline_indel[1],"\t",($seattleline_indel[2]-1),"\t",$seattleline_indel[2],"\n"); # make bed file for protein annotation and targetscan lookup

# replace 46Way NA with 0
$seattleline_indel[16] =~ s/^NA$/0/g;
# replace none with NULL
foreach(@seattleline_indel) {s/^none$/NULL/g;}
# replace unknown with NULL
foreach(@seattleline_indel) {s/^unknown$/NULL/g;}

# then annotate using ESP6500 hash
if (exists $ESP6500_indel_hash{$full_position_indel}) {push(@seattleline_indel, "$ESP6500_indel_hash{$full_position_indel}")}
else                       {push(@seattleline_indel, "0\t0\t0")}
# then annotate using 1KG_2011 hash
if (exists $KG2011_indel_hash{$full_position_indel}) {push(@seattleline_indel, "$KG2011_indel_hash{$full_position_indel}")}
else                       {push(@seattleline_indel, "0\t0\t0\t0\t0")}
# then annotate using gene_info hash
if (exists $geneinfo_hash{$chr_gene_tab_indel}) {push(@seattleline_indel, "$geneinfo_hash{$chr_gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL\tNULL")}
# then annotate using heart_exp hash
if (exists $heartexp_hash{$chr_gene_tab_indel}) {push(@seattleline_indel, "$heartexp_hash{$chr_gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL")}
# then annotate using aovalve_exp hash
if (exists $aovalveexp_hash{$chr_gene_tab_indel}) {push(@seattleline_indel, "$aovalveexp_hash{$chr_gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL")}
# then annotate using brain_exp hash
if (exists $brainexp_hash{$chr_gene_tab_indel}) {push(@seattleline_indel, "$brainexp_hash{$chr_gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL")}
# then annotate using omim hash
if (exists $omim_hash{$chr_gene_tab_indel}) {push(@seattleline_indel, "$omim_hash{$chr_gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL")}
# then annotate using clinvar hash
if (exists $clinvar_indel_hash{$chr_pos_tab_indel}) {push(@seattleline_indel, "$clinvar_indel_hash{$chr_pos_tab_indel}")}
else                       {push(@seattleline_indel, "NULL\tNULL")}
# then annotate using hgmd hash
if (exists $hgmd_indel_hash{$chr_pos_tab_indel}) {push(@seattleline_indel, "$hgmd_indel_hash{$chr_pos_tab_indel}")}
else                       {push(@seattleline_indel, "NULL")}
# then annotate using rvis hash
if (exists $rvis_hash{$gene_tab_indel}) {push(@seattleline_indel, "$rvis_hash{$gene_tab_indel}")}
else                       {push(@seattleline_indel, "NULL\tNULL")}
# then annotate depth hash
if (exists $indelmeta_hash{$position_indel}) {push(@seattleline_indel, "$indelmeta_hash{$position_indel}")}
else                       {push(@seattleline_indel, "NULL\tNULL\tNULL\tNULL\tNULL")}

# then annotate indelvcf_meta array
if (exists $indelgenotypemeta_hash{$position_indel}) {push(@indelvcf_meta, "$indelgenotypemeta_hash{$position_indel}")}
else                       {push(@indelvcf_meta, "METAFAIL\n")}

# determine zygosity
my @indel_zygosity = ('NULL', 'NULL', 'NULL', 'NULL');
my @indel_genotypes = split(",", $seattleline_indel[4]);
for (my $i = 0; $i <@indel_genotypes; $i++) {
	
	my @indel_alleles = split("/", $indel_genotypes[$i]);
	
    if (($indel_alleles[0]  eq $seattleline_indel[3])&&($indel_alleles[1] eq $seattleline_indel[3])) {push (@indel_zygosity, "0/0")}
    elsif (($indel_alleles[0] eq $seattleline_indel[5])&&($indel_alleles[1] eq $seattleline_indel[5])) {push (@indel_zygosity, "1/1")}
    elsif (($indel_alleles[0] eq $seattleline_indel[3])&&($indel_alleles[1] eq $seattleline_indel[5])) {push (@indel_zygosity, "0/1")}
    elsif (($indel_alleles[0] eq $seattleline_indel[5])&&($indel_alleles[1] eq $seattleline_indel[3])) {push (@indel_zygosity, "0/1")}
    elsif (($indel_alleles[0] =~ m/N/)&&($indel_alleles[1] =~ m/N/)) {push (@indel_zygosity, "./.")}
    elsif ($indel_alleles[0] !~ m/N/) {push (@indel_zygosity, "?/?")}
}
my $indel_catgenotypes = join (" ", @indel_zygosity);
my $indel_altcount = ($indel_catgenotypes =~ tr/1/1/);
my $indel_refcount = ($indel_catgenotypes =~ tr/0/0/);
my $indel_allelecount = ($indel_refcount + $indel_altcount);
if ($indel_refcount==0) {unshift (@indel_zygosity, "NULL", "NULL", "NULL");}
else                    {my $indel_allelefreq = sprintf("%.4f", ($indel_altcount / ($indel_altcount + $indel_refcount)));
    unshift (@indel_zygosity, ("$indel_altcount", "$indel_allelefreq", "$indel_allelecount"));}

print FINAL_INDEL_ANNO join ("\t", @seattleline_indel), "\n";
print SAMPLE_INDEL_GENOTYPES join ("\t", @indel_zygosity), "\n";
} # end of seattleseq <>

open(INDELVCF_META_FILE,">$input/indelvcf_meta_file");

my @meta_indelheader=();
for (my $i = 0; $i < @samplenames; $i++) {
    push(@meta_indelheader, "$samplenames[$i]_meta");
}

print INDELVCF_META_FILE join ("\t", @meta_indelheader), "\n";

print INDELVCF_META_FILE join("\n", @indelvcf_meta);

close SEATTLE_INDEL;
close FINAL_INDEL_ANNO;
close SAMPLE_INDEL_GENOTYPES;


########################## annotate with protein domains using bedtools ###########################

print STDERR "Creating bed file of indel variants\n";

open(BED_INDEL, ">$input/positions_indel.bed") or die "error reading positions_indel.bed for reading";

print BED_INDEL @bed_indel;

print STDERR "Finding indel varints residing in protein domains\n";

system("intersectBed -a $input/positions_indel.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/protein_domains.bed -wb > $input/domains_indel.bed");
system("intersectBed -a $input/positions_indel.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/protein_domains.bed -v -wb | cat $input/domains_indel.bed - | sort -k1,1V -k2,2g > $input/domains_indel.sorted.bed");

print STDERR "Done\n\n";
print STDERR "Retrieving protein domains. This may take a minute...\n";


open(DOMAINS_INDEL,"<$input/domains_indel.sorted.bed");
while (<DOMAINS_INDEL>) {
    chomp;
# fill in NULL values where there was no match
# get array from each line of domains
    my @domain_indel_line=split("\t");
    if ($domain_indel_line[6]) {} # if array [6] is not empty, do nothing
    else {(push(@domain_indel_line, "NULL", "NULL", "NULL", "NULL"))}

    open(FOLD_INDEL, ">>$input/folds_indel.txt") or die "error reading $input/folds_indel.txt for reading";
    print FOLD_INDEL $domain_indel_line[6],"\n";

}
close DOMAINS_INDEL;
close FOLD_INDEL;
print STDERR "Done\n\n";

############################ annotate with targetscan domains using bedtools ############################
print STDERR "Finding indel varints residing in microRNA binding sites\n";

system("intersectBed -a $input/positions_indel.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/targetscan.bed -wb > $input/targetscan_indel.bed");
system("intersectBed -a $input/positions_indel.bed -b /home/groups/cardio/References/Annotations/Exome_seq_pipeline/targetscan.bed -v -wb | cat $input/targetscan_indel.bed - | sort -k1,1V -k2,2g > $input/targetscan_indel.sorted.bed");

print STDERR "Done\n\n";
print STDERR "Retrieving target scan sites. This may take a minute...\n";

open(TARGETSCAN_INDEL,"<$input/targetscan_indel.sorted.bed");
while (<TARGETSCAN_INDEL>) {
    chomp;
# fill in NULL values where there was no match
# get array from each line of domains
    my @targetscan_indel_line=split("\t");
    if ($targetscan_indel_line[6]) {} # if array [6] is not empty, do nothing
    else {(push(@targetscan_indel_line, "NULL", "NULL", "NULL", "NULL", "NULL"))}

    open(MIRNA_INDEL, ">>$input/mirna_indel_sites.txt") or die "error reading $input/mirna_indel_sites.txt for reading";
    print MIRNA_INDEL "$targetscan_indel_line[6]\t$targetscan_indel_line[7]\n";

}
close TARGETSCAN_INDEL;
close MIRNA_INDEL;
print STDERR "Done\n\n";

########################## create final indel annotation file ################################

system("echo \"Chr\tNucleotide\tRef\tAlt\tGene\tGene_name\tGene_description\tHeart_fkpm\tAortic_valve_fkpm\tBrain_fkpm\tTranscript_ID\tFunctionGVS\tFunctionDBSNP\trs_ID\tESP6500_EA\tESP6500_AA\tESP_All\t1KG_AF\t1KG_AMR\t1KG_EAS\t1KG_SAS\t\t1KG_AFR\t1KG_EUR\tExAC_count\t46way_Phast\tOMIM\tRVIS\tRVIS_percent\tClinvar_ID\tClinvar_prediction\tHGMD_ID\tKEGG_pathway\tProtein_domain\tMiRNA_site\tMiRNA_score\tQD\tMQ\tMQ0\tDP_ave\tDP_range\tAlt_count\tAlt_freq\tAllele_count\tSanger_sequenced\tSequenced_by\tDate\tComments\" > $input/$cohort/Variants/$cohort.annotated_indels.text");

system("paste $input/final_indel_anno $input/folds_indel.txt $input/mirna_indel_sites.txt | awk -v OFS=\"\t\" -F\"\t\" '{print\$2,\$3,\$4,\$6,\$21,\$47,\$48,\$49,\$50,\$51,\$8,\$9,\$10,\$11,\$38,\$39,\$40,\$41,\$42,\$43,\$44,\$45,\$46,\$56,\$17,\$52,\$57,\$58,\$53,\$54,\$55,\$32,\$64,\$65,\$66,\$59,\$60,\$61,\$62,\$63}' >> $input/$cohort/Variants/$cohort.annotated_indels.text");


system("paste $input/$cohort/Variants/$cohort.annotated_indels.text $input/sample_indel_genotypes $input/indelvcf_meta_file > $input/$cohort/Variants/$cohort.annotated.indels.txt");


print STDERR "Completed Exome analysis\n";
print STDERR "Indexed bam and vcf files are in $input/$cohort\n";
print STDERR "Final annotated snps are in $input/$cohort/Variants/$cohort.annotated_snps.txt\n";
print STDERR "Final annotated indels are in $input/$cohort/Variants/$cohort.indels.txt\n";
clock();


unlink ("$input/final_indel_anno");
unlink ("$input/sample_indel_genotypes");
unlink ("$input/$cohort/Variants/$cohort.annotated_indels.text");
#unlink ("$input/$cohort/Variants/$cohort.indels.annotated.txt");
unlink ("$input/indelvcf_meta_file");
unlink ("$input/domains_indel.bed");
unlink ("$input/domains_indel.sorted.bed");
unlink ("$input/folds_indel.txt");
unlink ("$input/targetscan_indel.bed");
unlink ("$input/targetscan_indel.sorted.bed");
unlink ("$input/mirna_indel_sites.txt");
unlink ("$input/positions_indel.bed");

################## SUBROUTINES ###################

sub clock {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	print "*** $theTime ***\n\n";
}

sub usage {
	print "\n\n\n\t---------EXOME SEQUENCING ANALYSIS PIPELINE------------------\n\n";
	print "\tUsage:\tExome_cohort_v1.pl [-fastq -cohort]\n";
	print "\t-fastq - is required. The path to the directory containing the fastq files. NB: do not add trailing / \n";
	print "\t-cohort - is required. Define a name for the final VCF file / \n";
	print "\tUse this script to align and annotate a cohort of exome sequencing data. Requires paired end fastq.gz compressed reads\n";
	print "\tReads must be stored in your /path/2/Rawdata folder and labelled:\n\n";
	print "\t\tsamplename1.f.fastq.gz\n";
	print "\t\tsamplename1.r.fastq.gz\n";
	print "\t\tsamplename2.f.fastq.gz\n";
	print "\t\tsamplename2.r.fastq.gz\n";
	print "\t\tetc..  \n\n";
	print "\tReplace samplename with a unique ID (e.g. blood code)\n";
	print "\tUses the Illumina TruSeq Exome intervals and SeattleSeq Annotation website for annotation of snps and indels\n\n";
	print "\tContact r.bagnall\@centenary.org.au\n";
	print "\t-------------------------------------------------------------\n\n";
	exit;
}

sub bwa_mem {
	my ($current_read_pair) = shift(@_);
	my @read_pair = split(" ", $current_read_pair);
	my $current_samplename = substr $read_pair[0], ($offset + 9), -11; # i.e. get IO2 from $input/Rawdata/IO2.f.fastq.gz
	my $RG ='@RG'; #need this to be able to print out @RG

    print "\n\n*** Aligning $current_samplename reads ***\n";
	clock();
    
        my $bwa_mem = system("bwa mem -PM -t 6 -R '$RG\tID:$current_samplename\tSM:$current_samplename\tPL:ILLUMINA' $path2ref $current_read_pair | samtools view -Su - | novosort - --threads 6 --removeDuplicates --ram 10G --output $path2bam/$current_samplename.badheader.bam");
    
        my $reheader = system("samtools view -H $path2bam/$current_samplename.badheader.bam | sed -e 's/CL:bwa.*//' | samtools reheader - $path2bam/$current_samplename.badheader.bam > $path2bam/$current_samplename.sorted.bam");
}

sub sortbam {
	my ($current_bam) = shift(@_);
	my $current_name = substr $current_bam, ($offset + 10), -13; # eg $input/tempBAMs/IO2.unsorted.bam or $input/tempBAMs/IO2.dedupped.bam
	print "\n\n*** Sorting $current_name bam file ***\n";
	clock();
	my $sortbam = system("samtools sort -o $path2bam/$current_name.sorted -T $path2bam $current_bam");
	unlink($current_bam) or die "Could not delete $current_bam\n"; # delete the unsorted bamfile
}

sub indexbam {
	my ($current_bam) = shift(@_);	
	my $current_name = substr $current_bam, ($offset + 10), -4; # eg $input/tempBAMs/IO2.sorted.bam
	print "\n\n*** Indexing $current_name file ***\n";
	clock();
	my $indexbam = system("samtools index $current_bam");
}

sub realigner_target_creator {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.sorted.bam
    my $current_name = substr $current_bam, ($offset + 10), -11;
	print "\n\n*** Creating $current_name target indel sites ***\n";
	clock();
	my $realigner_target_creator = system("java -jar $path2gatk -T RealignerTargetCreator -R $path2ref -o $path2bam/$current_name.indel.intervals -I $current_bam -L $path2mergedexome50 -known $path2indel1kg -known $path2indelmills");
}
	
sub indel_realigner {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.sorted.bam
    my $current_name = substr $current_bam, ($offset + 10), -11;
	print "\n\n*** Realigning $current_name bam file ***\n";
	clock();
	my $local_realignment = system("java -jar $path2gatk -T IndelRealigner -R $path2ref -o $path2bam/$current_name.realigned.bam -I $current_bam -L $path2mergedexome50 -known $path2indel1kg -known $path2indelmills -targetIntervals $path2bam/$current_name.indel.intervals");
}

sub base_recalibration {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.realigned.bam
	my $current_name = substr $current_bam, ($offset + 10), -14;
	print "\n\n*** Performing pre base quality score recalibration for $current_name ***\n";
	clock();
	my $bqsr_pre = system("java -jar $path2gatk -T BaseRecalibrator -R $path2ref -I $current_bam -knownSites $path2SNP -knownSites $path2indel1kg -knownSites $path2indelmills -L $path2mergedexome50 -o $path2bam/$current_name.prerecal_data.grp");
	print "\n\n*** Performing post base quality score recalibration for $current_name ***\n";
	clock();
    my $bqsr_post = system("java -jar $path2gatk -T BaseRecalibrator -R $path2ref -I $current_bam -knownSites $path2SNP -knownSites $path2indel1kg -knownSites $path2indelmills -L $path2mergedexome50 -BQSR $path2bam/$current_name.prerecal_data.grp -o $path2bam/$current_name.postrecal_data.grp");
    print "\n\n*** Printing recalibration plots for $current_name ***\n";
	clock();
    my $plots = system("java -jar $path2gatk -T AnalyzeCovariates -R $path2ref -before $path2bam/$current_name.prerecal_data.grp -after $path2bam/$current_name.postrecal_data.grp -plots $path2bam/$current_name.postrecal.pdf");
	print "\n\n*** Printing recalibrated reads for $current_name ***\n";
	clock();
	my $bqsr = system("java -jar $path2gatk -T PrintReads -R $path2ref -I $current_bam -BQSR $path2bam/$current_name.prerecal_data.grp -L $path2mergedexome50 -o $path2bam/$current_name.recalibrated.bam");

	unlink ($current_bam); # cleanup: delete bam files
	unlink ("$input/tempBAMs/$current_name.realigned.bai"); # cleanup: delete bam index files
}

sub reduce_reads {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.recalibrated.bam
	my $current_name = substr $current_bam, ($offset + 10), -17;
	print "\n\n*** Reducing reads for $current_name ***\n";
	clock();
	my $reduce_reads = system("java -jar $path2gatk -T ReduceReads -R $path2ref -I $current_bam -L $path2mergedexome50 -o $path2bam/$current_name.reduced.bam");
}

sub coverage {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.recalibrated.bam
	my $current_name = substr $current_bam, ($offset + 10), -17;
	print "\n\n*** Calculating coverage for $current_name ***\n";
	clock();
	my $coverage = system("samtools view -b $current_bam | coverageBed -abam stdin -b $path2targetregions -d | awk -F\"\t\" '{print\$5}' > $input/tempBAMs/$current_name.coverage");
    my $coverage1 = system("Rscript /home/groups/cardio/Applications/Rscripts/exome_coverage.R $input/tempBAMs/$current_name.coverage");
    unlink ("$input/tempBAMs/$current_name.coverage");
    print "\n\n*** Finished $current_name ***\n\n";
}

sub unified_genotyper {
	my ($current_bam_list) = shift(@_);
	# get names of samples
	my @recalibrated_bams = glob("$path2bam/*recalibrated.bam");
	my @filenames = ();
		# loop through array of bamfiles and extract filename 
		for (my $i = 0; $i < @recalibrated_bams; $i++) {
		my $filename = substr $recalibrated_bams[$i], ($offset + 10), -17;
		push(@filenames, "$filename\n");
		}
	open(SAMPLE_NAMES, ">$path2vcf/$cohort.samples.txt") or die "error creating $cohort.samples.txt file: <$!>\n";
	print SAMPLE_NAMES @filenames;
	print "\n\n*** Starting unified genotyper ***\n";
 	clock();
	my $unified_genotyper = system("java -jar $path2gatk -T UnifiedGenotyper -nt 8 -R $path2ref -I $current_bam_list -o $path2vcf/$cohort.cohort.raw.vcf -L $path2mergedexome10 -G StandardAnnotation -stand_emit_conf 10.0 -stand_call_conf 20.0 -dcov 200 -l INFO -rf BadCigar -glm BOTH")
}	

sub haplotype_caller {
	my ($current_bam_list) = shift(@_);
	# get names of samples
	my @recalibrated_bams = glob("$path2bam/*recalibrated.bam");
	my @filenames = ();
    # loop through array of bamfiles and extract filename
    for (my $i = 0; $i < @recalibrated_bams; $i++) {
		my $filename = substr $recalibrated_bams[$i], ($offset + 10), -17;
		push(@filenames, "$filename\n");
    }
	open(SAMPLE_NAMES, ">$input/$cohort.samples.txt") or die "error creating $cohort.samples.txt file: <$!>\n";
	print SAMPLE_NAMES @filenames;
	print "\n\n*** Starting haplotype caller ***\n";
 	clock();
	my $haplotype_caller = system("java -jar $path2gatk -T HaplotypeCaller -R $path2ref -I $current_bam_list --dbsnp $path2SNP -o $path2vcf/$cohort.cohort.raw.vcf -L $path2mergedexome10 -stand_emit_conf 10.0 -stand_call_conf 20.0 -rf BadCigar")
}

sub vqsr {
	my ($current_vcf) = shift(@_);
	print "\n\n*** Running VQSR for SNPs in $cohort.raw.vcf file ***\n";
	clock();
	my $vqsr_snps = system("java -jar $path2gatk -T VariantRecalibrator -R $path2ref -L $path2mergedexome10 -input $current_vcf -percentBad 0.01 -minNumBad 1000 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $path2hapmap -resource:omni,known=false,training=true,truth=true,prior=12.0 $path2omni -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $path2SNP -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -mode SNP -recalFile $path2vcf/$cohort.SNP.recal -tranchesFile $path2vcf/$cohort.SNP.tranches -rscriptFile $path2vcf/$cohort.SNP.plots.R");
	
	print "\n\n*** Running VQSR known INDELS in $cohort.raw.vcf file ***\n";
	clock();

	my $vqsr_indels = system("java -jar $path2gatk -T VariantRecalibrator -R $path2ref -L $path2mergedexome10 -input $current_vcf --maxGaussians 4 -percentBad 0.01 -minNumBad 1000 -resource:mills,known=false,training=true,truth=true,prior=12.0 $path2indelmills -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $path2SNP -an InbreedingCoeff -an ReadPosRankSum -an FS -an MQRankSum -mode INDEL -recalFile $path2vcf/$cohort.INDEL.recal -tranchesFile $path2vcf/$cohort.INDEL.tranches -rscriptFile $path2vcf/$cohort.INDEL.plots.R");
	
	print "\n\n*** Applying recalibration to SNPs in $cohort.raw.vcf file ***\n";
	clock();
	my $apply_recal_snps = system("java -jar $path2gatk -T ApplyRecalibration -R $path2ref -L $path2mergedexome10 -input $current_vcf -tranchesFile $path2vcf/$cohort.SNP.tranches -recalFile $path2vcf/$cohort.SNP.recal -o $path2vcf/$cohort.recalibrated.SNP.filtered.vcf --ts_filter_level 99.0 -mode SNP");
	
	print "\n\n*** Applying recalibration to INDELs in $cohort.recalibrated.SNP.filtered.vcf file ***\n";
	clock();
	my $apply_recal_indels = system("java -jar $path2gatk -T ApplyRecalibration -R $path2ref -L $path2mergedexome10 -input $path2vcf/$cohort.recalibrated.SNP.filtered.vcf -tranchesFile $path2vcf/$cohort.INDEL.tranches -recalFile $path2vcf/$cohort.INDEL.recal -o $path2vcf/$cohort.recalibrated.SNP.INDEL.filtered.vcf --ts_filter_level 99.0 -mode INDEL");
}

sub average{
    my($data) = @_;
    if (not @$data) {
        die("Empty array\n");
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}

sub vcf_meta{
    print STDERR "Retrieving meta data and full genotypes from VCF file\n";
    my (@current_vcf) = shift(@_);
    my @meta_anno=();
    my @full_genotypes=();
    
    open(CURRENTVCF, "$current_vcf[0]");
    for (<CURRENTVCF>){
        chomp;
        if  (/^\#/) {next;}
            
            my @vcf_complete=split("\t"); # get array from each line of vcf
        push (@meta_anno, "$vcf_complete[0]\t$vcf_complete[1]\t$vcf_complete[3]\t");
        push (@full_genotypes, "$vcf_complete[0]\t$vcf_complete[1]\t$vcf_complete[3]@");
        my @DP=();
        my @required_meta = ('0', '0', '0', '0', '0', '0');
        # while reading along each line of vcf, starting from the 10th element
        for (my $i = 9; $i <@vcf_complete; $i++) {
            
            # if vcf is not a ./.
            if ($vcf_complete[$i] !~ m/\.\/\./) {
                # split the meta data of the current element
                my @meta = split(":", $vcf_complete[$i]);
                
                #push the 3rd element (DP) into a new array called DP
                push (@DP, $meta[2]);
                
            }
            my @vcf_meta=split(";", $vcf_complete[7]);
            
            for (my $i = 0; $i < @vcf_meta; $i++) {
                if ($vcf_meta[$i] =~ m/^AC=/) {splice @required_meta, 0, 1, substr($vcf_meta[$i], 3)};
                if ($vcf_meta[$i] =~ m/^AF=/) {splice @required_meta, 1, 1, substr($vcf_meta[$i], 3)};
                if ($vcf_meta[$i] =~ m/^AN=/) {splice @required_meta, 2, 1, substr($vcf_meta[$i], 3)};
                if ($vcf_meta[$i] =~ m/^QD=/) {splice @required_meta, 3, 1, substr($vcf_meta[$i], 3)};
                if ($vcf_meta[$i] =~ m/^MQ=/) {splice @required_meta, 4, 1, substr($vcf_meta[$i], 3)};
                if ($vcf_meta[$i] =~ m/^MQ0=/) {splice @required_meta, 5, 1, substr($vcf_meta[$i], 4)};
            }
        }
        push (@meta_anno, (join("\t", @required_meta)),"\t");
        my $cat_full_genotypes = join("\t", (splice @vcf_complete, 9));
        push (@full_genotypes, "$cat_full_genotypes\n");
        my $DP_sum;
        foreach (@DP) { $DP_sum += $_};
        
        my $aveDP = sprintf("%.1f", ($DP_sum/@DP));
        my $minDP = min @DP;
        my $maxDP = max @DP;
        
        push (@meta_anno, "$aveDP\t$minDP-$maxDP\n");
    }
    open(VCFMETA, ">>$input/VCFMETA_anno") or die "error creating VCFMETA_anno file: <$!>\n";
    print VCFMETA @meta_anno;
    close VCFMETA;
    open(VCFFULL, ">>$input/VCFFULL_anno") or die "error creating VCFFULL_anno file: <$!>\n";
    print VCFFULL (@full_genotypes);
    close VCFFULL;
    print STDERR "Done\n\n";
}

sub yn{
    
    print "\n\n*** Do you wish to continue with genotyping (y/n)? ***\n";
    
    sleep (10);
    
    my $yesno;
    
    $yesno=<STDIN>;
    
    chomp $yesno;
    if ($yesno eq "n") {
        cleanup();
        exit 0;
    }
}

sub cleanup{
    
    #create Variants folder using $cohort, loop through array of vcf files and move them
    
    my @final_vcfs = glob("$path2vcf/*");
    mkdir ("$input/$cohort") or die "Unable to create new $cohort directory: <$!>\n";
    mkdir ("$input/$cohort/Variants") or die "Unable to create new $cohort Variants directory: <$!>\n";
    mkdir ("$input/$cohort/Bamfiles") or die "Unable to create new $cohort Bamfiles directory: <$!>\n";
    #@#mkdir ("$input/$cohort/ReducedBamfiles") or die "Unable to create new $cohort ReducedBamfiles directory: <$!>\n";
    mkdir ("$input/$cohort/Logs") or die "Unable to create new $cohort Logs directory: <$!>\n";
    
    for (my $i = 0; $i <@final_vcfs; $i++) {
        my $current_vcf_file = substr $final_vcfs[$i], ($offset + 10);
        my $new_file_extension = "$input/$cohort/Variants/$current_vcf_file";
        move($final_vcfs[$i], $new_file_extension) or die "Unable to move $final_vcfs[$i] to $new_file_extension: <$!>\n";
    }
    
    my @logs = glob("$path2log/*");
    
    for (my $i = 0; $i <@logs; $i++) {
        my $current_log = substr $logs[$i], ($offset + 10);
        my $new_file_extension = "$input/$cohort/Logs/$current_log";
        move($logs[$i], $new_file_extension) or die "Unable to move $logs[$i] to $new_file_extension: <$!>\n";
    }
    
    
    #loop through array of reduced.bam files files and move them
    
    #@#my @reduced_bams = glob("$path2bam/*reduced.bam");
    #@#for (my $i = 0; $i <@reduced_bams; $i++) {
    #@#    my $current_reducedbam_file = substr $reduced_bams[$i], ($offset + 10);
    #@#    my $new_file_extension = "$input/$cohort/ReducedBamfiles/$current_reducedbam_file";
    #@#    move($reduced_bams[$i], $new_file_extension) or die "Unable to move $reduced_bams[$i] to $new_file_extension: <$!>\n";
    #@#}
    
    #loop through array of remaining bam files and index files and move them
    
    my @final_bams = glob("$path2bam/*");
    
    for (my $i = 0; $i <@final_bams; $i++) {
        my $current_bam_file = substr $final_bams[$i], ($offset + 10);
        my $new_file_extension = "$input/$cohort/Bamfiles/$current_bam_file";
        move($final_bams[$i], $new_file_extension) or die "Unable to move $final_bams[$i] to $new_file_extension: <$!>\n";
    }
    
    rmdir "$path2bam" or die "Unable to delete $path2bam folder: <$!>\n";
    rmdir "$path2vcf" or die "Unable to delete $path2vcf folder: <$!>\n";
    rmdir "$path2log" or die "Unable to delete $path2log folder: <$!>\n";
}