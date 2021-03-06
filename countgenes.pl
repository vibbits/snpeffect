#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It is intended to be run as the last step of the SNP analysis pipeline
# to make am overview of which genes were found.

# count mutations (= non-comment lines in input VCF file)
$Nmutations = `grep -v \'^#\' in.vcf | wc -l`;
$Nmutations = $Nmutations * 1; # remove file name from wc output

# count mutations that were deemed suited for further investigation
# (use code from script makeSNPlist4PolyPhen.pl)
open MUT, 'mutated_sequences.fa' or die "cannot open mutated_sequences.fa\n";
while (<MUT>) {
  if (not /^>\d+\t[\w\.]+ Variant (\w+:\d+)-\d+ Ref:([A-Z]) Alt:([A-Z]) HGVS/) {
    die "error in mutated_sequences.fa :\n$_\n";
  }
  $missencemutations{"$1_$2_$3"} = 1;
  <MUT>; # skip sequence
}
close MUT;
$Nsnp = keys %missencemutations;

# make hash with genes found by SnpEff and count them
open IN, 'snpEff_genes.txt' or die "cannot open snpEff_genes.txt\n";
<IN> ; <IN>; # skip header
while (<IN>) {
  $_ =~ /^([^\t]+)\t/;
  $mutated_genes{$1} = 1;
}
close IN;
$Nmutated_genes = keys %mutated_genes;

# make hash with genes for which there is at least one missence mutation
# in not too long protein and count them
open IN, 'SEQANALreport.tab' or die "cannot open SEQANALreport.tab\n";
<IN>; # skip header
while (<IN>) {
  @fields = split;
  $snp_genes{$fields[3]} = 1;
}
close IN;
$Nsnp_genes = keys %snp_genes;

# make hash with genes for which there is at least one missence mutation
# in a known structure and count them
open IN, 'FoldXreport.tab' or die "cannot open FoldXreport.tab\n";
<IN>; # skip header
while (<IN>) {
  @fields = split;
  $toFoldX_genes{$fields[3]} = 1;
}
close IN;
$NtoFoldX_genes = keys %toFoldX_genes;

# adapt the values of %mutated_genes to :
# 1 : mutation but no missence SNP
# 2 : missence SNP but no structure
# 3 : structure
for $key (keys %mutated_genes) {
  if ($snp_genes{$key}) { $mutated_genes{$key}++ }
  if ($toFoldX_genes{$key}) { $mutated_genes{$key}++ }
}

# print report
open OUT, '>finalreport.txt';
print OUT "total number of mutations : $Nmutations\n";
print OUT "total number of genes with variants : $Nmutated_genes\n";
print OUT "\n  We exclude from analysis mutations that are not SNP's, ORF's that are poorly\n  defined (internal stop codons and/or missing parts), SNP's that fall outside\n  an ORF, SNP's that cause a silent mutation or a nonsense mutation and\n  proteins that are too long for analysis.\n\n";
print OUT "number of investigated SNP's : $Nsnp\n";
print OUT "number of genes with SNP : $Nsnp_genes\n";
print OUT "number of genes with SNP in known structure : $NtoFoldX_genes\n";

print OUT "\n";
print OUT "genes for which there are mutations but not SNP's suited for analysis :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 1) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there are SNP's suited for analysis but no related known\n  structure in the PDB files :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 2) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there is a related known structure in the\n  PDB files :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 3) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n";
