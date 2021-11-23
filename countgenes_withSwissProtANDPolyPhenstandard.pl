#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It is intended to be run as the last step of the SNP analysis pipeline
# to make am overview of which genes were found.

# count mutations (= non-comment lines in input VCF file)
$Nmutations = `grep -v \'^#\' in.vcf | wc -l`;
$Nmutations = $Nmutations * 1; # remove file name from wc output

# count mutations that were deemed suited for further investigation
$Nsnp = `wc -l SNP4PolyPhen.list`;
$Nsnp = $Nsnp * 1; # remove file name from wc output

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

# make hash with genes that have no UniProt standard
open IN, 'variants_not_investigated.txt'
  or die "cannot open variants_not_investigated.txt\n";
while (<IN>) {
  if (/ : ([^ ]+) does not have a UniProt standard/ and not exists $snp_genes{$1}) {
    $rejected_genes{$1} = 1;
  }
}
close IN;

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

# make hash with genes for which we have reported a PolyPhen standard
open IN, 'SEQANALreport_withPolyPhen.tab' or
  die "cannot open SEQANALreport_withPolyPhen.tab\n";
<IN>; # skip header
while (<IN>) {
  @fields = split;
  $withstandard_genes{$fields[3]} = 1;
}
close IN;
$Nwithstandard_genes = keys %withstandard_genes;

# make hash with genes for which we have reported a PolyPhen standard and a structure
open IN, 'FoldXreport_withPolyPhenonly.tab' or
  die "cannot open FoldXreport_withPolyPhenonly.tab\n";
<IN>; # skip header
while (<IN>) {
  @fields = split;
  $withstandardandstructure_genes{$fields[3]} = 1;
}
close IN;
$Nwithstandardandstructure_genes = keys %withstandardandstructure_genes;

# adapt the values of %mutated_genes to :
# 1 : mutation but no missence SNP
# 2 : missence SNP but no UniProt standard
# 3 : no structure or PolyPhen standard
# 4 : structure but no PolyPhen standard
# 5 : PolyPhen standard but no structure
# X : structure and PolyPhen standard but not in same transcript
#   DOES NOT HAPPEN since we only invertigate one transcript per gene
# 6 : structure and PolyPhen standard in same transcript
for $key (keys %mutated_genes) {
  if ($rejected_genes{$key}) { $mutated_genes{$key}++ }
  if ($snp_genes{$key}) { $mutated_genes{$key} += 2 }
  if ($toFoldX_genes{$key}) { $mutated_genes{$key}++ }
  if ($withstandard_genes{$key}) { $mutated_genes{$key} += 2 }
#  if ($withstandardandstructure_genes{$key}) { $mutated_genes{$key} += 1 }
}
$Nwithstandardandstructure_genes = keys %withstandardandstructure_genes;

# print report
open OUT, '>finalreport.txt';
print OUT "total number of mutations : $Nmutations\n";
print OUT "total number of genes with variants : $Nmutated_genes\n";
print OUT "\n  We exclude from analysis mutations that are not SNP's, ORF's that are poorly\n  defined (internal stop codons and/or missing parts), SNP's that fall outside\n  an ORF, SNP's that cause a silent mutation or a nonsense mutation and\n  proteins that are too long for analysis.\n\n";
print OUT "  To save CPU time we further restrict for each mutation the analysis to\n  one transcript that translates into a UniProt \"standard sequence\".\n\n";
print OUT "number of investigated SNP's : $Nsnp\n";
print OUT "number of genes with SNP : $Nsnp_genes\n";
print OUT "number of genes with SNP in known structure : $NtoFoldX_genes\n";
print OUT "number of genes with SNP in PolyPhen standard transcript : $Nwithstandard_genes\n";
print OUT "number of genes with SNP in PolyPhen standard transcript with known structure : $Nwithstandardandstructure_genes\n";

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
print OUT "genes for which there is no UniProt standard transcript :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 2) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there are SNP's suited for analysis but no related known\n  structure in the SwitchLab PDB files and no PolyPhen standard transcript\n  in the PolyPhen database :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 3) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there is only a related known structure in the SwitchLab\n  PDB files :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 4) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there is only a PolyPhen standard transcript in the\n  PolyPhen database :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 5) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n\n";
print OUT "genes for which there is an SNP suited for analysis in a PolyPhen standard\n  transcript for which there is a related known structure in the SwitchLab\n  PDB files :\n";
$n = 0;
for $key (sort keys %mutated_genes) {
  if ($mutated_genes{$key} == 6) {
    print OUT " $key";
    $n++;
    if ($n % 5 == 0) { print OUT "\n" } # new line each 5 genes
  }
}
print OUT "\n";
