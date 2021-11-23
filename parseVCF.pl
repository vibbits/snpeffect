#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It takes as input the file out.vcf made by SnpEff and the file
# variants.tab made by parseSnpEfffasta.pl and makes a file SNPAAchange.vcf
# in VCF format with only the records selected for further analysis.

# make a hash (%Var) with as key the position of the variant on the genome
# and as value the number of the variant
open VAR, 'variants.tab' or "die cannot open variants.tab\n";
while (<VAR>) {
  if (not /^(\d+)\t(\w+):(\d+)-(\d+)\t[A-Z]+\/[A-Z]+\t.+\t.+\t[\w\.]+\t[A-Z]\t\d+\t[A-Z]$/) {
    die "error in file variants.tab :\n$_\n";
  }
  $Nvar = $1;
  $chromosome = $2;
  $chrompos = $3;
  #$chromposend = $4;
  #if ($chrompos != $chromposend) {
  #  die "variants.tab contains a non-SNP mutation :\n$_\n";
  #} erroneous test because some 1 AA mutations involve 2 nucleotide mutations
  $Var{"$chromosome\t$chrompos"} .= ",var_$Nvar";
}
close VAR;
foreach $key (keys %Var) {
  # remove the ',' at the begin of the string
  $Var{$key} =~ s/^,//;
}


# make a VCF file with only lines corresponding to variants in %Var
open IN, 'out.vcf' or "die cannot open out.vcf\n";
open OUT, '>SNPAAchange.vcf';
while (<IN>) {
  if (/^#/) { # reading header
    print OUT;
    if (/^##INFO=<ID=NMD,Number/) {
      print OUT "##INFO=<ID=SNPVAR,Number=.,Type=String,Description=\"Numbers selected variants\">\n";
    }
  } else {
    # note that the header does not contain tabs and hence is not modified
    @fields = split /\t/;
      # fields are CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ...
    $fields[0] =~ s/^chr//;
      # chromosome name in out.vcf can start with "chr", name in out.fa does not
    $key = "$fields[0]\t$fields[1]";
    if (exists $Var{$key}) {
      $fields[7] = "SNPVAR=$Var{$key};$fields[7]";
      $line = join "\t", @fields;
      print OUT $line;
    }
  }
}
close IN;
close OUT;
