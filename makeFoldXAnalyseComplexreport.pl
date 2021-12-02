#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of the FoldX AnalyseComplex command and makes
# the final report.
# It takes as input variants4foldX.tab and AnalyseComplex_meandiff.tab
# and writes as output FoldXreport_AnalyseComplex.tab (which does not
# contain entries for each variant since not all proteins have interacting
# chains)
open IN, 'AnalyseComplex_meandiff.tab' or die "cannot open AnalyseComplex_meandiff.tab\n";
while (<IN>) {
  chomp;
  if (not /^(.+\.pdb\t[A-Za-z0-9,;]+)\t(.\t.\t.+)$/) {
    die "error in file AnalyseComplex_meandiff.tab :\n$_\n";
  } else {
    $key = $1;
    $rest_of_line = $2;
    $datablock{$key} .= "R_E_C_O_R_D$rest_of_line";
  }
}
close IN;

open IN, 'variants4foldX.tab' or die "cannot open variants4foldX.tab\n";
open OUT, '>FoldXreport_AnalyseComplex.tab';
print OUT "variant\tchromosomal location\tmutation\tgene ID\tgene name\ttranscript ID\tposition mutation\tPDB file\tmutationstring\tChain1\tChain2\tdelta IntraclashesChain1\tdelta IntraclashesChain2\tdelta Interaction Energy\tdelta StabilityChain1\tdelta StabilityChain2";
while (<IN>) {
  chomp;
  if (not /^\d+\t\w+:\d+-\d+\t[A-Z]+\/[A-Z]+\t.+\t.+[\w\.]+\t\d+\t(.+.pdb)\t(.*)$/) {
    die "error in file variants4foldX.tab :\n$_\n";
  } else {
    $PDBfile = $1;
    $mutationstring = $2;
    $key = "${PDBfile}\t${mutationstring}";
    if (exists $datablock{$key}){
      $datablock = $datablock{$key};
      $datablock =~ s/R_E_C_O_R_D/\nvar_$_\t/g;
      print OUT $datablock;
    }
  }
}
print OUT "\n";
close IN;
