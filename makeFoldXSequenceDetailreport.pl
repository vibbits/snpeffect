#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of the FoldX AnalyseComplex command and makes
# the final report.
# It takes as input variants4foldX.tab and SequenceDetail_mutlines.tab
# and writes as output FoldXreport_SequenceDetail.tab (which will have as
# many variants as variants4foldX.tab but can have several records
# for one variant if the protein has multiple chains).
open IN, 'SequenceDetail_mutlines.tab' or die "cannot open SequenceDetail_mutlines.tab\n";
while (<IN>) {
  chomp;
  if (not /^([A-Za-z0-9,;]+)\t(RepairPDB_\d+\.pdb)\t(.+)$/) {
    die "error in file AnalyseComplex_meandiff.tab :\n$_\n";
  } else {
    $mutationstring = $1;
    $PDBfile = $2;
    $key = "$PDBfile\t$mutationstring";
    $rest_of_line = $3;
    $datablock{$key} .= "R_E_C_O_R_D$rest_of_line";
      # dumb tag R_E_C_O_R_D is to be replaced later
  }
}
close IN;
#print $datablock{"RepairPDB_28403.pdb\tFA113C;"};

open IN, 'variants4foldX.tab' or die "cannot open variants4foldX.tab\n";
open OUT, '>FoldXreport_SequenceDetail.tab';
print OUT "variant\tchromosomal location\tmutation\tgene ID\tgene name\ttranscript ID\tposition mutation\tPDB file\tamino acid\tchain\tnumber\tomega angle\tphi angle\tpsi angle\tSecondary Structure\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\tElectrostatics\tSolvation Polar\tSolvation Hydrophobic\tVan der Waals clashes\tentropy sidechain\tentropy mainchain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional clash\tbackbone clash\thelix dipole\twater bridge\tdisulfide\telectrostatic kon\tpartial covalent bonds\tenergy Ionisation\tEntropy Complex\theteroBackHbond\theteroSideHbond\tsidechain burial\tmainchain burial\tsidechain Occ\tmainchain Occ\tindex";
while (<IN>) {
  chomp;
  if (not /^\d+\t\w+:\d+-\d+\t[A-Z]+\/[A-Z]+\t.+\t.+[\w\.]+\t\d+\t(RepairPDB_\d+.pdb)\t(.*)$/) {
    die "error in file variants4foldX.tab :\n$_\n";
  } else {
    $PDBfile = $1;
    $mutationstring = $2;
    $key = "${PDBfile}\t${mutationstring}";
    if (exists $datablock{$key}){
      $_ =~ s/\t[A-Za-z0-9,]+;$//;
        # we do not need the mutationstring since there are already
        # amino acid and chain fields
      $datablock = $datablock{$key};
      $datablock =~ s/R_E_C_O_R_D/\nvar_$_\t/g;
      print OUT $datablock;
    } else {
      die "pair $key not found in SequenceDetail_mutlines.tab\n";
    }
  }
}
print OUT "\n";
close IN; close OUT;
