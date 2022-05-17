#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of the FoldX AnalyseComplex command and makes
# the final report.
# It takes as input variants4foldX.tab and SequenceDetail_mutlines.tab
# and writes as output FoldXreport_SequenceDetail.tab (which will have as
# many variants as variants4foldX.tab but can have several records
# for one variant if the protein has multiple chains).

$PDBFILEFOLDER = $ARGV[0];

open IN, 'SequenceDetail_mutlines.tab' or die "cannot open SequenceDetail_mutlines.tab\n";
while (<IN>) {
  chomp;
  if (not /^([A-Za-z0-9,;]+)\t(.+\.pdb)\t(\w+)\t([A-Z]+)\t([0-9]+)\t(.+)$/) {
    die "error in file SequenceDetail_mutlines.tab :\n$_\n";
  } else {
    $mutationstring = $1;
    $PDBfile = $2;
    $PDBfilepath = "${PDBFILEFOLDER}/${PDBfile}";
    $aminoacid = $3;
    $chain = $4;
    $number = $5;
    $tempfactor = fetch_temp_factor($PDBfilepath, $chain, $number);
    #print "pLDDT is: ${tempfactor}\n";
    $key = "$PDBfile\t$mutationstring";
    print "${key}\n";
    # add temp factor at end of the line
    $rest_of_line = "$6\t${tempfactor}";
    $datablock{$key} .= "R_E_C_O_R_D$rest_of_line";
      # dumb tag R_E_C_O_R_D is to be replaced later
  }
}
close IN;
#print $datablock{"AF-O43869-F1-model_v1_Repair.pdb\tIA142L;"};

open IN, 'variants4foldX.tab' or die "cannot open variants4foldX.tab\n";
open OUT, '>FoldXreport_SequenceDetail.tab';
print OUT "variant\tchromosomal location\tmutation\tgene ID\tgene name\ttranscript ID\tposition mutation\tPDB file\tomega angle\tphi angle\tpsi angle\tSecondary Structure\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\tElectrostatics\tSolvation Polar\tSolvation Hydrophobic\tVan der Waals clashes\tentropy sidechain\tentropy mainchain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional clash\tbackbone clash\thelix dipole\twater bridge\tdisulfide\telectrostatic kon\tpartial covalent bonds\tenergy Ionisation\tEntropy Complex\theteroBackHbond\theteroSideHbond\tsidechain burial\tmainchain burial\tsidechain Occ\tmainchain Occ\tindex\tpLDDT";
while (<IN>) {
  chomp;
  if (not /^\d+\t\w+:\d+-\d+\t[A-Z]+\/[A-Z]+\t.+\t.+[\w\.]+\t\d+\t(.+.pdb)\t(.*)$/) {
    die "error in file variants4foldX.tab :\n$_\n";
  } else {
    $PDBfile = $1;
    $mutationstring = $2;
    $key = "${PDBfile}\t${mutationstring}";
    print "${key}\n";
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

sub setflag {
  # to indicate that script has terminated, notwithstanding it crashed or
  # ended normally ; useful if script is started under a masterscript
  open FLAG, '>flag.PDBTempFactor';
  close FLAG;
}

sub fetch_temp_factor {
  my ($PDBfile, $chain, $number) = @_;
  # since we are dealing with structures from AlphaFold, we only take the last match and its temp factor
  open INPDB, "${PDBfile}";
  while (<INPDB>) {
      chomp;
      my $line = $_;
      if ($line =~ m/^ATOM\s+\d+\s+\w+\s+(\w{3})\s+(\w)\s?(\s?\s?\d+)\s+(\S+\.\S+)\s+(\S+\.\S+)\s+(\S+\.\S+)\s+.+\..+\..+/ig) {
         $id = $3;
         $id =~ s/^\s+//;
         if ($2 eq $chain & $id eq $number & length($line) == 66 ) {
                $tempfactor = substr($line, -6);
                $tempfactor =~ s/^\s+//;
           } 
      }
  } 
  if (not $tempfactor) {
     &setflag;
     die "Amino acid ${number} of ${chain} not found in ${PDBfile}\n";
  }
  close INPDB;

  return $tempfactor;
}
