#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of the FoldX BuilModel command and makes
# the final report. It also adds some info about chains.
# It takes as input variants4foldX.tab, BuildModel_lastline.tab and
# interactingchains.tab. It writes as output FoldXreport.tab
# (which will have as many records as variants4foldX.tab)

$PDBFILEFOLDER = $ARGV[0];

open IN, 'BuildModel_lastline.tab' or die "cannot open BuildModel_lastline.tab\n";
while (<IN>) {
  chomp;
  if (not /^([A-Za-z0-9,]+;)\t(.+)_1\t(.+)$/) {
    die "error in file BuildModel_lastline.tab :\n$_\n";
  }
  $mutationstring = $1;
  $PDBref = $2;
  $key = "$PDBref.pdb\t$mutationstring";
  $BuilModeloutput{$key} = $3;
}
close IN;

open IN, 'interactingchains.tab' or die "cannot open interactingchains.tab\n";
while (<IN>) {
  chomp;
  if (not /^(.+\.pdb\t[A-Za-z0-9,]+;)\t(.+\t.+\t[-0-9\.]+)$/) {
    die "error in file interactingchains.tab :\n$_\n";
  }
  $chainsinfo{$1} = $2;
}

open IN, 'variants4foldX.tab' or die "cannot open variants4foldX.tab\n";
open OUT, '>FoldXreport.tab';
print OUT "variant\tchromosomal location\tmutation\tgene ID\tgene name\ttranscript ID\tlocation on transcript\tPDB file\tmutationstring\tchains\tinteracting chains\tdelta interaction energy\tSD\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\tElectrostatics\tSolvation Polar\tSolvation Hydrophobic\tVan der Waals clashes\tentropy sidechain\tentropy mainchain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional clash\tbackbone clash\thelix dipole\twater bridge\tdisulfide\telectrostatic kon\tpartial covalent bonds\tenergy Ionisation\tEntropy Complex\n";
while (<IN>) {
  chomp;
  if (not /^\d+\t\w+:\d+-\d+\t[A-Z]+\/[A-Z]+\t.+\t.+[\w\.]+\t\d+\t(.+.pdb)\t(.*)$/) {
    die "error in file variants4foldX.tab :\n$_\n";
  } else {
    $PDBfile = $1;
    $mutationstring = $2;
#    $PDBfilepath = "${PDBFILEFOLDER}/${PDBfile}";
#    $aminoacid = $3;
#    $chain = $4;
#    $number = $5;
#    $tempfactor = fetch_temp_factor($PDBfilepath, $chain, $number);
    $key = "${PDBfile}\t${mutationstring}";
    if (not exists $BuilModeloutput{$key}){
      die "pair $key not found in BuildModel_lastline.tab\n";
    }
    if (not exists $chainsinfo{$key}){
      die "pair $key not found in interactingchains.tab\n";
    }
    print OUT "var_$_\t$chainsinfo{$key}\t$BuilModeloutput{$key}\n";
  }
}
close IN; close OUT;

sub setflag {
  # to indicate that script has terminated, notwithstanding it crashed or
  #   # ended normally ; useful if script is started under a masterscript
  open FLAG, '>flag.PDBTempFactor';
  close FLAG;
}
  #
sub fetch_temp_factor {
  my ($PDBfile, $chain, $number) = @_;
  #           # since we are dealing with structures from AlphaFold, we only take the last match and its temp factor
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
