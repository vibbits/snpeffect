#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It lets FoldX investigate the effect of mutations on protein structure.
# It can best be run on the cluster
# It takes as input :
#   - the file protvariants4foldX.tab with metadata produced by
#     parseBLASToutput.pl and parsevariants4foldX.pl
#   - files FoldXoptions, testFoldXcommands and FoldXcommands4AnalyseComplex
#   - the PDB files in /switchlab/group/robkan/Homology_Human/RepairPDBs
# FoldX works only properly if it is run in the folder where the PDB file
# is. There it will also write the output. This folder must contain
# rotabase.txt (or a logical link to it).
# By default it deletes the unneeded FoldX output because we would otherwise
# fill the disk for anything but a tiny data set. It extracts the information
# that we need for the further steps of the analysis and stores them in
# temporary files :
#   - BuildModel_lastline.tab
#   - SequenceDetail_mutlines.tab
#   - AnalyseComplex_meandiff.tab
#   - interactingchains.tab

use File::Copy;

$delete = 1; # delete FoldX output except the needed files.
$Nruns = 3; # number of FoldX BuildModel runs, as set in "FoldXoptions"

$FoldX = $ARGV[0];
#$FoldX = '/switchlab/group/tools/FoldX_2015/FoldX';
if (not $FoldX or $FoldX !~ /FoldX$/) {
  &setflag;
  die "script submit2FoldX.pl needs FoldX as first argument\n"
}
$PDBdir = $ARGV[1];
#$PDBdir = '/switchlab/group/robkan/Homology_Human/RepairPDBs';
if (not $PDBdir) {
  &setflag;
  die "script submit2FoldX.pl needs as second argument folder with PDB files\n"
}

# test if we have everything needed
if (not -e 'protvariants4foldX.tab')
  { &setflag ; die "protvariants4foldX.tab not found\n" }
if (not -e '../FoldXoptions') { &setflag ; die "FoldXoptions not found\n" }
if (not -e '../FoldXcommands') { &setflag ; die "FoldXcommands not found\n" }
if (not -e '../FoldXcommands4AnalyseComplex')
  { &setflag ; die "FoldXcommands4AnalyseComplex not found\n" }
if (not -e "$FoldX") { &setflag ; die "$FoldX not found\n" }
if (not -e 'rotabase.txt') { &setflag ; die "rotabase.txt found\n" }

%aa1 = ('ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N', 'ASP' => 'D', 'CYS' => 'C',
  'PHE' => 'F', 'GLN' => 'Q', 'GLU' => 'E', 'GLY' => 'G', 'HIS' => 'H',
  'ILE' => 'I', 'LEU' => 'L', 'LYS' => 'K', 'MET' => 'M', 'PRO' => 'P',
  'SER' => 'S', 'THR' => 'T', 'TRP' => 'W', 'TYR' => 'Y', 'VAL' => 'V',
  # unusal amino acids supported by foldX
  'H1S' => 'H', 'H2S' => 'H', 'H3S' => 'H', 'HYP' => 'P', 'M3L' => 'K',
  'MLY' => 'K', 'MLZ' => 'K', 'PTR' => 'T', 'SEP' => 'S', 'TPO' => 'Y',
  'TYS' => 'Y');


# make arrays with respectively PDB file and mutation
open IN, 'protvariants4foldX.tab';
while (<IN>) {
  if (not /^(RepairPDB_\d+.pdb)\t(.*)$/) {
    &setflag;
    die "error in file protvariants4foldX.tab :\n$_\n";
  } else {
    $N++;
    $PDBfile[$N] = $1;
    $mutationstring[$N] = $2;
  }
}
close IN;

# run FoldX and store useful output in temporary files
open BUILDMODEL, '>BuildModel_lastline.tab';
open SEQUENCEDETAIL, '>SequenceDetail_mutlines.tab';
open ANALYSECOMPLEX, '>AnalyseComplex_meandiff.tab';
open CHAINS, '>interactingchains.tab';
for ($i = 1 ; $i <= $N ; $i++) {
  open OUT, '>individual_list.txt';
  print OUT "$mutationstring[$i]\n";
  close OUT;
  $localPDBfile = "$PDBfile[$i]";
  if (not -e "$PDBdir/$PDBfile[$i]") {
    &setflag;
    die "file $PDBdir/$PDBfile[$i] does not exist\n";
  }
  copy("$PDBdir/$PDBfile[$i]", $localPDBfile);
  $cmd = "$FoldX -manual $localPDBfile ../FoldXoptions ../FoldXcommands >> LogFile";
  system $cmd;

  $specifier = $localPDBfile; $specifier =~ s/.pdb$//;
    # specifier is local name PDB file without extension .pdb

  # parse BuildModel output
  open IN, "Average_BuildModel_$specifier.fxout";
  @lines = <IN>;
  $lastline = $lines[$#lines];
  print BUILDMODEL "$mutationstring[$i]\t$lastline";
  close IN;

  # parse SequenceDetail output
  $mutationstring = $mutationstring[$i];
  chop $mutationstring; # remove trainling ';'
  @mutations = split /,/, $mutationstring;
  @mutations = sort @mutations;
    # sort to have the chains in the order A B C...
    # However, some entries have chains that have a non A-Z symbol
    #   and/or are not in order, hence we cannot parse just once
    #   through the file but better do it separately for each chain
  foreach $m (@mutations) {
    $m =~ /^([A-Z])(.)(\d+)[A-Z]$/;
    $AA = $1 ; $CHAIN = $2 ; $POS = $3;
    open IN, "SequenceDetail_$specifier.fxout";
    $parsingheader = 1;
    while (<IN>) {
      if ($parsingheader) {
        if (/^Pdb\tamino acid/) {
          $parsingheader = 0;
        }
      } elsif (not /^$/) { # note : there is an empty line at end file
        if ($_ !~ /^\w+\.pdb\t(\w+)\t/) {
          &setflag;
          die "error in SequenceDetail_$specifier.fxout line\n $_";
        }
        if (exists $aa1{$1}) {
          if ($_ !~ /^\w+\.pdb\t(\w+)\t(.)\t(-?\d+)/) {
            # note : some PDB files have negative pos numbers
            &setflag;
            die "error in SequenceDetail_$specifier.fxout line\n $_";
          }
          $aa = $aa1{$1} ; $chain = $2 ; $pos = $3;
          if ($POS == $pos and $CHAIN eq $chain) {
            if ($AA ne $aa) {
              &setflag;
              die "amino acid mismatch in SequenceDetail_$specifier.fxout line\n $_";
            }
            print SEQUENCEDETAIL "$mutationstring[$i]\t";
            print SEQUENCEDETAIL;
          }
        }
      }
    }
    close IN;
  }

  # parse AnalyseComplex output to check which chains interact with each
  # other. Note that if the PDB file contains only one chain there will be no
  # AnalyseComplex output.
  print CHAINS "$PDBfile[$i]\t$mutationstring[$i]";
  %chains = () ; %interactingchains = ();
  if (-e "Interface_Residues_AnalyseComplex_$specifier.fxout") {
    open IN, "Interface_Residues_AnalyseComplex_$specifier.fxout";
    $parsingheader = 1;
    while (<IN>) {
      if ($parsingheader) {
        if (/^$specifier\.pdb interface residues between (.) and (.)$/) {
          $parsingheader = 0;
          $chain1 = $1 ; $chain2 = $2;
          $chains{$chain1} = 1 ; $chains{$chain2} = 1;
          $residueline = <IN>;
          if ($residueline !~ /^$/) {
            $interactingchains{$chain1} = 1 ; $interactingchains{$chain2} = 1;
            $interactingchainpairs{"$chain1\t$chain2"} = 1;
          }
        }
      } else {
        if (/^$specifier\.pdb interface residues between (.) and (.)$/) {
          $chain1 = $1 ; $chain2 = $2;
          $chains{$chain1} = 1 ; $chains{$chain2} = 1;
          $residueline = <IN>;
          if ($residueline !~ /^$/) {
            $interactingchains{$chain1} = 1 ; $interactingchains{$chain2} = 1;
            $interactingchainpairs{"$chain1\t$chain2"} = 1;
          }
        } else {
          &setflag;
          die "error in line \n$_ of Interface_Residues_AnalyseComplex_$specifier.fxout\n";
        }
      }
    }
    close IN;
    $chains = join '', sort keys %chains;
    $interactingchains = join '', sort keys %interactingchains;
  } else {
    $chains = 'A' ; $interactingchains = '';
    print CHAINS "\tA\t.\t.\n";
  }

  # if there are interacting chains, then run FoldX AnalyseComplex on the
  # 3 alternative structures computed by FoldX BuildModel and compute
  # the average change in interaction energy between the wild type and
  # the mutant.
  %chaininteraction_records = ();
  if ($interactingchains) {
    print CHAINS "\t$chains\t$interactingchains";
    $total_interactionenergy_change = 0;
    for ($j = 0 ; $j < $Nruns ; $j++) {
      &run_analyse_complex('', $j, 1);
      &run_analyse_complex('WT_', $j, -1);
    }
    for $k (sort keys %chaininteraction_records) {
      print ANALYSECOMPLEX "$PDBfile[$i]\t$mutationstring[$i]\t";
      print ANALYSECOMPLEX $k;
      for ($J = 2 ; $J <= 6 ; $J++) {
        printf ANALYSECOMPLEX '%s%.4f', "\t",
          $chaininteraction_records{"$k"}[$J] / 3;
      }
      print ANALYSECOMPLEX "\n";
      $total_interactionenergy_change += $chaininteraction_records{"$k"}[4] / 3;
    }
    printf CHAINS  '%s%.4f%s', "\t", $total_interactionenergy_change, "\n";
  } elsif ($chains ne 'A') {
    print CHAINS "\t$chains\t.\t.\n";
  }

  # delete local PDB file and FoldX output
  if ($delete) {
    opendir DIR, '.';
    @FoldXfiles =  readdir DIR;
    foreach $f (@FoldXfiles) {
      if ($f =~ /${specifier}\D/) {
        unlink $f;
      }
    }
    closedir DIR;
  } else {
    unlink $localPDBfile;
  }
}

# terminate
unlink 'individual_list.txt', 'LogFile';
&setflag;

sub setflag {
  # to indicate that script has terminated, notwithstanding it crashed or
  # ended normally ; useful if script is started under a masterscript
  open FLAG, '>flag.FoldX';
  close FLAG;
}

sub run_analyse_complex {
  my ($prefix, $n, $x) = @_;
  $cmd = "$FoldX -manual ${prefix}${specifier}_1_$n.pdb ../FoldXoptions ../FoldXcommands4AnalyseComplex >> LogFile";
  system $cmd;

  open IN, "Summary_AnalyseComplex_${prefix}${specifier}_1_$n.fxout";
  my $parsingheader = 1;
  while (<IN>) {
    if ($parsingheader) {
      if (/^Pdb\tGroup1\tGroup2\t/) {
        $parsingheader = 0;
      }
    } else {
      if (/^${prefix}${specifier}_1_$n.pdb\t(.)\t(.)\t([\d\-\.e]+)\t([\d\-\.e]+)\t([\d\-\.e]+)\t([\d\-\.e]+)\t([\d\-\.e]+)$/) {
        # Group1 Group2 IntraclashesGroup1 IntraclashesGroup2 Interaction_Energy StabilityGroup1 StabilityGroup2
        if (exists $interactingchainpairs{"$1\t$2"}) {
          $chaininteraction_records{"$1\t$2"}[0] = $1;
          $chaininteraction_records{"$1\t$2"}[1] = $2;
          $chaininteraction_records{"$1\t$2"}[2] += $3 * $x;
          $chaininteraction_records{"$1\t$2"}[3] += $4 * $x;
          $chaininteraction_records{"$1\t$2"}[4] += $5 * $x;
          $chaininteraction_records{"$1\t$2"}[5] += $6 * $x;
          $chaininteraction_records{"$1\t$2"}[6] += $7 * $x;
        }
      } else {
        &setflag;
        die "error in line \n$_ of ${prefix}Summary_AnalyseComplex_${specifier}_1_$n.fxout\n";    }
    }
  }
  close IN;
}
