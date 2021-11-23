#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It searches for (predicted) transmembrane domains in protein sequences,
# using the TMHMM software. As part of the SNP pipeline, it also checks
# whether the mutation fall in a TM and whether the mutated protein
# domain has TMs.
#
# NOTE : the TMHMM software must be customized. In tmhmm and tmhmmformat.pl
# the PATH to perl must be set corrctly and in tmhmm $opt_short must be set
# to 1 and $opt_plot to 0.
#
# It takes as input :
# - variants.tab
# - reference_sequences_nonredundant.fa
# - mutated_domains.tab, with the position of the domain that contains
#   the mutation
# It produces as output :
# - TMHMMoutput.tab
# - TMcheck.tab (to be incorporated into final report)

$tmhmm = $ARGV[0];
#$tmhmm = '/switchlab/group/guybot/tmhmm-2.0c/bin/tmhmm';
if (not $tmhmm or $tmhmm !~ /[^\/]\/tmhmm$/) {
  die "script findTMs.pl needs tmhmm as first argument\n"
}

# run TMHMM on a nonredundant set of sequences and make a hash with the
# transmembrane topology string in the output
open REF, 'reference_sequences_nonredundant.fa' or die "cannot open reference_sequences_nonredundant.fa\n";
open OUT, '>TMHMMoutput.tab';
while(<REF>) {
  if (not /^>\d+\t([\w\.]+)/) {
    die "$_ does not look like fastA description line for transcript ID\n";
  }
  $transcriptID = $1;
  open SEQ, '>tempseq.fa';
  print SEQ ">$transcriptID\n";
  $line = <REF>;
  print SEQ $line;
  close SEQ;

  $outputline = `$tmhmm tempseq.fa`;
  print OUT $outputline;
  $outputline =~ /Topology=(.+)$/;
  $topologystring{$transcriptID} = $1;
}
close REF ; close OUT;
unlink 'tempseq.fa';

# read in information about domain that contains mutation
open IN, 'mutated_domains.tab' or die "cannot open mutated_domains.tab\n";
while (<IN>) {
  if (not /^(\d+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)$/) {
    die "error in file mutated_domains.tab :\n$_\n";
  } else { # $1 is Nvar
    $begindomain{$1} = $4;
    $enddomain{$1} = $5;
  }
}
close IN;

# parse the TMHMM output and make TM checklist
open VAR, 'variants.tab' or "die cannot open variants.tab\n";
open CHECK, '>TMcheck.tab';
while (<VAR>) {
  if (not /^(\d+)\t(\w+:\d+-\d+)\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
    die "error in file variants.tab :\n$_\n";
  }
  $Nvar = $1;
  $transcriptID = $5;
  $seqpos = $7;

  print CHECK "$Nvar";

  @beginTM = () ; @endTM = ();
  while ($topologystring{$transcriptID} =~ /(\d+)-(\d+)/g) {
    push @beginTM, $1;
    push @endTM, $2;
  }
  $NTMs = @beginTM;

  # check whether the protein has a TM
  if ($NTMs > 0) {
    print CHECK "\t+";
  } else {
    print CHECK "\t-";
  }
  # check  whether the mutation falls in a protein domain that has TMs
  if ($NTMs > 0 and exists $begindomain{$Nvar}) {
    $domain_has_TM = 0;
    for ($i = 0 ; $i < $NTMs ; $i++) {
      if ($beginTM[$i] >= $begindomain{$Nvar} and $endTM[$i] <= $enddomain{$Nvar}) {
        $domain_has_TM = 1;
      }
    }
    if ($domain_has_TM){
      print CHECK "\t+";
    } else {
      print CHECK "\t-";
    }
  } else {
    print CHECK "\t-";
  }
  # check whether the mutation falls in a TM
  if ($NTMs > 0) {
    $mutation_is_in_TM = 0;
    for ($i = 0 ; $i < $NTMs ; $i++) {
      if ($seqpos >= $beginTM[$i] and $seqpos <= $endTM[$i]) {
        $mutation_is_in_TM = 1;
      }
    }
    if ($mutation_is_in_TM) {
      print CHECK "\t+\n";
    } else {
      print CHECK "\t-\n";
    }
  } else {
    print CHECK "\t-\n";
  }
}

