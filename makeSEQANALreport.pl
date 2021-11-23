#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It makes the final report of a SNP on protein structure and aggregation
# study (excluding analyses involving FoldX which can only be run if the
# structure is known). It puts the outputs of Gene3D, AGADIR, TMHMM and
# the mutation effect predictors in one file and, while doing so, makes a
# second pass through the TANGO and WALTZ output to  add info about missing
# regions. It also computes the change in TANGO and WALTZ scores between
# reference and mutant.
# It takes as input :
# - mutated_domains.txt
# - domain_differences_AGADIR.tab
# - region_differences_AGADIR.tab
# - domain_regions_differences_AGADIR.tab
# - the files with TANGO and WALTZ output for reference and mutant
# - TMcheck.tab
# - SIFT_PROVEAN.tab
# It can optionally add a column with the UniProt standard sequence ID

$SwissProtstandard = $ARGV[0];
#$SwissProtstandard = '../scripts/NM_AC_ID.tab';

# make hash with domain info in function of variant number
# the missing numbers correspond to variants without domain info
open IN, 'mutated_domains.tab' or die "cannot open mutated_domains.tab\n";
while (<IN>) {
  if (not /^(\d+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)$/) {
    die "error in file mutated_domains.tab :\n$_\n";
  } else { # $1 is Nvar
    $domainID = $2;
    $familyID = $3;
    $begindomain = $4;
    $enddomain = $5;
    $domain_info{$1} = "$domainID\t$familyID\t$begindomain\t$enddomain";
    $position_fromdomaintable{$1} = $6;
    $domain_length{$1} = $enddomain - $begindomain + 1;
  }
}
close IN;

# make hash with PROVEAN and SIFT predictions in function
# of transcript ID and mutation
open IN, 'SIFT_PROVEAN.tab' or die "cannot open SIFT_PROVEAN.tab\n";
while (<IN>) {
  if (not /^([\w\.]+)\t([A-Z]\d+[A-Z])\t([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^t]+)\n$/) {
    die "error in file SIFT_PROVEAN.tab :\n$_\n";
  } else {
    $transcriptID = $1;
    $mutationstring = $2;
    $predictionstring = $3;
    $key = "$transcriptID\t$mutationstring";
    $prediction{$key} = $predictionstring;
  }
}
close IN;

# if needed, make hashes with gene name and UniProt standard sequence ID
# in function of key the RefSeq NP_* ID, this to add a UniProt ID column
if ($SwissProtstandard) {
  open IN, $SwissProtstandard or die "cannot open file $SwissProtstandard\n";
  while (<IN>) {
    chomp;
    $_ =~ /^(.+)\t(.+)\t(.+)$/;
    $standard4refseqID{$1} = $2;
    $genefield = $3;
    @items = split '/', $genefield;
    for $item (@items) {
      $standard4gene{$item} = $2;
    }
  }
  close IN;
}

open VAR, 'variants.tab' or die "cannot open variants.tab\n";
open DOMD, 'domain_differences_AGADIR.tab'
  or die "cannot open domain_differences_AGADIR.tab\n";
open REGD, 'region_differences_AGADIR.tab'
  or die "cannot open region_differences_AGADIR.tab\n";
open DOMREGD, 'domain_regions_differences_AGADIR.tab'
  or die "cannot open domain_regions_differences_AGADIR.tab\n";
open TANGOREF, 'reference_TANGO.tab' or die "cannot open reference_TANGO.tab\n";
open WALTZREF, 'reference_WALTZ.tab' or die "cannot open reference_WALTZ.tab\n";
open TANGOMUT, 'mutated_TANGO.tab' or die "cannot open mutated_TANGO.tab\n";
open WALTZMUT, 'mutated_WALTZ.tab' or die "cannot open mutated_WALTZ.tab\n";
open TM, 'TMcheck.tab' or die "cannot open TMcheck.tab\n";
open OUT, '>SEQANALreport.tab';

# print header
print OUT "variant\tchromosomal location\tmutation\tgene ID\tgene name\t";
print OUT "transcript ID\treference aa\tposition in transcript\tmutant aa\t";
print OUT "CATH domain ID\tCATH family ID\tbegin domain\tend domain\t";
print OUT "protein has TM\tdomain has TM\tmutation in TM\t";
print OUT "TANGO score reference\tTANGO score mutant\t";
print OUT "delta TANGO score\t";
print OUT "WALTZ score reference\tWALTZ score mutant\t";
print OUT "delta WALTZ score\t";
print OUT "TANGO score domain reference\tTANGO score domain mutant\t";
print OUT "delta TANGO score domain\tnormalized TANGO score domain mutant\t";
print OUT "WALTZ score domain reference\tWALTZ score domain mutant\t";
print OUT "delta WALTZ score domain\t";
print OUT "TANGO score regions in domain reference\tTANGO score regions in domain mutant\t";
print OUT "delta TANGO score regions in domain\t";
print OUT "WALTZ score regions in domain reference\tWALTZ score regions in domain mutant\t";
print OUT "delta WALTZ score regions in domain\t";
print OUT "position TANGO region reference\tTANGO region reference\tTANGO score region reference\t";
print OUT "position TANGO region mutant\tTANGO region mutant\tTANGO score region mutant\t";
print OUT "delta TANGO score region\t";
print OUT "position TANGO region highest score in domain\tTANGO region highest score in domain\thighest TANGO score in domain\t";
print OUT "position WALTZ region reference\tWALTZ region reference\tWALTZ score region reference\t";
print OUT "position WALTZ region mutant\tWALTZ region mutant\tWALTZ score region mutant\t";
print OUT "delta WALTZ score region\t";
print OUT "SIFT prediction\tSIFT score\tSIFT median\tSIFT Nsequences\t";
print OUT "PROVEAN prediction\tPROVEAN score";
if ($SwissProtstandard) { print OUT "\tUniProt standard sequence" }
print OUT "\n";

# parse through files and put data together
while(<VAR>) {
  chomp;
  if (not /^(\d+)\t(\w+:\d+-\d+)\t[A-Z]+\/[A-Z]+\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
  # Nvar chrom_pos geneID genename seqID refAA seqpos altAA
  die "error in file variants.tab :\n$_\n";
  }
  $Nvar = $1;
  $seqpos = $7;
  $geneID = $4;
  $transcriptID = $5;
  $predictionkey = "$5\t$6$7$8";
    # needed to add SIFT PROVEAN output
  print OUT "var_$_";

  if (exists $position_fromdomaintable{$Nvar}) {
    if ($seqpos != $position_fromdomaintable{$Nvar}) {
      die "variants.tab and mutated_domains.tab position mismatch for variant $Nvar\n";
    }
    print OUT "\t$domain_info{$Nvar}";
  } else {
    print OUT "\t.\t.\t.\t.";
  }

  $line = <TM>;
  if ($line !~ /^(\d+)\t(.\t.\t.)$/) {
    die "error in file TMcheck.tab :\n$line";
  }
  if ($1 != $Nvar) {
    die "variant number mistmatch between variants.tab $Nvar and TMcheck.tab $1\n";
  }
  print OUT "\t$2";

  $line = <REGD>;
  if ($line !~ /^(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([\d\.]+)\t([A-Z\.]+)\t([-e\d\.]+)\t([\d\.]+)\t([A-Z\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([\d\.]+)\t([A-Z\.]+)\t([-e\d\.]+)\t([\d\.]+)\t([A-Z\.]+)\t([-e\d\.]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\n$/) { die "error in region_differences_AGADIR.tab line\n$line\n" }
  if ($Nvar != $1) { die "$Nvar in variants.tab does not have a corresponding line in region_differences_AGADIR.tab\n" }
  $tangorefscore = $2 ; $tangomutscore = $3;
  $tangorefposition = $4 ; $tangorefpeptide = $5 ; $tangorefregionscore = $6;
  $tangomutposition  = $7 ; $tangomutpeptide = $8 ; $tangomutregionscore = $9;
  $waltzrefscore = $10 ; $waltzmutscore = $11;
  $waltzrefposition = $12 ; $waltzrefpeptide = $13 ; $waltzrefregionscore = $14;
  $waltzmutposition = $15 ; $waltzmutpeptide = $16 ; $waltzmutregionscore = $17;
  $Ntangoref = $18 ; $Ntangomut = $19 ; $Nwaltzref = $20 ; $Nwaltzmut = $21;

  $line = <DOMREGD>;
  if ($line !~ /^(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([^\t]+)\t([^\t]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([\d\.]+)\t([A-Z\.]+)\t([-e\d\.]+)\n$/) { die "error in domain_regions_differences_AGADIR.tab line\n$line\n" }
  #if ($line !~ /^(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([^\t]+)\t([^\t]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\t(\d+)\t([A-Z]+)\t([-e\d\.]+)\n$/) { die "error in domain_regions_differences_AGADIR.tab line\n$line\n" }
  if ($Nvar != $1) { die "$Nvar in variants.tab does not have a corresponding line in domain_regions_differences_AGADIR.tab\n" }
  if ($2 ne $tangorefscore or $3 ne $tangomutscore or $4 ne $waltzrefscore or $5 ne $waltzmutscore) { die "line\n$line\nin domain_regions_differences_AGADIR.tab does not contain the same global TANGO/WALTZ score as the corresponding line in region_differences_AGADIR.tab\n"}
  $domainID = $6 ; $familyID = $7;
  $tangorefregionSscore = $8 ; $tangomutregionSscore = $9;
  $waltzrefregionSscore = $10 ; $waltzmutregionSscore = $11;
  $tangomutMAXregionSposition = $12;
  $tangomutMAXregionS = $13;
  $tangomutMAXregionSscore = $14;

  $line = <DOMD>;
  if ($line !~ /^(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\n$/) { die "error in domain_differences_AGADIR.tab line\n$line\n" }
  if ($Nvar != $1) { die "$Nvar in variants.tab does not have a corresponding line in domain_differences_AGADIR.tab\n" }
  $tangorefdomainscore = $2 ; $waltzrefdomainscore = $3;
  $tangomutdomainscore = $4 ; $waltzmutdomainscore = $5;

  # start the second pass throug the AGADIR output
  if ($tangorefpeptide eq '.' and $tangomutpeptide ne '.' and $Ntangoref > 0) {
    until ($N == $Nvar) {
      $N = 0 ;
      $line = <TANGOREF>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
    }
    while ($N == $Nvar) {
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 }
        $peptide =~ s/\t//g;
        $BEGIN = $tangomutposition;
        $END = $BEGIN + length($tangomutpeptide) - 1;
        if ($BEGIN <= $end and $END >= $begin) {
          $tangorefposition = $begin;
          $tangorefpeptide = $peptide;
          $tangorefregionscore = $score;
        }
        $line = <TANGOREF>;
      } else {
        $N = 0;
      }
    }
  }

  if ($tangorefpeptide ne '.' and $tangomutpeptide eq '.' and $Ntangomut > 0) {
    until ($N == $Nvar) {
      $N = 0 ;
      $line = <TANGOMUT>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
    }
    while ($N == $Nvar) {
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 }
        $peptide =~ s/\t//g;
        $BEGIN = $tangorefposition;
        $END = $BEGIN + length($tangorefpeptide) - 1;
        if ($BEGIN <= $end and $END >= $begin) {
          $tangomutposition = $begin;
          $tangomutpeptide = $peptide;
          $tangomutregionscore = $score;
        }
        $line = <TANGOMUT>;
      } else {
        $N = 0;
      }
    }
  }

  if ($waltzrefpeptide eq '.' and $waltzmutpeptide ne '.' and $Nwaltzref > 0) {
    until ($N == $Nvar) {
      $N = 0 ;
      $line = <WALTZREF>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
    }
    while ($N == $Nvar) {
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 }
        $peptide =~ s/\t//g;
        $BEGIN = $waltzmutposition;
        $END = $BEGIN + length($waltzmutpeptide) - 1;
        if ($BEGIN <= $end and $END >= $begin) {
          $waltzrefposition = $begin;
          $waltzrefpeptide = $peptide;
          $waltzrefregionscore = $score;
        }
        $line = <WALTZREF>;
      } else {
        $N = 0;
      }
    }
  }

  if ($waltzrefpeptide ne '.' and $waltzmutpeptide eq '.' and $Nwaltzmut > 0) {
    until ($N == $Nvar) {
      $N = 0 ;
      $line = <WALTZMUT>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
    }
    while ($N == $Nvar) {
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 }
        $peptide =~ s/\t//g;
        $BEGIN = $waltzrefposition;
        $END = $BEGIN + length($waltzrefpeptide) - 1;
        if ($BEGIN <= $end and $END >= $begin) {
          $waltzmutposition = $begin;
          $waltzmutpeptide = $peptide;
          $waltzmutregionscore = $score;
        }
        $line = <WALTZMUT>;
      } else {
        $N = 0;
      }
    }
  }

  # write the TANGO and WALTZ output, make sure that real numbers
  #   are written as decimals with 2 digits
  printf OUT '%s%.2f%s%.2f', "\t", $tangorefscore, "\t", $tangomutscore;
  printf OUT '%s%.2f', "\t", $tangomutscore - $tangorefscore;
  printf OUT '%s%.2f%s%.2f', "\t", $waltzrefscore, "\t", $waltzmutscore;
  printf OUT '%s%.2f', "\t", $waltzmutscore - $waltzrefscore;
  if ($tangorefdomainscore eq '.') {
    print OUT "\t.\t.\t.\t.";
  } else {
    printf OUT '%s%.2f%s%.2f', "\t", $tangorefdomainscore, "\t", $tangomutdomainscore;
    printf OUT '%s%.2f', "\t", $tangomutdomainscore - $tangorefdomainscore;
    printf OUT '%s%.2f', "\t", $tangomutdomainscore / $domain_length{$Nvar};
  }
  if ($waltzrefdomainscore eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    printf OUT '%s%.2f%s%.2f', "\t", $waltzrefdomainscore, "\t", $waltzmutdomainscore;
    printf OUT '%s%.2f', "\t", $waltzmutdomainscore - $waltzrefdomainscore;
  }
  if ($tangorefregionSscore eq '.') {
    print OUT "\t.";
    $TANGOREFREGIONSSCORE = 0;
      # replace nonexisting value by 0 to make computation difference possible
  } else {
    printf OUT '%s%.2f', "\t", $tangorefregionSscore;
    $TANGOREFREGIONSSCORE = $tangorefregionSscore;
  }
  if ($tangomutregionSscore eq '.') {
    print OUT "\t.";
    $TANGOMUTREGIONSSCORE = 0;
  } else {
    printf OUT '%s%.2f', "\t", $tangomutregionSscore;
    $TANGOMUTREGIONSSCORE = $tangomutregionSscore;
  }
  if ($tangorefregionSscore eq '.' and $tangomutregionSscore eq '.') {
    print OUT "\t.";
  } else {
    printf OUT '%s%.2f', "\t", $TANGOMUTREGIONSSCORE - $TANGOREFREGIONSSCORE;
  }
  if ($waltzrefregionSscore eq '.') {
    print OUT "\t.";
    $WALTZREFREGIONSSCORE = 0;
  } else {
    printf OUT '%s%.2f', "\t", $waltzrefregionSscore;
    $WALTZREFREGIONSSCORE = $waltzrefregionSscore;
  }
  if ($waltzmutregionSscore eq '.') {
    print OUT "\t.";
    $WALTZMUTREGIONSSCORE = 0;
  } else {
    printf OUT '%s%.2f', "\t", $waltzmutregionSscore;
    $WALTZMUTREGIONSSCORE = $waltzmutregionSscore;
  }
  if ($waltzrefregionSscore eq '.' and $waltzmutregionSscore eq '.') {
    print OUT "\t.";
  } else {
    printf OUT '%s%.2f', "\t", $WALTZMUTREGIONSSCORE - $WALTZREFREGIONSSCORE;
  }
  if ($tangorefpeptide eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    print OUT "\t$tangorefposition\t$tangorefpeptide";
    printf OUT '%s%.2f', "\t", $tangorefregionscore;
  }
  if ($tangomutpeptide eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    print OUT "\t$tangomutposition\t$tangomutpeptide";
    printf OUT '%s%.2f', "\t", $tangomutregionscore;
  }
  if ($tangorefpeptide eq '.' or $tangomutpeptide eq '.') {
    print OUT "\t.";
  } else {
    printf OUT '%s%.2f', "\t", $tangomutregionscore - $tangorefregionscore;
  }
  if ($tangomutMAXregionSposition eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    print OUT "\t$tangomutMAXregionSposition\t$tangomutMAXregionS";
    printf OUT '%s%.2f', "\t", $tangomutMAXregionSscore;
  }
  if ($waltzrefpeptide eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    print OUT "\t$waltzrefposition\t$waltzrefpeptide";
    printf OUT '%s%.2f', "\t", $waltzrefregionscore;
  }
  if ($waltzmutpeptide eq '.') {
    print OUT "\t.\t.\t.";
  } else {
    print OUT "\t$waltzmutposition\t$waltzmutpeptide";
    printf OUT '%s%.2f', "\t", $waltzmutregionscore;
  }
  if ($waltzrefpeptide eq '.' or $waltzmutpeptide eq '.') {
    print OUT "\t.";
  } else {
    printf OUT '%s%.2f', "\t", $waltzmutregionscore - $waltzrefregionscore;
  }

  # add SIFT and PROVEAN output
  if (exists $prediction{$predictionkey}) {
    print OUT "\t$prediction{$predictionkey}";
  } else {
    print OUT "\t.\t.\t.\t.\t.";
  }

  # eventually add standard UniProt/SwissProt sequence ID
  if ($SwissProtstandard) {
    $transcriptID =~ s/\.\d+$//; # remove version number
    if (exists $standard4refseqID{$transcriptID}) {
      print OUT "\t$standard4refseqID{$transcriptID}";
    } elsif (exists $standard4gene{$geneID}) {
      print OUT "\t~$standard4gene{$geneID}";
    } else {
      print OUT "\t.";
    }
  }

  print OUT "\n";
}
