#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of AGADIR. It takes a input :
# - the 6 files written by submit2AGADIR.pl
# - the file variants.tab with information about the mutation
# It produces as output region_differences.AGADIR.tab with info about the
# effect  of the mutation on the TANGO and WALTZ predictions of regions prone
# to aggegation or amyloid formation. Note that it makes one line for each
# variant. If the mutation is not present in a predicted region the output
# file will contain only information about the global TANGO or WALTZ score.
# A special case is when the mutation causes a region to shift or shrink
# rather then to appear or disappear so that the mutation comes to lie
# outside ; to handle these cases the script makefinalreport.pl will
# do a second pass to still add the missing region.

open VAR, 'variants.tab' or die "cannot open variants.tab\n";
open AGAREF, 'reference_AGADIR_summary.tab'
  or die "cannot open reference_AGADIR_summary.tab\n";
open TANGOREF, 'reference_TANGO.tab' or die "cannot open reference_TANGO.tab\n";
open WALTZREF, 'reference_WALTZ.tab' or die "cannot open reference_WALTZ.tab\n";
open AGAMUT, 'mutated_AGADIR_summary.tab'
  or die "cannot open mutated_AGADIR_summary.tab\n";
open TANGOMUT, 'mutated_TANGO.tab' or die "cannot open mutated_TANGO.tab\n";
open WALTZMUT, 'mutated_WALTZ.tab' or die "cannot open mutated_WALTZ.tab\n";
open OUT, '>region_differences_AGADIR.tab';

while (<VAR>) {
  if (not /^(\d+)\t(\w+:\d+-\d+)\t([A-Z]+\/[A-Z]+)\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
    die "error in file variants.tab :\n$_\n";
  }
  $Nvar = $1;
  $refAA = $7;
  $POS = $8;
  $mutAA = $9;

  # move down in reference_AGADIR.tab
  <AGAREF> ; <AGAREF> ;
  $line = <AGAREF>;
  if ($line !~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\n/) { die "error in reference_AGADIR_summary.tab line :\n$line" }
  if ($1 != $Nvar) { die "$Nvar in variants.tab does not have a corresponding line in reference_AGADIR_summary.tab\n" }
  $Ntangoref = $2 ; $Nwaltzref = $3 ;
  $tangorefscore = $5 ; $waltzrefscore = $6;

  # move down in mutated_AGADIR.tab
  <AGAMUT> ; <AGAMUT> ;
  $line = <AGAMUT>;
  if ($line !~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\n/) { die "error in mutated_AGADIR_summary.tab line :\n$line" }
  if ($1 != $Nvar) { die "$Nvar in variants.tab does not have a corresponding line in mutated_AGADIR_summary.tab\n" }
  $Ntangomut = $2 ; $Nwaltzmut = $3 ;
  $tangomutscore = $5 ; $waltzmutscore = $6;

  # handle TANGO output for the reference sequence
  # NOTE : it might be better to use a subroutine instead of copying
  #   this code 3 times, but here we sacrifice program short code
  #   for avoiding complexity
  $tangorefposition = '.'; $tangorefpeptide = '.'; $tangorefregionscore = '.';
  if ($Ntangoref > 0) { # 0 TANGO windows means no entry in the file
    until ($N == $Nvar) { # move down in reference_TANGO.tab
      $N = 0 ;
      $line = <TANGOREF>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
      if (eof TANGOREF and $N != $Nvar) {
        die "cannot find variant ${N} while parsing through reference_TANGO.tab\n";
      }
    }
    while ($N == $Nvar) { # treat entry, search region that has mutation
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
          # peptide is region with preceding and following aa
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 } # in case region begins with first aa
        $peptide =~ s/\t//g;
        if ($POS >= $begin and $POS <= $end) {
          if (substr($peptide, $POS - $begin, 1) ne $refAA) {
            die "mutation in $Nvar $peptide not found\n";
          }
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

  # handle TANGO output for the mutated sequence
  $tangomutposition = '.'; $tangomutpeptide = '.'; $tangomutregionscore = '.';
  if ($Ntangomut > 0) { # 0 TANGO windows means no entry in the file
    until ($N == $Nvar) { # move down in mutated_TANGO.tab
      $N = 0 ;
      $line = <TANGOMUT>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
      if (eof TANGOMUT and $N != $Nvar) {
        die "cannot find variant ${N} while parsing through mutated_TANGO.tab\n";
      }
    }
    while ($N == $Nvar) { # treat entry, search region that has mutation
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
          # peptide is region with preceding and following aa
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 } # in case region begins with first aa
        $peptide =~ s/\t//g;
        if ($POS >= $begin and $POS <= $end) {
          if (substr($peptide, $POS - $begin, 1) ne $mutAA) {
            die "mutation in $Nvar $peptide not found\n";
          }
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

  # handle WALTZ output for the reference sequence
  $waltzrefposition = '.'; $waltzrefpeptide = '.'; $waltzrefregionscore = '.';
  if ($Nwaltzref > 0) { # 0 WALTZ windows means no entry in the file
    until ($N == $Nvar) { # move down in reference_WALTZ.tab
      $N = 0 ;
      $line = <WALTZREF>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
      if (eof WALTZREF and $N != $Nvar) {
        die "cannot find variants ${N} while parsing through reference_WALTZ.tab\n";
      }
    }
    while ($N == $Nvar) { # treat entry, search region that has mutation
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
          # peptide is region with preceding and following aa
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 } # in case region begins with first aa
        $peptide =~ s/\t//g;
        if ($POS >= $begin and $POS <= $end) {
          if (substr($peptide, $POS - $begin, 1) ne $refAA) {
            die "mutation in $Nvar $peptide not found\n";
          }
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

  # handle WALTZ output for the mutated sequence
  $waltzmutposition = '.'; $waltzmutpeptide = '.'; $waltzmutregionscore = '.';
  if ($Nwaltzmut > 0) { # 0 WALTZ windows means no entry in the file
    until ($N == $Nvar) { # move down in mutated_WALTZ.tab
      $N = 0 ;
      $line = <WALTZMUT>;
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1;
      }
      if (eof WALTZMUT and $N != $Nvar) {
        die "cannot find variant ${N} while parsing through mutated_WALTZ.tab\n";
      }
    }
    while ($N == $Nvar) { # treat entry, search region that has mutation
      if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
        $N = $1 ; $peptide = $3 ; $score = $4;
          # peptide is region with preceding and following aa
        $begin = $2 ; $end = $begin + $5;
        if ($begin == 0) { $begin = 1 } # in case region begins with first aa
        $peptide =~ s/\t//g;
        if ($POS >= $begin and $POS <= $end) {
          if (substr($peptide, $POS - $begin, 1) ne $mutAA) {
            die "mutation in $Nvar $peptide not found\n";
          }
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

  print OUT "$Nvar\t$tangorefscore\t$tangomutscore\t";
  print OUT "$tangorefposition\t$tangorefpeptide\t$tangorefregionscore\t";
  print OUT "$tangomutposition\t$tangomutpeptide\t$tangomutregionscore\t";
  print OUT "$waltzrefscore\t$waltzmutscore\t";
  print OUT "$waltzrefposition\t$waltzrefpeptide\t$waltzrefregionscore\t";
  print OUT "$waltzmutposition\t$waltzmutpeptide\t$waltzmutregionscore\t";
  print OUT "$Ntangoref\t$Ntangomut\t$Nwaltzref\t$Nwaltzmut\n";
    # this last line is only useful for the second parsing, performed
    # to add info about missing regions, and will not appear in final output
}
