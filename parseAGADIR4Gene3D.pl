#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of AGADIR. It takes a input :
# - the 6 files written by submit2AGADIR.pl
# - the file mutated_domains.tab written by finddomains.pl, which contains
#   information about the protein domain that contains the mutation.
# It produces as output domain_regions_differences_AGADIR.tab with info
# about the effect of the mutation on the TANGO and WALTZ predictions of
# regions prone to aggegation or amyloid formation, summed for all the regions
# that are inside a domain predicted by Gene3D. Note that it makes one line
# for each variant. If the mutation is not present in a predicted region the
# output file will contain only information about the global TANGO or WALTZ
# score.

# read in information about domain that contains mutation
open IN, 'mutated_domains.tab' or die "cannot open mutated_domains.tab\n";
while (<IN>) {
  if (not /^(\d+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)$/) {
    die "error in file mutated_domains.tab :\n$_\n";
  } else { # $1 is Nvar
    $domainID{$1} = $2;
    $familyID{$1} = $3;
    $begindomain{$1} = $4;
    $enddomain{$1} = $5;
  }
}

open AGAREF, 'reference_AGADIR_summary.tab'
  or die "cannot open reference_AGADIR_summary.tab\n";
open TANGOREF, 'reference_TANGO.tab' or die "cannot open reference_TANGO.tab\n";
open WALTZREF, 'reference_WALTZ.tab' or die "cannot open reference_WALTZ.tab\n";
open AGAMUT, 'mutated_AGADIR_summary.tab'
  or die "cannot open mutated_AGADIR_summary.tab\n";
open TANGOMUT, 'mutated_TANGO.tab' or die "cannot open mutated_TANGO.tab\n";
open WALTZMUT, 'mutated_WALTZ.tab' or die "cannot open mutated_WALTZ.tab\n";
open OUT, '>domain_regions_differences_AGADIR.tab';

while (<AGAREF>) {
  if (not /^(TANGO|name)/) {
    if (not /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\n/) { die "error in reference_AGADIR.tab line $_" }
    $Nvar = $1 ; $Ntangoref = $2 ; $Nwaltzref = $3 ;
    $tangorefscore = $5 ; $waltzrefscore = $6;
    $SNP_IN_DOMAIN = 0;
    if (exists $begindomain{$Nvar}) {
      $SNP_IN_DOMAIN = 1;
      $BEGIN = $begindomain{$Nvar};
      $END = $enddomain{$Nvar};
    }

    # move down in mutated_AGADIR.tab
    <AGAMUT> ; <AGAMUT> ;
    $line = <AGAMUT>;
    if ($line !~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t([-e\d\.]+)\t([-e\d\.]+)\t([-e\d\.]+)\n/) { die "error in mutated_AGADIR.tab line $_" }
    if ($1 != $Nvar) { die "$Nvar in reference_AGADIR.tab does not have a corresponding line in mutated_AGADIR.tab\n" }
    $Ntangomut = $2 ; $Nwaltzmut = $3 ;
    $tangomutscore = $5 ; $waltzmutscore = $6;

    # handle TANGO output for the reference sequence
    # NOTE : it might be better to use a subroutine instead of copying
    #   this code 3 times, but here we sacrifice program short code
    #   for avoiding complexity
    $tangorefregionSscore = '.';
    if ($SNP_IN_DOMAIN and $Ntangoref > 0) {
      # 0 TANGO windows means no entry in the file.
      # We only treat cases where the mutation is in a protein domain and
      # there are predicted aggregration-prone regions.
      $tangorefregionSscore = 0;
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
      while ($N == $Nvar) { # treat entry, search regions inside domain
        if ($line =~ /^(\d+)\t(\d+)\t[A-Z]?\t[A-Z]+\t[A-Z]?\t([\d\.]+)\t(\d+)$/) {
          $N = $1 ; $score = $3;
          $begin = $2 ; $end = $2 + $4;
          if ($begin == 0) { $begin = 1 }
          if ($begin >= $BEGIN and $end <= $END) {
            $tangorefregionSscore += $score;
          }
          $line = <TANGOREF>;
        } else {
          $N = 0;
        }
      }
      if ($tangorefregionSscore == 0) { $tangorefregionSscore = '.' }
    }

    # handle TANGO output for the mutated sequence
    # also search the region with the highest score
    $tangomutregionSscore = '.';
    $tangomutMAXregioninfo = ".\t.\t.";
    if ($SNP_IN_DOMAIN and $Ntangomut > 0) {
      # 0 TANGO windows means no entry in the file.
      # We only treat cases where the mutation is in a protein domain and
      # there are predicted aggragration-prone regions.
      $tangomutregionSscore = 0;
      $tangomutMAXregionSscore = 0;
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
      while ($N == $Nvar) { # treat entry, search regions inside domain
        if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
          $N = $1 ; $peptide = $3 ; $score = $4;
          $begin = $2 ; $end = $begin + $5;
          if ($begin == 0) { $begin = 1 }
          $peptide =~ s/\t//g;
          if ($begin >= $BEGIN and $end <= $END) {
            $tangomutregionSscore += $score;
            if ($score > $tangomutMAXregionSscore) {
              $tangomutMAXregionSscore = $score;
              $tangomutMAXregioninfo = "$begin\t$peptide\t$score";
            }
          }
          $line = <TANGOMUT>;
        } else {
          $N = 0;
        }
      }
      if ($tangomutregionSscore == 0) {
        $tangomutregionSscore = '.';
        $tangomutMAXregioninfo = ".\t.\t.";
      }
    }

    # handle WALTZ output for the reference sequence
    $waltzrefregionSscore = '.';
    if ($SNP_IN_DOMAIN and $Nwaltzref > 0) {
      # 0 WALTZ windows means no entry in the file.
      # We only treat cases where the mutation is in a protein domain and
      # there are predicted aggragration-prone regions.
      $waltzrefregionSscore = 0;
      until ($N == $Nvar) { # move down in reference_WALTZ.tab
        $N = 0 ;
        $line = <WALTZREF>;
        if ($line =~ /^(\d+)\t(\d+)\t([A-Z]?\t[A-Z]+\t[A-Z]?)\t([\d\.]+)\t(\d+)$/) {
          $N = $1;
        }
        if (eof WALTZREF and $N != $Nvar) {
          die "cannot find variant ${N} while parsing through reference_WALTZ.tab\n";
        }
      }
      while ($N == $Nvar) { # treat entry, search regions inside domain
        if ($line =~ /^(\d+)\t(\d+)\t[A-Z]?\t[A-Z]+\t[A-Z]?\t([\d\.]+)\t(\d+)$/) {
          $N = $1 ; $score = $3;
          $begin = $2 ; $end = $2 + $4;
          if ($begin == 0) { $begin = 1 }
          if ($begin >= $BEGIN and $end <= $END) {
            $waltzrefregionSscore += $score;
          }
          $line = <WALTZREF>;
        } else {
          $N = 0;
        }
      }
      if ($waltzrefregionSscore == 0) { $waltzrefregionSscore = '.' }
    }

    # handle WALTZ output for the mutated sequence
    $waltzmutregionSscore = '.';
    if ($SNP_IN_DOMAIN and $Nwaltzmut > 0) {
      # 0 WALTZ windows means no entry in the file.
      # We only treat cases where the mutation is in a protein domain and
      # there are predicted aggragration-prone regions.
      $waltzmutregionSscore = 0;
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
      while ($N == $Nvar) { # treat entry, search regions inside domain
        if ($line =~ /^(\d+)\t(\d+)\t[A-Z]?\t[A-Z]+\t[A-Z]?\t([\d\.]+)\t(\d+)$/) {
          $N = $1 ; $score = $3;
          $begin = $2 ; $end = $2 + $4;
          if ($begin == 0) { $begin = 1 }
          if ($begin >= $BEGIN and $end <= $END) {
            $waltzmutregionSscore += $score;
          }
          $line = <WALTZMUT>;
        } else {
          $N = 0;
        }
      }
      if ($waltzmutregionSscore == 0) { $waltzmutregionSscore = '.' }
    }

    print OUT "$Nvar\t$tangorefscore\t$tangomutscore\t";
    print OUT "$waltzrefscore\t$waltzmutscore\t";
    if (exists $domainID{$Nvar}) {
      print OUT "$domainID{$Nvar}\t$familyID{$Nvar}\t";
    } else {
      print OUT ".\t.\t";
    }
    print OUT "$tangorefregionSscore\t$tangomutregionSscore\t";
    print OUT "$waltzrefregionSscore\t$waltzmutregionSscore\t";
    print OUT "$tangomutMAXregioninfo\n";
  }
}
