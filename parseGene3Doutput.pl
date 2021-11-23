#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the output of submit2Gene3D.pl and produces the file
# mutated_domains.tab with data to be incorparated in the final report.
# It also writes a log variants_without_domain_info.txt with info
# why for some variants no domain info could be provided.

# make hash of sequences that have no domain info
open IN, 'no_domain_info.txt'
  or die "cannot open no_domain_info.txt\n";
while (<IN>) {
  if (/^([\w\.]+) : (no Gene3D domains .*)$/) {
    $nodomaininfo{$1} = $2;
  } else {
    die "problem with no_domain_info.txt line :\n$_\n";
  }
}
close IN;

# make hash with entries in Gene3D output
open IN, 'seqs.domains'
  or die "cannot open seqs.domains\n";
<IN>; # skip first header
while (<IN>) {
  if (/^#domain_id,cath-superfamily/) {
    &treat_datablock;
    $datablock = '';
  } else {
    $datablock .= $_;
  }
}
close IN;
&treat_datablock; # last datablock

open IN, 'variants.tab' or die "cannot open variants.tab\n";
open OUT, '>mutated_domains.tab';
open LOG, '>variants_without_domain_info.txt';
while (<IN>) {
  if (not /^(\d+)\t(\w+:\d+-\d+)\t([A-Z]+\/[A-Z]+)\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
    # Nvar chrom_pos geneID genename seqID refAA seqpos altAA
    die "error in file variants.tab :\n$_\n";
  }
  $Nvar = $1;
  $seqID = $6;
  $seqpos = $8;

  if (exists $datablock{$seqID}) {
    $founddomain = 0;
    $domainID = '' ;
    $CATHfam = '' ; %foundCATHfam = ();
    @datalines = split /\n/, $datablock{$seqID};
    foreach $dataline (@datalines) {
      @fields = split /,/, $dataline;
        # domain ID, CATH family, Nvar, Gene3D ID, score, match range, .....
      $range = $fields[5];
      $range =~ /(\d+)-(\d+)/;
      $begin = $1 ; $end = $2;
      if ($seqpos >= $begin and $seqpos <= $end) {
        if ($founddomain) { # have already parsed line with same range
          $domainID .= ";$fields[0]";
          if (not $foundCATHfam{$fields[1]}) {
            $CATHfam .= ";$fields[1]" ; $foundCATHfam{$fields[1]} = 1;
          }
        } else {
          $founddomain = 1;
          $domainID = $fields[0];
          $CATHfam = $fields[1] ; $foundCATHfam{$fields[1]} = 1;
          $BEGIN = $begin ; $END = $end;
        }
      }
    }
    if ($founddomain) {
      print OUT "$Nvar\t$domainID\t$CATHfam\t$BEGIN\t$END\t$seqpos\n";
    } else {
      print LOG "var_$Nvar ($seqID) : no SNP in domains\n";
    }
  } elsif (exists $nodomaininfo{$seqID}) {
    print LOG "var_$Nvar ($seqID) : $nodomaininfo{$seqID}\n";
  } else {
    die "transcript $seqID cannot be found in seqs.domains or no_domain_info.txt\n";
  }
}
close IN ; close OUT ; close LOG;

########################################################

sub treat_datablock {
  if (not $datablock =~ /^[^,]+,[^,]+,([\w\.]+),/) {
    die "problem with data block :\n$datablock\n";
  } else {
    $datablock{$1} = $datablock;
  }
}
