#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It adds the PolyPhen predictions to the final report and uses
# the existence of a PolyPhen prediction to select from the final
# report as well as from the FoldX predictions a representative
# transcript for each gene

$correspondencefile = $ARGV[0];
#$correspondencefile = '../scripts/NP_NM.tab';
if (not $correspondencefile) {
  die "you must give the file NP_NM.tab as first argument\n";
}

# make a hash (%RefSeqmrnaID) with RefSeq protein ID as key and corresponding
# RefSeq mRNA ID as value
open IN, $correspondencefile or die "cannot open $correspondencefile\n";
while (<IN>) {
  $_ =~ /(NP_.+)\t(NM_.+)$/;
  $RefSeqmrnaID{$1} = $2;
}
close IN;

# make a hash (%RefSeqmrnaID) with RefSeq mRNA ID as key and the
# PolyPhen.predictions fields you want to write in the final report as value
# Also parse SNP4PolyPhen.list to find which variants are not in
# PolyPhen.predictions because PolyPhen mapsnps.pl failed to find them
open PRED, 'PolyPhen.predictions.tab'
  or die "cannot open PolyPhen.predictions.tab\n";
open LIST, 'SNP4PolyPhen.list' or die "cannot open SNP4PolyPhen.list\n";
open LOG, '>variants_without_PolyPhen_standard.txt';
<PRED>; # skip header
while (<PRED>) {
  chomp;
  s/ +//g; # remove blanks, leave only tabs
  @fields = split /\t/;
  $fields[11] =~ s/damaging/ damaging/; # restore lost space
  $IDstring = $fields[55];
  $IDstring =~ s/^#//;
  @IDitems = split /\|/, $IDstring;
  $chrompos = $IDitems[0];
  $mutation = $IDitems[1];
  $mutation =~ s/(.)(.)/$1\/$2/;
  until ($line eq "$chrompos\t$mutation") {
    $line = <LIST>;
    chomp $line;
    if ($line ne "$chrompos\t$mutation") {
      $line =~ s/\t/ /;
      print LOG "$line : PolyPhen found no UniProt protein\n";
    }
    if (eof LIST and $line ne "$chrompos\t$mutation") {
      die "cannot find $chrompos $mutation in SNP4PolyPhen.list\n";
    }
  }
  if (exists $IDitems[4]) {
    $RefSeqproteinID = $IDitems[4];
    if (exists $RefSeqmrnaID{$RefSeqproteinID}) {
      $key = "$chrompos\t$RefSeqmrnaID{$RefSeqproteinID}\t$fields[2]\t$fields[1]\t$fields[3]";
        # key is like : chr10:101715198 NM_015221 P 678 H
      $value = "$fields[11]\t$fields[12]\t$fields[14]\t$fields[15]\t$fields[16]\t$fields[17]\t$fields[18]\t$fields[O]\t$fields[1]";
      $PolyPhenfields{$key} = $value;
    } else {
      print LOG "$chrompos $mutation $RefSeqproteinID : no corresponding mRNA in RefSeq\n";
    }
  } else {
    print LOG "$chrompos $mutation $fields[0] : PolyPhen found no RefSeq protein \n";
  }
}
close LIST ; close PRED ; close LOG;

# parse SEQANALreport.tab and make SEQANALreport_withPolyPhen.tab, adding
# PolyPhen fields bases on commom RefSeq mRNA ID
# make also a hash (%variants_with_standard) with the numbers of the
# selected variants
open IN, 'SEQANALreport.tab' or die "cannot open SEQANALreport.tab\n";
open OUT, '>SEQANALreport_withPolyPhen.tab';
$line = <IN>; # read header
chomp $line;
print OUT $line;
print OUT "\tPolyPhen prediction\tbased on\tPolyPhen-2 class\tPolyPhen-2 prob\tPolyPhen-2 FPR\tPolyPhen-2 TPR\tPolyPhen-2 FDR\tUniProt ID\tUniProt position\n";
while (<IN>) {
  chomp;
  if (not /^var_\d+\t(\w+:\d+)-\d+\t/) {
    die "error in SEQANALreport.tab :\n$_\n";
  }
  $chrompos = "chr$1"; # Polyphen uses chromosome names that start with chr
  @fields = split;
  $Nvar = $fields[0];
  $transcriptID = $fields[5];
  $transcriptID =~ s/\.\d+$//; # PolyPhen does not have the version extension
  $refaa = $fields[6];
  $pos = $fields[7];
  $mutaa = $fields[8];
  $key = "$chrompos\t$transcriptID\t$refaa\t$pos\t$mutaa";
  if (exists $PolyPhenfields{$key}) {
    print OUT "$_\t$PolyPhenfields{$key}\n";
    $variants_with_standard{$Nvar} = 1;
  }
}
close IN ; close OUT;

# parse FoldXreport.tab and select variants_with_standard
# NOTE that we only need this to count the genes that have both a structure
# in the Switchlab database and a standard in the PolyPhen database, the
# output will be deleted later.
open IN, 'FoldXreport.tab' or die "cannot open FoldXreport.tab\n";
open OUT, '>FoldXreport_withPolyPhenonly.tab';
$line = <IN> ; print OUT $line; # copy header
while (<IN>) {
  $_ =~ /^(var_\d+)\t/;
  if ($variants_with_standard{$1}) {
    print OUT;
  }
}
close IN ; close OUT;
