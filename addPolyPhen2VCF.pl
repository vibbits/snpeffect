#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It extracts the Polyphen prediction from the file
# SEQANALreport_withPolyPhen.tab and adds it to the VCF file

%aa3 = ( 'A' => 'Ala', 'C' => 'Cys', 'D' => 'Asp', 'E' => 'Glu', 'F' => 'Phe',
  'G' => 'Gly', 'H' => 'His', 'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu',
  'M' => 'Met', 'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg',
  'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp', 'Y' => 'Tyr'  );

open IN, 'SEQANALreport_withPolyPhen.tab' or die "cannot open SEQANALreport_withPolyPhen.tab\n";
  # read the header line and determine which field contains the PolyPhen
  # prediction (can vary because info about standard SwissProt ID used
  # to restrict output might or might not be present).
@fields = split /\t/, <IN>;
$i = 0;
foreach $item (@fields) {
  if ($item eq 'PolyPhen prediction') {
    $I = $i;
  }
  $i++;
}
while (<IN>) {
  @fields = split /\t/;
  $key = "$fields[1]\t$fields[6]\t$fields[8]";
  if (exists $fields{$key}){
    $fields{$key}[3] .= ",$fields[55]"; # PolyPhen prediction
  } else {
    $fields{$key}[0] = $fields[1]; # chromosomal location
    $fields{$key}[1] = $fields[6]; # original aa
    $fields{$key}[2] = $fields[8]; # mutated aa
    $fields{$key}[3] = $fields[$I]; # PolyPhen prediction
  }
}
close IN;

for $key (keys %fields) {
  @fields = @{ $fields{$key} };
  $reduced_key = $fields[0];
  $fields[3] = &unique($fields[3]);
  $mutationstring = "$aa3{$fields[1]}2$aa3{$fields[2]}";
  $Polyphen_VCFitems{$reduced_key} .= "$mutationstring|$fields[3],";
}

open VCFIN, 'SNPpipelinereport.vcf' or die "cannot open SNPpipelinereport.vcf\n";
open VCFOUT, '>SNPpipelinereport_withPolyPhen.vcf';
while (<VCFIN>) {
  if (/^#/) { # reading header
    print VCFOUT;
    if (/^##INFO=<ID=PROVEAN,Number/) {
      print VCFOUT "##INFO=<ID=POLYPHEN,Number=.,Type=String,Description=\"PROVEAN output: ' mutation | PolyPhen prediction(s) ' \">\n";
    }
  } else {
    @fields = split /\t/;
    $key = "$fields[0]:$fields[1]-$fields[1]";
    if (exists $Polyphen_VCFitems{$key}) {
      $item = $Polyphen_VCFitems{$key}; $item =~ s/,$//;
      $fields[7] .= ";POLYPHEN=$item";
      $line = join "\t", @fields;
      print VCFOUT $line;
    } else {
      print VCFOUT;
    }
  }
}
close VCFIN ; close VCFOUT;

sub unique {
  # takes as input a set of comma-separated items and replaces it by
  # a set of unique slash-separeted items
  my ($set) = @_;
  my ($item, $uniqueset, @items, %items);
  @items = split /,/, $set;
  for $item (@items) {
    $items{$item} = 1;
  }
  for $item (sort keys %items) {
    $uniqueset .= "/$item";
  }
  $uniqueset =~ s/^\///;
  return $uniqueset;
}
