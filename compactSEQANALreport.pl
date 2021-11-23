#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# The purpose is to compact the output of the various tools in the
# SNPpipeline that operate on the sequence (not the structure). It is
# however not possible to make a nonredundant set of the data in SEQANALreport
# in any sensible way. Therefore this script makes a selection of the data
# deemed most important (the user can always refer to the complete datasets).
# It also adds data to the VCF file.
# It takes as input :
#   - SEQANALreport.tab
#   - withPDB.vcf
# It writes :
#   - SEQANALreport_compact.tab
#   - SNPpipelinereport.vcf
# NOTE : if the pipeline is run with SwissProtstandard (only analyze
# for each gene the sequence that corresponds to the UniProt standard
# sequence) this compacting will only contribute little (only reduncant
# entries for a variant if there are overlapping genes or if UniProt
# has more than 1 standard because the alternative transcripts are too
# different). The script does detect automatically the presence of the
# UniProt standard in SEQANALreport.tab and includes it in the VCF output.

%aa3 = ( 'A' => 'Ala', 'C' => 'Cys', 'D' => 'Asp', 'E' => 'Glu', 'F' => 'Phe',
  'G' => 'Gly', 'H' => 'His', 'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu',
  'M' => 'Met', 'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg',
  'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp', 'Y' => 'Tyr' );

# parse SEQANAL report and make hash of array %fields
# each array contains a selected set of fields and the hash key is unique,
# based on chromosome location + original aa + mutant aa
# (several overlapping genes and their transripts can have the same key)
# We make a selection of fields for the final compact report
# For the "max delta" fields we search the extreme value
open IN, 'SEQANALreport.tab' or die "cannot open SEQANALreport.tab\n";
  # skip header line,
  #   but do test if there is a "UniProt standard sequence" field
$line = <IN>;
if ($line =~ /UniProt standard sequence/) { $SwissProtstandard = 1 }
while (<IN>) {
  chomp;
  @fields = split /\t/;
  $key = "$fields[1]\t$fields[6]\t$fields[8]";
  if (exists $fields{$key}){
    $fields{$key}[0] .= ",$fields[0]"; # variant
    $fields{$key}[2] .= ",$fields[2]"; # mutation
    $fields{$key}[3] .= ",$fields[3]"; # gene ID
    $fields{$key}[4] .= ",$fields[4]"; # gene name
    $fields{$key}[5] .= ",$fields[5]"; # transcript ID
    $fields{$key}[7] .= ",$fields[7]"; # position mutaton in transcript
    if ($fields[9] ne '.') {
      $fields{$key}[9] .= ",$fields[9]"; # CATH domain
    } # note that you can have . or ;-separated list
    $fields{$key}[10] .= ",$fields[15]"; # mutation in TM?
    $fields{$key}[11] = &extreme($fields{$key}[11], $fields[18]);
      # max delta TANGO score
    $fields{$key}[12] = &extreme($fields{$key}[12], $fields[21]);
      # max delta WALTZ score
    $fields{$key}[13] = &extreme($fields{$key}[13], $fields[24]);
      # max delta TANGO score domain
 
    $fields{$key}[14] = &extreme($fields{$key}[14], $fields[25]);
      # max normalized TANGO score domain

    $fields{$key}[15] = &extreme($fields{$key}[15], $fields[28]);
      # max delta WALTZ score domain
    $fields{$key}[16] = &extreme($fields{$key}[16], $fields[31]);
      # max delta TANGO score regions in domain
    $fields{$key}[17] = &extreme($fields{$key}[17], $fields[34]);
      # max delta WALTZ score regions in domain
    $fields{$key}[18] = &extreme($fields{$key}[18], $fields[41]);
      # max delta TANGO score region
    $fields{$key}[19] = &extreme($fields{$key}[19], $fields[51]);
      # max delta WALTZ score region
    $fields{$key}[20] .= ",$fields[52]"; # SIFT prediction
    $fields{$key}[21] .= ",$fields[56]"; # PROVEAN prediction
    if ($SwissProtstandard) {
      $fields{$key}[22] .= ",$fields[58]";
    }
  } else {
    $fields{$key}[0] = $fields[0]; # variant
    $fields{$key}[1] = $fields[1]; # chromosomal location
    $fields{$key}[2] = $fields[2]; # mutation
    $fields{$key}[3] = $fields[3]; # gene ID
    $fields{$key}[4] = $fields[4]; # gene name
    $fields{$key}[5] = $fields[5]; # transcript ID
    $fields{$key}[6] = $fields[6]; # original aa
    $fields{$key}[7] = $fields[7]; # position mutaton in transcript
    $fields{$key}[8] = $fields[8]; # mutated aa
    $fields{$key}[9] = $fields[9]; # CATH domain
    $fields{$key}[10] = $fields[15]; # mutation in TM?
    $fields{$key}[11] = $fields[18]; # max delta TANGO score
    $fields{$key}[12] = $fields[21]; # max delta WALTZ score
    $fields{$key}[13] = $fields[24]; # max delta TANGO score domain

    $fields{$key}[14] = $fields[25]; # max normalized TANGO score domain

    $fields{$key}[15] = $fields[28]; # max delta WALTZ score domain
    $fields{$key}[16] = $fields[31]; # max delta TANGO score regions in domain
    $fields{$key}[17] = $fields[34]; # max delta WALTZ score regions in domain
    $fields{$key}[18] = $fields[41]; # max delta TANGO score region
    $fields{$key}[19] = $fields[51]; # max delta WALTZ score region
    $fields{$key}[20] = $fields[52]; # SIFT prediction
    $fields{$key}[21] = $fields[56]; # PROVEAN prediction
    if ($SwissProtstandard) {
      $fields{$key}[22] = $fields[58];
    }
  }
}

# make new hash %reportitems with the lines for the compact report and
# as key just the chromosome location.
# Make also hashes for the VCF INFO field
for $key (keys %fields) {
  @fields = @{ $fields{$key} };
  $reduced_key = $fields[1];
  $fields[2] = &unique_withslash($fields[2]);
  $fields[3] = &unique($fields[3]);
    # combine gene ID's in nonredundant set
  $fields[4] = &unique($fields[4]);
  $fields[9] =~ s/;/,/g;
    # remove ; from CATH ID field, ; accepted in VCF format
  $fields[9] = &unique($fields[9]);
  $fields[9] =~ s/^\.\///;
  $fields[10] = &unique($fields[10]);
  $fields[20] = &unique($fields[20]);
  $fields[21] = &unique($fields[21]);
  if ($SwissProtstandard) {
    $fields[22] = &unique($fields[22]);
  }
  $reportitems{$reduced_key} .= join "\t", @fields;
  $reportitems{$reduced_key} .= "\n";

  $mutationstring = "$aa3{$fields[6]}2$aa3{$fields[8]}";
  if ($fields[9] ne '.') {
    $CATH_VCFitems{$reduced_key} .= "$fields[9],";
  }
  $TMHMM_VCFitems{$reduced_key} .= "$fields[10],";
  $AGADIR_VCFitems{$reduced_key} .= "$mutationstring|$fields[11]|$fields[12]|$fields[13]|$fields[14]|$fields[15]|$fields[16]|$fields[17]|$fields[18]|$fields[19],";
  $SIFT_VCFitems{$reduced_key} .= "$mutationstring|$fields[20],";
  $PROVEAN_VCFitems{$reduced_key} .= "$mutationstring|$fields[21],";
  if ($SwissProtstandard) {
    $SwissProtstandard_VCFitems{$reduced_key} .= "$fields[22],";
  }
}

# make compact report and add extra items to VCF file
open VCFIN, 'withPDB.vcf' or die "cannot open withPDB.vcf\n";
open OUT, '>SEQANALreport_compact.tab';
open VCFOUT, '>SNPpipelinereport.vcf';
print OUT "variant(s)\tchromosomal location\tmutation\tgene ID(s)\tgene name(s)\ttranscript ID(s)\treference aa\tposition(s) in transcript\tmutant aa\tCATH domain ID(s)\tmutation in TM\t";
print OUT "max delta TANGO score\tmax delta WALTZ score\tmax delta TANGO score domain\tmax normalized TANGO score domain\tmax delta WALTZ score domain\tmax delta TANGO score regions in domain\tmax delta WALTZ score regions in domain\tmax delta TANGO score region\tmax delta WALTZ score region\t";
print OUT "SIFT prediction(s)\tPROVEAN prediction(s)";
if ($SwissProtstandard) {
  print OUT "\tUniProt standard sequence(s)";
}
print OUT "\n";
while (<VCFIN>) {
  if (/^#/) { # reading header
    print VCFOUT;
    if (/^##INFO=<ID=FOLDX,Number/) {
      print VCFOUT "##INFO=<ID=CATH_GENE3D,Number=.,Type=String,Description=\"Gene3D output: ' CATH domain ID(s) ' \">\n";
      print VCFOUT "##INFO=<ID=TMHMM,Number=1,Type=String,Description=\"TMHMM output: ' mutation in TM ' \">\n";
      print VCFOUT "##INFO=<ID=TANGO_WALTZ,Number=.,Type=String,Description=\"AGADIR TANGO and WALTZ output: ' mutation | max delta TANGO score | max delta WALTZ score | max delta TANGO score domain | max normalized TANGO score domain | max delta WALTZ score domain | max delta TANGO score regions in domain | max delta WALTZ score regions in domain | max delta TANGO score region | max delta WALTZ score region ' \">\n";
      print VCFOUT "##INFO=<ID=SIFT,Number=.,Type=String,Description=\"SIFT output: ' mutation | SIFT prediction(s) ' \">\n";
      print VCFOUT "##INFO=<ID=PROVEAN,Number=.,Type=String,Description=\"PROVEAN output: ' mutation | PROVEAN prediction(s) ' \">\n";
      if ($SwissProtstandard) {
        print VCFOUT "##INFO=<ID=UNIPROT,Number=.,Type=String,Description=\"UniProt standard sequence\">\n";
      }
    }
  } else {
    # note that the header does not contain tabs and hence is not modified
    @fields = split /\t/;
      # fields are CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ...
    $fields[0] =~ s/^chr//;
      # chromosome name in out.vcf can start with "chr", name in out.fa does not
    $key = "$fields[0]:$fields[1]-$fields[1]";
    if (exists $reportitems{$key}) {
      print OUT $reportitems{$key};

      if (exists $CATH_VCFitems{$key}) {
        $item = $CATH_VCFitems{$key};
        $item =~ s/\//,/g; $item = &unique($item);
        $fields[7] .= ";CATH_GENE3D=$item";
      }
      $item = $TMHMM_VCFitems{$key};
      $item =~ s/\//,/g ; $item = &unique($item);
      $fields[7] .= ";TMHMM=$item";
      $item = $AGADIR_VCFitems{$key} ; $item =~ s/,$//;
      $fields[7] .= ";TANGO_WALTZ=$item";
      $item = $SIFT_VCFitems{$key} ; $item =~ s/,$//;
      $fields[7] .= ";SIFT=$item";
      $item = $PROVEAN_VCFitems{$key} ; $item =~ s/,$//;
      $fields[7] .= ";PROVEAN=$item";
      if ($SwissProtstandard) {
        $item = $SwissProtstandard_VCFitems{$key} ; $item =~ s/,$//;
        if ($item ne '.') {
          $fields[7] .= ";UNIPROT=$item";
        }
      }
      $line = join "\t", @fields;
      print VCFOUT $line;
    } else {
      print VCFOUT;
    }
  }
}
close OUT, VCFIN, VCFOUT;

sub unique {
  # takes as input a set of comma-separated items and replaces it by
  # a set of unique slash-separated items
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

sub unique_withslash {
  # takes as input a set of comma-separated items and replaces it by
  # a set of unique slash-separated items, but takes account that
  # the items have a slash themselves
  my ($set) = @_;
  my ($item, $uniqueset, @items, %items);
  @items = split /,/, $set;
  for $item (@items) {
    $items{$item} = 1;
  }
  for $item (sort keys %items) {
    $uniqueset .= "//$item";
  }
  $uniqueset =~ s/^\/\///;
  return $uniqueset;
}

sub extreme {
  # returns which of both has highest absolute value
  my ($firstvalue, $secondvalue) = @_;
  if ($firstvalue eq '.' or $secondvalue ne '.' and abs($secondvalue) > abs($firstvalue)) {
    return  $secondvalue;
  } else {
    return $firstvalue;
  }
}
