#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# Its purpose is to compact the output of FoldX, so as to have no
# redundant information on different lines. It also adds FoldX output
# to the VCF file with the mutations selected for further analysis.
# It takes as input :
#   - FoldXreport.tab
#   - SNPAAchange.vcf
# It writes :
#   - FoldXreport_compact.tab
#   - withPDB.vcf

# parse FoldX report and make hash of arrays %fields
# each array contains the fields of the report and the hash key is unique,
# based on chromosome location + PDB file + mutation string
# If there are several items for the same key (because different transcripts)
# we transform the noncommon fields into a comma-separated list
open IN, 'FoldXreport.tab' or die "cannot open FoldXreport.tab\n";
# skip header line
<IN>;
# note : fields are :
# 0:var 1:chromloc 2:mut 3+4:gene 5:trID 6:trloc 7:PDB 8:mutst
while (<IN>) {
  @fields = split /\t/;
  $key = "$fields[1]\t$fields[7]\t$fields[8]\n";
  if (exists $fields{$key}) {
    $fields{$key}[0] .= ",$fields[0]";
    $fields{$key}[2] .= ",$fields[2]";
    $fields{$key}[3] .= ",$fields[3]";
    $fields{$key}[4] .= ",$fields[4]";
    $fields{$key}[5] .= ",$fields[5]";
    $fields{$key}[6] .= ",$fields[6]";
  } else {
    $fields{$key} = [ @fields ];
  }
}
close IN;

# make new hash %reportitems with the lines for the compact report and
# as key just the chromosome location.
# make also hash %VCFitems for the VCF INFO field
# note that if you have different aa mutations for the same site you
#   get multiple lines in one item because the last field ends with a newline
for $key (keys %fields) {
  @fields = @{ $fields{$key} };
  $reduced_key = $fields[1];
  $fields[2] = &unique_withslash($fields[2]);
  $fields[3] = &unique($fields[3]);
  $fields[4] = &unique($fields[4]);
  $reportitems{$reduced_key} .= join "\t", @fields;
  $VCFitems{$reduced_key} .= join "|", @fields[7..$#fields];
}

# make compact report and add FOLDX item to VCF file.
open VCFIN, 'SNPAAchange.vcf' or die "cannot open SNPAAchange.vcf\n";
open OUT, '>FoldXreport_compact.tab';
open VCFOUT, '>withPDB.vcf';
print OUT "variant(s)\tchromosomal location\tmutation\tgene ID\tgene name\ttranscript ID(s)\tlocation(s) on transcript\tPDB file\tmutationstring\tchains\tinteracting chains\tdelta interaction energy\tSD\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\tElectrostatics\tSolvation Polar\tSolvation Hydrophobic\tVan der Waals clashes\tentropy sidechain\tentropy mainchain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional clash\tbackbone clash\thelix dipole\twater bridge\tdisulfide\telectrostatic kon\tpartial covalent bonds\tenergy Ionisation\tEntropy Complex\n";
while (<VCFIN>) {
  if (/^#/) { # reading header
    print VCFOUT;
    if (/^##INFO=<ID=SNPVAR,Number/) {
      print VCFOUT "##INFO=<ID=FOLDX,Number=.,Type=String,Description=\"FoldX output: 'PDB file | mutationstring | chains | interacting chains | delta interaction energy | SD | total energy | Backbone Hbond | Sidechain Hbond | Van der Waals | Electrostatics | Solvation Polar | Solvation Hydrophobic | Van der Waals clashes | entropy sidechain | entropy mainchain | sloop_entropy | mloop_entropy | cis_bond | torsional clash | backbone clash | helix dipole | water bridge | disulfide | electrostatic kon | partial covalent bonds | energy Ionisation | Entropy Complex' \">\n";
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
      $FOLDXinfo = $VCFitems{$key};
      $FOLDXinfo =~ s/[,]/../g;
      $FOLDXinfo =~ s/[;]//g;
        # no , ; = allowed in VCF INFO items
      $FOLDXinfo =~ s/\n/,/g;
        # replace newlines by comma's to have all on one line
      $FOLDXinfo =~ s/,$//;
        # remove last comma
      $fields[7] = "$fields[7];FOLDX=$FOLDXinfo";
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
  # a set of unique slash-separeted items
  my ($set) = @_;
  my ($item, $uniqueset, @items, %items);
  @items = split /,/, $set;
  for $item (@items) {
    if (not exists $items{$item}) {
      $items{$item} = 1;
      $uniqueset .= "/$item";
    }
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
