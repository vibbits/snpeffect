#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the SnpEff output derived file mutated_sequences.fa and
# extracts data (location on genome of particular nucleotide change)
# that can serve as input for PolyPhen. It will be necessary to remove
# redundancy from the output.

open MUT, 'mutated_sequences.fa' or die "cannot open mutated_sequences.fa\n";
open OUT, '>SNP4PolyPhen_unsorted.list';
while (<MUT>) {
  if (not /^>\d+\t[\w\.]+ Variant (\w+:\d+)-\d+ Ref:([A-Z]) Alt:([A-Z]) HGVS/) {
    die "error in mutated_sequences.fa :\n$_\n";
  }
  print OUT "chr$1\t$2/$3\n";
    # PolyPhen wants chromosome names starting with chr
  <MUT>; # skip sequence
}
