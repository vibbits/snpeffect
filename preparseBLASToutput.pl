#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It preprocesses the BLAST output before it is given to parseBLASToutput,pl
# It does 2 things :
#   - add info about variants that were not in the nonredundant set
#     of reference sequences
#   - remove lines of text that are not used by parseBLASToutput.pl,
#     making the file much smaller, even if it contains more entries
#  It takes as input :
#    - sequence_identities.tab
#    - out.fmt7_blastp
#  It makes as output BLASTfile_preprocessed

$BLASTfile = 'out.fmt7_blastp';

# make a hash with the identical sequence correspondances
open IN, 'sequence_identities.tab'
  or die "cannot open sequence_identities.tab\n";
while (<IN>) {
  if (/^(\d+)\t(.+)$/) {
    @correspondingIDs = split /\t/, $2;
    $correspondingIDs{$1} = [@correspondingIDs];
  }
}
close IN;

open IN, "$BLASTfile" or die "cannot open $BLASTfile\n";
open OUT, '>BLASTfile_preprocessed';
$firstentry = 1;
while (<IN>) {
  if (/^# BLASTP/) {
    if ($firstentry) { $firstentry = 0 } else { &treatentry }
    $entry = '';
    $entry = $_;
  } elsif (/^# Query: (\d+)/) {
    $Nvar_entry = $1;
    $entry .= $_;
  } elsif (/^(\d+)\t(\d+)_(\w)\t([0-9\.]+\t(.+)\n$)/) {
      # variant number ; PDB number ; chain ; %ident ; rest of line
    if ($1 != $Nvar_entry) { die "var number mismatch at $_ \n"; }
    if ($4 == 100.0) {
      $entry .= $_;
    }
  }
}
&treatentry;
close IN ; close OUT;

#######################################################################

sub treatentry {
  print OUT $entry;
  if (exists $correspondingIDs{$Nvar_entry}) {
    foreach $n (@{$correspondingIDs{$Nvar_entry}}) {
      $copy_of_entry = $entry;
      $copy_of_entry =~ s/# Query: $Nvar_entry/# Query: $n/;
      $copy_of_entry =~ s/\n$Nvar_entry/\n$n/g;
      print OUT $copy_of_entry;
    }
  }
}
