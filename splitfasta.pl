#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It splits a fastA file in a number of equal parts (except for the last one)
# Its purpose is to be able to perform a BLAST search in parts

$fastafile = 'reference_sequences_nonredundant.fa';
$N = $ARGV[0]; # how many parts to make
if (not $N) { $N = 10 }

$Nlines = `wc -l $fastafile`;
$Nlines = $Nlines * 1; # remove file name from wc output
if ($Nlines % 2) {
  die "$fastafile has an uneven number of lines !\n";
}
$Nseq = $Nlines / 2;
if ($N > $Nseq) {
  die "there are less than $N variants, do not run in parallel\n";
}

$chunksize = int($Nlines / $N);
if ($chunksize % 2) { # need even number of lines per file (def + seq pairs)
  $chunksize--;
}
$Nline = 1;
open IN, $fastafile;
for ($i = 1 ; $i < $N ; $i++) {
  $outfilename = $fastafile;
  $outfilename =~ s/.fa$/.$i.fa/;
  open OUT, ">$outfilename";
  for ($j = $Nline ; $j <= $chunksize * $i ; $j++) {
    $line = <IN>;
    print OUT $line;
    $Nline++;
  }
  close OUT;
}

# write the last chunk
$outfilename = $fastafile;
$outfilename =~ s/.fa$/.$N.fa/;
open OUT, ">$outfilename";
for ($j = $Nline ; $j <= $Nlines ; $j++) {
  $line = <IN>;
  print OUT $line;
  $Nline++;
}
close OUT;
close IN;
