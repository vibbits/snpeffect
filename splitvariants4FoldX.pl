#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It splits the file protvariants4foldX.tab in a number of equal parts
#   (except for the last one) and puts them in different folders
# Its purpose is to be able to run foldX in parallel on the cluster.

$variantfile = 'protvariants4foldX.tab';
$rotabase = $ARGV[0];
if (not $rotabase or $rotabase !~ /rotabase.txt$/) {
  die "script splitvariants.pl needs rotabase.txt as first argument\n"
}
#$rotabase = '/home/guybot/rotabase.txt';
$N = $ARGV[1]; # how many parts to make
if (not $N) { $N = 20 }

$Nlines = `wc -l $variantfile`;
$Nlines = $Nlines * 1; # remove file name from wc output
if ($Nlines < $N) {
  die "there are less than $N variants, do not run in parallel\n";
}

for ($i = 1 ; $i <= $N ; $i++) {
  mkdir "FOLDXrunfolder_${i}";
  symlink $rotabase, "FOLDXrunfolder_${i}/rotabase.txt";
}

$chunksize = int($Nlines / $N);
$Nline = 1;
open IN, $variantfile;
for ($i = 1 ; $i < $N ; $i++) {
  $outfilename = "FOLDXrunfolder_${i}/$variantfile";
  open OUT, ">$outfilename";
  for ($j = $Nline ; $j <= $chunksize * $i ; $j++) {
    $line = <IN>;
    print OUT $line;
    $Nline++;
  }
  close OUT;
}

# write the last chunk
$outfilename = "FOLDXrunfolder_${N}/$variantfile";
open OUT, ">$outfilename";
for ($j = $Nline ; $j <= $Nlines ; $j++) {
  $line = <IN>;
  print OUT $line;
  $Nline++;
}
close OUT;
close IN;
