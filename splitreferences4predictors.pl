#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It splits the file reference_sequences_nonredundant.fa in a number of
# equal parts (except for the last one) and puts them in different folders.
# It also makes from variants.tab files with protein name - mutation string
# pairs and puts them in the folders next to the corresponding sequences.
# Its purpose is to be able to run SIFT and PROVEAN in parallel on the
# cluster.

$fastafile = 'reference_sequences_nonredundant.fa';
$N = $ARGV[0]; # how many parts to make
if (not $N) { $N = 20 }

$Nlines = `wc -l $fastafile`;
$Nlines = $Nlines * 1; # remove file name from wc output
if ($Nlines % 2) {
  die "$fastafile has an uneven number of lines !\n";
}
$Nseqs = $Nlines / 2;
if ($N > $Nseqs) {
  die "there are less than $N variants, do not run in parallel\n";
}

for ($i = 1 ; $i <= $N ; $i++) {
  mkdir "runfolder4predictors_${i}";
  mkdir "runfolder4predictors_${i}/PROVEAN_TMP";
  mkdir "runfolder4predictors_${i}/SIFT_TMP";
}

# make hash of hashes with sequences IDs and corresponding mutations
open IN, 'variants.tab' or die "cannot open variants.tab\n";
while (<IN>) {
  if (not /^(\d+)\t(\w+:\d+-\d+)\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
    die "error in file variants.tab :\n$_\n";
  }
  $transcriptID = $5;
  $refAA = $6;
  $POS = $7;
  $mutAA = $8;
  $mutations{$transcriptID}{"$refAA$POS$mutAA"} = 1;
}
close IN;

$chunksize = int($Nlines / ($N * 2) );
  # chunksize is number of sequences (defline+seqline) per chunk
$Nseq = 1;
open IN, $fastafile;
for ($i = 1 ; $i < $N ; $i++) {
  $sequences_out = "runfolder4predictors_${i}/$fastafile";
  $variants_out = "runfolder4predictors_${i}/variants";
  open SEQOUT, ">$sequences_out";
  open OUT, ">$variants_out";
  for ($j = $Nseq ; $j <= $chunksize * $i ; $j++) {
    $line = <IN>;
    $line =~ /^>\d+\t(.+)\n$/;
    $transcriptID = $1;
    $line =~ s/^>\d+\t/>/;
    print SEQOUT $line;
    $line = <IN>;
    print SEQOUT $line;
    foreach $var (keys  $mutations{$transcriptID} ) {
      print OUT "$transcriptID\t$var\n";
    }
    $Nseq++;
  }
  close SEQOUT;
  close OUT;
}

# write the last chunk
$sequences_out = "runfolder4predictors_${N}/$fastafile";
$variants_out = "runfolder4predictors_${N}/variants";
open SEQOUT, ">$sequences_out";
open OUT, ">$variants_out";
for ($j = $Nseq ; $j <= $Nseqs ; $j++) {
  $line = <IN>;
  $line =~ /^>\d+\t(.+)\n$/;
  $transcriptID = $1;
  $line =~ s/^>\d+\t/>/;
  print SEQOUT $line;
  $line = <IN>;
  print SEQOUT $line;
  foreach $var (keys  $mutations{$transcriptID} ) {
    print OUT "$transcriptID\t$var\n";
  }
  $Nseq++;
}
close SEQOUT;
close OUT;
close IN;
