#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It splits the file reference_sequences_nonredundant.fa in a number of
# equal parts (except for the last one) and puts them in different folders
# Its purpose is to be able to run the Gene3D pipeline in parallel
# on the cluster, while obtaining domain data for the different
# reference sequences, which will have to be assembled in the final report

$fastafile = 'reference_sequences_nonredundant.fa';
$gene3d = $ARGV[0];
if (not $gene3d or $gene3d !~ /gene3d_hmmsearch$/) {
  die "script splitreferences4Gene3D.pl
needs folder gene3d_hmmsearch as first argument\n"
}
#$gene3d = '/home/guybot/gene3d_hmmsearch';
$N = $ARGV[1]; # how many parts to make
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
  mkdir "Gene3Drunfolder_${i}";
  symlink "$gene3d/discontinuous", "Gene3Drunfolder_${i}/discontinuous";
  symlink "$gene3d/cath_release", "Gene3Drunfolder_${i}/cath_release";
}

$chunksize = int($Nlines / ($N * 2) );
  # chunksize is number of sequences (defline+seqline) per chunk
$Nseq = 1;
open IN, $fastafile;
for ($i = 1 ; $i < $N ; $i++) {
  $outfilename = "Gene3Drunfolder_${i}/$fastafile";
  open OUT, ">$outfilename";
  for ($j = $Nseq ; $j <= $chunksize * $i ; $j++) {
    $line = <IN>;
    $line =~ s/^>\d+\t/>/;
    print OUT $line;
    $line = <IN>;
    print OUT $line;
    $Nseq++;
  }
  close OUT;
}

# write the last chunk
$outfilename = "Gene3Drunfolder_${N}/$fastafile";
open OUT, ">$outfilename";
for ($j = $Nseq ; $j <= $Nseqs ; $j++) {
  $line = <IN>;
  $line =~ s/^>\d+\t/>/;
  print OUT $line;
  $line = <IN>;
  print OUT $line;
  $Nseq++;
}
close OUT;
close IN;
