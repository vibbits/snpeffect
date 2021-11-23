#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It splits the files variants.tab, reference_sequences.fa and
# mutated_sequences.fa in a number of equal parts (except for the last one)
# and puts them in different folders.  Its purpose is to be able to run the
# AGADIR pipeline in parallel on the cluster, while obtaining protein
# aggregation data, which will have to be assembled in the final report.

$variantsfile = 'variants.tab';
$referencefile = 'reference_sequences.fa';
$mutationfile = 'mutated_sequences.fa';
$optionsfile = $ARGV[0];
$N = $ARGV[1]; # how many parts to make
if (not $N) { $N = 10 }

# check that the files have the right size
$Nlines_variantsfile = `wc -l $variantsfile`;
$Nlines_variantsfile = $Nlines_variantsfile * 1;
$Nlines_referencefile = `wc -l $referencefile`;
$Nlines_referencefile = $Nlines_referencefile * 1;
$Nlines_mutationfile = `wc -l $mutationfile`;
$Nlines_mutationfile = $Nlines_mutationfile * 1;
if ($Nlines_referencefile != $Nlines_mutationfile) {
  die "$referencefile and $mutationfile do not have the same lenght !\n";
}
if ($Nlines_referencefile % 2) {
  die "$referencefile has an uneven number of lines !\n";
}
if ($Nlines_referencefile != 2 * $Nlines_variantsfile) {
  die "$referencefile does not seem to have 2 lines (def + seq) per variant !\n";
}
if ($N > $Nlines_variantsfile) {
  die "there are less than $N variants, do not run in parallel\n";
}

for ($i = 1 ; $i <= $N ; $i++) {
  mkdir "AGADIRrunfolder_${i}";
  symlink $optionsfile, "AGADIRrunfolder_${i}/Options.txt";
}

$chunksize = int($Nlines_variantsfile / $N);
$Nline = 1;
open VAR, $variantsfile;
open REF, $referencefile;
open MUT, $mutationfile;
for ($i = 1 ; $i < $N ; $i++) {
  open OUTVAR, ">AGADIRrunfolder_${i}/$variantsfile";
  open OUTREF, ">AGADIRrunfolder_${i}/$referencefile";
  open OUTMUT, ">AGADIRrunfolder_${i}/$mutationfile";
  for ($j = $Nline ; $j <= $chunksize * $i ; $j++) {
    $line = <VAR>;
    print OUTVAR $line;
    $line = <REF>;
    $line =~ s/^>(\d+)\t.*$/>\1/;
    print OUTREF $line;
    $line = <REF>;
    print OUTREF $line;
    $line = <MUT>;
    $line =~ s/^>(\d+)\t.*$/>\1/;
    print OUTMUT $line;
    $line = <MUT>;
    print OUTMUT $line;
    $Nline++;
  }
  close OUTVAR ; close OUTREF ; close OUTMUT;
}

# write the last chunk
open OUTVAR, ">AGADIRrunfolder_${N}/$variantsfile";
open OUTREF, ">AGADIRrunfolder_${N}/$referencefile";
open OUTMUT, ">AGADIRrunfolder_${N}/$mutationfile";
for ($j = $Nline ; $j <= $Nlines_variantsfile ; $j++) {
  $line = <VAR>;
  print OUTVAR $line;
  $line = <REF>;
  $line =~ s/^>(\d+)\t.*$/>\1/;
  print OUTREF $line;
  $line = <REF>;
  print OUTREF $line;
  $line = <MUT>;
  $line =~ s/^>(\d+)\t.*$/>\1/;
  print OUTMUT $line;
  $line = <MUT>;
  print OUTMUT $line;
  $Nline++;
}
close VAR ; close REF ; close MUT;
close OUTVAR ; close OUTREF ; close OUTMUT;
