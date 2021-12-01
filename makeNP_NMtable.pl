#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It performs the last step for making the NP_* NM_* correspondence table
# to be used for picking UniProt standard sequences as well as to integrate
# PolyPhen output.
#
# procedure (last done on 3O December 2020 with data from 29 December 2020)
# 1. download from ftp.ncbi.nlm.nih.gov
#    in refseq/H_sapiens/mRNA_Prot
#    human.?.protein.faa.gz and human.?.rna.fna.gz
#    (currently 10 files each, sequences in same order, although there
#     are also ncRNA's mixed with the mRNA's)
#     and gunzip
# 2. grep '^>NP' human.*.protein.faa > NP
#    grep '^>NM' human.*.rna.fna > NM
# 3. run this script

open NP, 'NP';
open NM, 'NM';
open OUT, '>NP_NM.tab';
while (<NP>) {
  $np_line = $_;
  $np_line =~ /^human\.\d+\.protein\.faa:>(NP_\d+)\.\d+ /;
  $NP_ID = $1;
  $nm_line = <NM>;
  $nm_line =~ /^human\.\d+\.rna\.fna:>(NM_\d+)\.\d+ /;
  $NM_ID = $1;
  print OUT "$NP_ID\t$NM_ID\n";
}
