#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It makes the NM_* "standard SwissProt AC" table
# to be used for selecting just one transcript per gene
#
# procedure (last done on 30 December with UniProt from 2 December 2020) :
# 1. download from ftp.uniprot.org
#    in pub/databases/uniprot/current_release/knowledgebase/
#          reference_proteomes/Eukaryota
#    UP000005640_9606.fasta.gz   (is for Homo sapiens)  and gunzip
#    -> 20600 entries
# 2. download from ftp.ncbi.nlm.nih.gov
#    in refseq/H_sapiens/mRNA_Prot
#    human.?.protein.faa.gz (currently 9 files)      and gunzip
#    -> 46121800 entries
# 3. run this script to get only entries like
#      NM_001008219       P04745/P0DTE7/P0DTE8    AMY1A/AMY1B/AMY1C
#    as well as a table with alternative transcripts
#      NM_001008219 NM_001008218
#      NM_001008219 NM_001008221
#      ...


open IN, 'NP_NM.tab' or die "cannot open NP_NM.tab\n";
while (<IN>) {
  $_ =~ /(.+)\t(.+)$/;
  $NM{$1} = $2;
}
close IN;

open OUT, '>temp_NM_alternatives.tab';
for $n (1..9) {
  open IN, "human.$n.protein.faa" or die "cannot open human.$n.protein.faa\n";
  $line = <IN>;
  $line =~ /^>([A-Z]P_\d+)/;
  $ID = $1;
  while (<IN>) {
    if (/^>([A-Z]P_\d+)/) {
      if (substr($ID,0,3) eq 'NP_') {
        if (exists ($NM{$ID})) {
          if (exists $RefSeq{$seq}) {
            print OUT "$RefSeq{$seq}\t$NM{$ID}\n";
          } else {
            $RefSeq{$seq} = $NM{$ID};
          }
        }
      }
      $seq = '';
      $ID = $1;
    } else {
      chomp;
      $seq .= $_;
    }
  }
  if ($ID =~ /^NP_/) {
    if (substr($ID,0,3) eq 'NP_') {
      if (exists ($NM{$ID})) {
        if (exists $RefSeq{$seq}) {
          print OUT "$RefSeq{$seq}\t$NM{$ID}\n";
        } else {
          $RefSeq{$seq} = $NM{$ID};
        }
      }
    }
  }
  close IN;
}
close OUT;

open IN, 'UP000005640_9606.fasta' or die "cannot open UP000005640_9606.fasta\n";
open LOG, '>no_entry_in_RefSeq.txt';
$line = <IN>;
$line =~ /^>(sp|tr)\|([^|]+)\|.+ GN=([^ ]+) /;
$UniProtID = $2;
$geneID = $3;
while (<IN>) {
  if (/^>(sp|tr)\|([^|]+)\|.+ GN=([^ ]+) /) {
    if (exists $RefSeq{$seq}) {
      if (exists $UniProtID{$RefSeq{$seq}}) {
        $UniProtID{$RefSeq{$seq}} .= "/$UniProtID";
        $geneID{$RefSeq{$seq}} .= "/$geneID";
      } else {
        $UniProtID{$RefSeq{$seq}} = $UniProtID;
        $geneID{$RefSeq{$seq}} = $geneID;
      }
    } else {
      print LOG "$UniProtID\t$geneID\n";
    }
    $seq = '';
    $UniProtID = $2;
    $geneID = $3;
  } else {
    chomp;
    $seq .= $_;
  }
}
if (exists $RefSeq{$seq}) {
  if (exists $UniProtID{$RefSeq{$seq}}) {
    $UniProtID{$RefSeq{$seq}} .= "/$UniProtID";
    $geneID{$RefSeq{$seq}} .= "/$geneID";
  } else {
    $UniProtID{$RefSeq{$seq}} = $UniProtID;
    $geneID{$RefSeq{$seq}} = $geneID;
  }
} else {
  print LOG "$UniProtID\t$geneID\n";
}
close IN ; close LOG;

open OUT, '>NM_AC_ID.tab';
for $key (keys %UniProtID) {
  print OUT "$key\t$UniProtID{$key}\t$geneID{$key}\n";
}
close OUT;

open IN, 'temp_NM_alternatives.tab';
open OUT, '>NM_alternatives.tab';
while(<IN>) {
  $_ =~ /^(.+)\t(.+)$/;
  if (exists $UniProtID{$1}) {
    print OUT;
  }
}

unlink 'temp_NM_alternatives.tab';
