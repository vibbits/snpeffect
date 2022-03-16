#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses a BLAST output and prepares what must serve as input for foldX.
# We only consider 100% identity matches of at least length 10. We take the PDB
# sequence that makes the best hit with 100% identity (which is not always the
# very first hit !) and collect all the hits against the sequence(s) from that
# DB entry.
# Note that the early version of the script has been split in 2, some
# of the functionality has been transferred to preparseBLASToutput.pl
# It takes a input :
#   - BLASTfile_preprocessed
#   - the file PDBpositions.fa, derived from the PDB databank, with info
#     about the positions in the PDB files corresponding to the subject
#   - the file variants.tab made from the VSC file, with info about the
#     position of the mutation in the query
# It produces as output :
#   - a file variants4foldX_unsorted.tab with data about which variants must be
#     submitted to FoldX for further analysis, including the string
#     to be put in the individual_list.txt input file
#   - a file variants_without_structure_info_unsorted.txt with why
#     some variants cannot be investigated

#$BLASTfile = 'out.fmt7_blastp';
$BLASTfile = 'BLASTfile_preprocessed';
$seqfile = $ARGV[0];
if (not $seqfile) { $seqfile = 'PDBsequences.fa' }
$posfile = $ARGV[1];
if (not $posfile) { $posfile = 'PDBpositions.fa' }
$minhitlength = 10;

#open OUT, '>variants4foldX.tab';
#open LOG, '>variants_without_structure_info.txt';
open OUT, '>variants4foldX_unsorted.tab';
open LOG, '>variants_without_structure_info_unsorted.txt';

# make hashes with information anout the variants that must be
# investigated and were BLASTED against the SwitchLab PDB sequencess
open IN, 'variants.tab' or die "cannot open variants.tab\n";
while (<IN>) {
  if (not /^(\d+)\t(\w+:\d+-\d+)\t([A-Z]+\/[A-Z]+)\t(.+)\t(.+)\t([\w\.]+)\t([A-Z])\t(\d+)\t([A-Z])$/) {
    # Nvar chrom_pos mutation geneID genename seqID refAA seqpos altAA
    die "error in file variants.tab :\n$_\n";
  } else { # $1 is Nvar
    $chrom_loc[$1] = "$2\t$3\t$4\t$5";
    $seqID[$1] = $6;
    $refAA[$1] = $7;
    $seqpos[$1] = $8;
    $altAA[$1] = $9;
  }
}
close IN;

# make a hash of arrays with the numbering of the amino acids in the 
# SwitchLab PDB files, for the purpose of making a correct mutation string
# for FoldX
open IN, "$posfile" or die "cannot open $posfile\n";
while (<IN>) {
  if (/^>(.+_\w)/) {
    $NPDB_chain = $1;
  } elsif (/^\d/) {
    chomp; # to avoid having the end-of-line incorparated into the last number
    @positions = split / /;
    unshift @positions, 0; # so that $positions[1] (not 0) is the first
    for ($i = 1 ; $i <= $#positions ; $i++) {
      $NPDB_chain_pos{$NPDB_chain}[$i] = $positions[$i];
    }
  }
}
close IN;

# make a hash with the SwitchLab PDB sequences, for the purpose of testing
# whether the sequence does not differ from the SnpEff sequence in the
# position of the mutation
open IN, "$seqfile" or die "cannot open $seqfile\n";
while (<IN>) {
  if (/^>(.+_\w)/) {
    $NPDB_chain = $1;
  } elsif (/^[A-Z]/) {
    chomp;
    $NPDB_chain_seq{$NPDB_chain} = $_;
  }
}
close IN;

open IN, "$BLASTfile" or die "cannot open $BLASTfile\n";
$firstentry = 1;
$line = 0;
while (<IN>) {
  $line += 1;
  if (/^# BLASTP/) {
    if ($firstentry) { $firstentry = 0 } else { 
        print "treating entry ${line} - ${NPDB}:${mutationstring}\n"; 
       #&treatentry; 
       }
    $ismatch = 0;
  } elsif (/^# Query: (\d+)/) {
    $Nvar = $1;
  } elsif (/^(\d+)\t(.+)_(\w)\t([0-9\.]+\t(.+)\n$)/) {
      # variant number ; PDB number ; chain ; %ident ; rest of line
    if ($1 != $Nvar) { die "var number mismatch at $Nvar $1\n"; }
     
    if ( ($4 == 100.0) && ($2 eq $NPDB) ) {
      $ismatch = 1;
      $NPDB = $2 ; $chain = $3 ; $restofline = $5;
      print "match ${line} -  ${NPDB} - $mutationstring:  ${restofline}\n"; 
      $newmutationstring = &treathit;
      if ($newmutationstring) { print "new mutationstring: ${newmutationstring}\n"; &treatentry; } 
    }
    else {
      if ($4 == 100.0) {
        $ismatch = 1;
        $NPDB = $2 ; $chain = $3 ; $restofline = $5;
        print "else ${line} - ${NPDB}: ${restofline}\n"; 
        $mutationstring = &treathit;
        print "mutation string after hit: ${mutationstring}\n"; 
        if ($mutationstring) {
            &treatentry;
           }
      }
    }
  }
}

###########################################################################

sub treatentry {
  if ($ismatch) {
    if ($mutationstring) {
      chop $mutationstring; # remove trailing ,
      print OUT "$Nvar\t$chrom_loc[$Nvar]\t$seqID[$Nvar]\t$seqpos[$Nvar]\t$NPDB.pdb\t$mutationstring;\n";
    } elsif (not $aa_mismatch) {
      print LOG "var_$Nvar ($seqID[$Nvar] $seqpos[$Nvar]) has SNP outside PDB entry $NPDB.\n";
    }
  } else {
    print LOG "var_$Nvar ($seqID[$Nvar]) gave no BLAST hit with 100% identity.\n";
  }
  $mutationstring = ''; # reset mutationstring of file for next entry
  $aa_mismatch = 0;
}

sub treathit {
  $restofline =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([^\t]+)\t/;
  # length ; Nmismatches ; Ngaps ; s.start ; s.end ; PDB.start ; PDB.end ; E()
  if ($seqpos[$Nvar] >= $4 and $seqpos[$Nvar] <= $5) {
    if ($1 < $minhitlength) {
      print LOG "var_$Nvar againts PDB entry $NPDB chain $chain has only a $1 length hit\n";
    } else {
      $NPDB_chain = "${NPDB}_${chain}";
      $seq_pos = $seqpos[$Nvar] + $6 - $4;
      $NPDB_chain_pos = $NPDB_chain_pos{$NPDB_chain}[$seq_pos];
      $seq_aa = substr($NPDB_chain_seq{$NPDB_chain},$seq_pos - 1,1);
      if ($seq_aa eq $refAA[$Nvar]) {
        print "mutatationstring: ${NPDB_chain} ${seq_pos}\n";
        $mutationstring = "$refAA[$Nvar]$chain$NPDB_chain_pos$altAA[$Nvar],";
      } else {
        print LOG "var_$Nvar : $refAA[$Nvar] on position $seqpos[$Nvar] of $seqID[$Nvar] does not match $seq_aa on position $NPDB_chain_pos of PDB entry $NPDB chain $chain\n";
        $aa_mismatch = 1;
      }
    }
  }
  return $mutationstring;
}
