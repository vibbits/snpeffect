#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It takes as input :
#   - the fastA file out.fa written by SnpEff
#   - the file snpEff_genes.txt made by SnpEff (needed to find the gene name)
#  A number of sequences are excluded :
#   - entries without information about the substitution or (bizarrily)
#     where the original and new aa are the same (noncoding RNA, mutation
#     outside coding sequence, silent mutation)
#   - mutations that involve creating or removing a stop codon
#   - mutations of initiator methionine (protein expression uncertain)
#   - entries where the reference and/or mutated sequence contain internal
#     stop codons (sequence maybe not expressed or not well known, anyway
#     cannot be handled by FoldX or TANGO)
#   - protein sequences longer than 10.000 aa (cannot be handled by
#     AGADIR)
# It writes as output :
#   - reference_sequences.faa : fastA file with the sequences as they are
#     in the reference genome
#   - reference_sequences_nonredundant.faa : a nonredundant version of this file
#   - sequence_identities.tab : info about the redundant sequences
#   - mutated_sequences.fa : a fastA file with the mutated sequences
#   - variants.tab : a file with info about the variants that need
#     further analysis
#   - variants_not_investigated.txt : info why some variants were excluded
#     (by default only those that are too long)
#  The selected variants get a number that follows their order in the
#  SnpEff output and roughly follows their position on the genome.
#
#  There is an option to reduce the output by investigating only "standard"
#  sequences. The files with info about the standard can be made with
#  makeNM_ACtable.pl. Note that the SnpEff database contains "intermediate"
#  entries like NM_001101663.4.2. These will necessarily be ignored.

%aa1 = ('Ala' => 'A', 'Arg' => 'R', 'Asn' => 'N', 'Asp' => 'D', 'Cys' => 'C',
  'Phe' => 'F', 'Gln' => 'Q', 'Glu' => 'E', 'Gly' => 'G', 'His' => 'H',
  'Ile' => 'I', 'Leu' => 'L', 'Lys' => 'K', 'Met' => 'M', 'Pro' => 'P',
  'Ser' => 'S', 'Thr' => 'T', 'Trp' => 'W', 'Tyr' => 'Y', 'Val' => 'V');

$MAXLEN = $ARGV[0];
if (not $MAXLEN) { $MAXLEN = 10000 }
$SwissProtstandard = $ARGV[1];
#$SwissProtstandard = '../scripts/NM_AC_ID.tab';
$translationequivalents = $ARGV[2];
#$translationequivalents = '../scripts/NM_alternatives.tab';
if (@ARGV == 2) {
  die "If you provide a table of standard protein sequences,\nyou must also provide a table of mRNA's that produce the same protein.\n";
}

# make a hash with the gene ID and name
open IN, 'snpEff_genes.txt' or die "cannot open file snpEff_genes.txt\n";
while (<IN>) {
  if (not /^#/) {
    if (not /^([^\t]+)\t([^\t]+)\t([^\t]+)\t/) {
      # GeneName  GeneId  transcriptId
      die "error in file snpEff_genes.txt :\n$_\n";
    }
    $genename = $1 ; $geneID = $2 ; $transcriptID = $3;
    $genename{$transcriptID} = $genename;
    $geneID{$transcriptID} = $geneID;
  }
}
close IN;

# make, if available, a hash with standard SwissProt entries
# and a hash with mRNA equivalences
if ($SwissProtstandard) {
  open IN, $SwissProtstandard or die "cannot open file $SwissProtstandard\n";
  while (<IN>) {
    chomp;
    $_ =~ /^(.+)\t(.+)\t(.+)$/;
    $standard4refseqID{$1} = $2;
    $genefield = $3;
    @items = split '/', $genefield;
    for $item (@items) {
      $standard4geneexists{$item} = 1;
    }
  }
  close IN;

  open IN, $translationequivalents
    or die "cannot open file $translationequivalents\n";
  while (<IN>) {
    chomp;
    $_ =~ /^(.+)\t(.+)$/;
    $altNM{$2} = $1;
  }
  close IN;
}


# parse the SnpEff output
$N = 0;
open IN, 'out.fa' or die "cannot open out.fa\n";
open OUTFA, '>reference_sequences.fa';
open OUTFA2, '>mutated_sequences.fa';
open OUTPOS, '>variants.tab';
open LOG, '>variants_not_investigated.txt';
$firstline = <IN>;
$firstline =~ /^>([^ ]+) Ref/;
$ID = $1;
$reading_refseq = 1;
while (<IN>) {
  chomp;
  #s/^(>NM_\d+\.\d+)\.\d+ /\1 /;
  if (/^>([^ ]+) Ref/) {
    &treat_sequence_pair;
    $ID = $1;
    $reading_refseq = 1;
    $refseq = '';
  } elsif (/>([^ ]+) Variant/) {
    if ($1 ne $ID) {
      die "ID mismatch for $_\n";
    }
    $vardef = $_ ; chomp $vardef ; $vardef =~ s/^>//;
    $reading_refseq = 0;
    $varseq = '';
  } elsif ($reading_refseq) {
    $refseq .= $_;
  } else { # reading variant seq
    $varseq .= $_;
  }
}
&treat_sequence_pair;
close IN ; close OUTFA ; close OUTFA2 ; close OUTPOS ; close LOG;

# make a nonredundant reference sequence file and a hash with correspondence
# between first sequence found and its copies
open IN, 'reference_sequences.fa';
open OUT, '>reference_sequences_nonredundant.fa';
while (<IN>) {
  $_ =~ /^>(\d+)\t(.+)$/;
  $defline = $_ ; $Nvar = $1 ; $ID = $2;
  $seq = <IN>;
  if (exists $Nprimary{$ID}) {
    $Nsecondaries{$ID} .= "\t$Nvar";
  } else {
    $Nprimary{$ID} = $Nvar;
    print OUT "$defline$seq";
  }
}
close IN, close OUT;

# make the table with correspondencies
open TABLE, '>sequence_identities_unsorted.tab';
for $k (keys %Nprimary) {
  print TABLE $Nprimary{$k}, $Nsecondaries{$k}, "\n";
}
close TABLE;

########################################################################

sub treat_sequence_pair {
  # reject entries with no info about amino acid change
  # or if ref and alt aa are the same
  if ($vardef =~ /^NM_[^ ]+ Variant ([^ ]+) Ref:([^ ]+) Alt:([^ ]+) HGVS\.p:p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])$/) {
    $chrom_pos = $1;
    $mutation = "$2/$3";
    if ($4 ne $6) {
      $refA = $aa1{$4} ; $pos = $5 ; $altA = $aa1{$6};
      $refseq =~ s/\*$//;
      $varseq =~ s/\*$//;
        # remove terminal * (= stop codon)
      if ($refseq !~ /[*?]/ and $varseq !~ /[*?]/) {
        # reject sequences with internal * or terminal ?
        $length = length $refseq;
        if ($length <= $MAXLEN) {
          # reject sequences longer than the limit of AGADIR
          if ($SwissProtstandard) {
            #$ID_noextension = $ID ; $ID_noextension =~ s/\.\d+$//; # Old syntax: this regex gives the string up until the last period.
            $ID_noextension = $ID ; $ID_noextension = substr($ID_noextension, 0, index($ID_noextension, ".")); # This gives the string up until the first period.
            $varkey = $vardef ; $varkey =~ s/^.+Variant /Variant/;
            #if ($ID eq "NM_001367552.1") {print "$ID\n$ID_noextension\n$standard4geneexists{$genename{$ID_noextension}}\n$standard4refseqID{$ID_noextension}\n"; die ;}
            if (exists $standard4refseqID{$ID_noextension}) {
              if ($already_handled{$varkey}) {
                print LOG "$chrom_pos : $ID is handled with an alternative transcript\n";
              } else {
                $N++;
                print OUTFA ">$N\t$ID\n$refseq\n";
                print OUTFA2 ">$N\t$vardef\n$varseq\n";
                print OUTPOS "$N\t$chrom_pos\t$mutation\t$geneID{$ID}\t$genename{$ID}\t$ID\t$refA\t$pos\t$altA\n";
                $already_handled{$varkey} = 1;
              }
            } elsif (exists $altNM{$ID_noextension}) {
              if ($already_handled{$varkey}) {
                print LOG "$chrom_pos : $ID produces the same protein as $altNM{$ID_noextension}\n";
              } else {
                $N++;
                print OUTFA ">$N\t$ID\n$refseq\n";
                print OUTFA2 ">$N\t$vardef\n$varseq\n";
                print OUTPOS "$N\t$chrom_pos\t$mutation\t$geneID{$ID}\t$genename{$ID}\t$ID\t$refA\t$pos\t$altA\n";
                $already_handled{$varkey} = 1;
              }
            } elsif ($standard4geneexists{$genename{$ID}}) {
              print LOG "$chrom_pos : $ID does not produce a UniProt standard protein\n";
            } else {
              print LOG "$chrom_pos : $ID : $genename{$ID} does not have a UniProt standard protein (or SnpEff and UniProt use a different name for the same gene)\n";
            }

          } else {
            $N++;
            print OUTFA ">$N\t$ID\n$refseq\n";
            print OUTFA2 ">$N\t$vardef\n$varseq\n";
            print OUTPOS "$N\t$chrom_pos\t$mutation\t$geneID{$ID}\t$genename{$ID}\t$ID\t$refA\t$pos\t$altA\n";
          }

        } else {
          print LOG "$chrom_pos : $ID has length $length > $MAXLEN\n";
        }
      }
    }
  }
}
