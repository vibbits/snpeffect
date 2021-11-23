#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It searches for (predicted) structural domains in protein sequences,
# using a HMMER version 3 search against CATH/Gene3D domain models.
# It takes as input reference_sequences_nonredundant.fa and
# produces as output seqs.domains with the domains found
# and a file no_domain_info.txt with the sequences
# for which no domains could be found.

$run_under_master_script = 1;

$python = $ARGV[0];
#$python = '/usr/bin/python';
if (not $python or $python !~ /[^\/].*python2?$/) {
  die "script finddomains.pl needs python as first argument\n"
}
$hmmsearch = $ARGV[1];
#$hmmsearch = '/home/guybot/HMMER3/bin/hmmsearch';
if (not $hmmsearch or $hmmsearch !~ /[^\/].*hmmsearch$/) {
  die "script finddomains.pl needs hmmsearch as second argument\n"
}
$gene3d = $ARGV[2];
#$gene3d = '/home/guybot/gene3d_hmmsearch';
if (not $gene3d or $gene3d !~ /[^\/].*gene3d_hmmsearch$/ or not -d $gene3d ) {
  die "script finddomains.pl needs folder gene3d_hmmsearch as third argument\n"
}

# test if we have everything needed
if (not -e 'reference_sequences_nonredundant.fa') { &setflag ;
  die "reference_sequences_nonredundant.fa not found\n" }
if (not -e 'discontinuous') { &setflag ; die "discontinuous not found\n" }
if (not -e 'cath_release') { &setflag ; die "cath_release not found\n" }

open REF, 'reference_sequences_nonredundant.fa';
open WARN, '>no_domain_info.txt';
system 'touch seqs.domains';

while(<REF>) {
  # put sequence in temporary file
  if (not /^>([\w\.]+)/) {
    &setflag;
    die "$_ does not look like fastA description line for transcript ID\n";
  }
  $transcriptID = $1;
  open SEQ, '>tempseq.fa';
  print SEQ $_;
  $sequence = <REF>;
  print SEQ $sequence;
  close SEQ;

  # search the sequence against the Markov models in Gene3D
  system "$hmmsearch -Z 10000000 --domE 0.001 --incdomE 0.001 -o tempseq.hmmsearch $gene3d/hmms/main.hmm tempseq.fa";
  unlink 'tempseq.fa';

  # make a domain description of the sequence by making a set of
  # non-overlapping best hits
  system "$gene3d/cath-resolve-hits --min-dc-hmm-coverage=80 --worst-permissible-bitscore 25 --output-hmmer-aln --input-format hmmsearch_out tempseq.hmmsearch > tempseq.domains";
  unlink 'tempseq.hmmsearch';

  # eventually add CATH family information to domain description
  if (-z 'tempseq.domains') {
    print WARN "$transcriptID : no Gene3D domains found\n";
  } else {
    system "$python $gene3d/assign_cath_superfamilies.py tempseq.domains";
    $Nlines_tempseqdomainscsv = `wc -l tempseq.domains.csv`;
    $Nlines_tempseqdomainscsv = $Nlines_tempseqdomainscsv * 1;
    if ($Nlines_tempseqdomainscsv == 1) {
      print WARN "$transcriptID : no Gene3D domains reported because assign_cath_superfamilies.py BUG\n";
    } else {
      system 'cat tempseq.domains.csv >> seqs.domains';
    }
    unlink 'tempseq.domains.csv';
  }
  unlink 'tempseq.domains';
}

# terminate
&setflag;

sub setflag {
  # to indicate that script has terminated, notwithstanding it crashed or
  # ended normally ; useful if script is started under a masterscript
  open FLAG, '>flag.Gene3D';
  close FLAG;
}
