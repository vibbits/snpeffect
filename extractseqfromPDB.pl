#!/usr/bin/perl
# This script has been written by BioinformaticsCore for switchLab
# It parses the PDB files in the specified folder
# and makes a fastA file with the sequences
# It leaves out the following :
#   - non-aa molecules/subunits (curently these are : ADP GDP GTP RET )
#   - chains that are entirely composed of alanine
#   - "chains" that have a redundant identifier, usually containing only
#      non-aa subunits or a short peptide
#   - chains with redundant residue numbering
# It also makes a file with the position numbers of the amino acids
# in the PDB file which will later be needed to find the position number
# of the mutation

%aa1 = ('ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N', 'ASP' => 'D', 'CYS' => 'C',
  'PHE' => 'F', 'GLN' => 'Q', 'GLU' => 'E', 'GLY' => 'G', 'HIS' => 'H',
  'ILE' => 'I', 'LEU' => 'L', 'LYS' => 'K', 'MET' => 'M', 'PRO' => 'P',
  'SER' => 'S', 'THR' => 'T', 'TRP' => 'W', 'TYR' => 'Y', 'VAL' => 'V',
   # unusal amino acids supported by foldX
   'H1S' => 'H', 'H2S' => 'H', 'H3S' => 'H', 'HYP' => 'P', 'M3L' => 'K',
   'MLY' => 'K', 'MLZ' => 'K', 'PTR' => 'T', 'SEP' => 'S', 'TPO' => 'Y',
   'TYS' => 'Y');
   
$path = $ARGV[0];
if ($path eq '') {
  warn "You need to specify the directory containing the PDBs\n";
}
   
opendir my $dir, $ARGV[0] or die "Cannot open directory";
my @files = readdir $dir;
closedir $dir;

open OUTFA, '>PDBsequences.fa';
open OUTNR, '>PDBpositions.fa';
foreach (@files) {
	if ($_ =~ /.*\.pdb/) {
		$_ =~ /(.*)\.pdb/;
		$filename = $1;
	    if (open IN, "$path/$_") {
	    $PDBnumerror = 0;
	    $treatedchains = '';
	    $line = '';
	    # search for first ATOM line and initialize
	    until ($line =~ /^ATOM +\d+ +.* +([A-Z][A-Z][A-Z]) +(\w) *(\d+)/) {
	      # $1 is aa 3 letter code, $2 is chain identifier, $3 is position
	      # note : some chain identifiers are not A-Z and sometimes there
	      #   is no space between the identifier and the position number
	      $line = <IN>;
	    }
	    $line =~ /^ATOM +\d+ +.* +([A-Z][A-Z][A-Z]) +(\w) *(\d+)/;
	    $AA3 = $1 ; $chain = $2 ; $pos = $3;
	    $currentchain = $chain;
	    $currentpos = $pos;
	    $header = ">${filename}_${currentchain}";
	    if (exists $aa1{$AA3}) {
	      $sequence = $aa1{$AA3};
	      $POSITIONS = $currentpos;
	    } else {
	      warn "unknown amino acid $AA3\n";
	      $sequence = '';
	      $POSITIONS = '';
	    }
	    while (<IN>) {
	      if (/^ATOM +\d+ +.* +([A-Z][A-Z][A-Z]) +(\w) *(\d+)/) {
	        $AA3 = $1 ; $chain = $2 ; $pos = $3;
	        if ($chain ne $currentchain) { # begin new chain
	          if ($treatedchains =~ /$currentchain/) {
	            if (length $sequence > 0 and $sequence =~ /[^A]/) {
	              warn "chain $currentchain of $_ has redundant label\n";
	            }
	          } else {
	            if (length $sequence > 0 and $sequence =~ /[^A]/
	                and not $PDBnumerror) {
	              print OUTFA "$header\n$sequence\n";
	              $POSITIONS =~ s/^ ?//;
	              print OUTNR "$header\n$POSITIONS\n";
	            }
	            $treatedchains .= $currentchain;
	          }
	          $PDBnumerror = 0;
	          $currentchain = $chain;
	          $currentpos = $pos;
	          $header = ">${filename}_${currentchain}";
	          if (exists $aa1{$AA3}) {
	            $sequence = $aa1{$AA3};
	            $POSITIONS = $currentpos;
	          } else {
	            warn "unknown amino acid $AA3\n";
	            $sequence = '';
	            $POSITIONS = '';
	          }
	        } else {
	          if ($pos > $currentpos) { # next amino acid
	            $currentpos = $pos;
	            if (exists $aa1{$AA3}) {
	              $sequence .= $aa1{$AA3};
	              $POSITIONS .= " $currentpos";
	            } else {
	              warn "unknown amino acid $AA3\n";
	            }
	          } elsif ($pos < $currentpos) { # erroneous PDB file
	            $PDBnumerror = 1;
	            warn "no consistently increasing aa numbers in $_ chain $currentchain\n";
	          }
	        }
	      }
	    }
	    close IN;
	    # handle last sequence in file
	    if ($treatedchains =~ /$currentchain/) {
	      if (length $sequence > 0 and $sequence =~ /[^A]/) {
	        warn "chain $currentchain of $_ has redundant label\n";
	      }
	    } else {
	      if (length $sequence > 0 and $sequence =~ /[^A]/
	          and not $PDBnumerror) {
	        print OUTFA "$header\n$sequence\n";
	        $POSITIONS =~ s/^ ?//;
	        print OUTNR "$header\n$POSITIONS\n";
	      }
	    }
	  } else {
	    warn "cannot open $_\n";
	  }
   }
   else {
   	    warn "$_ is not a PDB file\n";
   }
}
