#!/usr/bin/perl
# This script has been written by BioinformaticsCore for SwitchLab.
# The purpose is to analyze SNP's to see if they change the structure and
# tendency to aggregate of a protein.
# This is the master script that sends the complete workflow to the GRID.
# It takes as input a VCF file and produces output in tab-separated format.
# It needs, besides the scripts, the following software :
#   - SnpEff (and a SnpEff database with the same genome as that used
#       to generate the VCF file)
#   - FoldX (and the a databse containing PDB files)
#   - AGADIR
#   - Gene3D and HMMER3 (and the included database of CATH protein
#     domain Markov models)
#   - TMHMM
#   - PROVEAN, CD-HIT, BLAST version 2.4.0+ and the August 2011 release
#     of the NCBI NR protein databank (from PROVEAN site)
#   - SIFT and the 11-Jan-2011 release of UniRef90

# To use this script, you need to make a folder/directory, put in that folder
# a file called in.vcf and then with that folder as working directory
# execute the command :
#   qsub -cwd -b y ./masterscript.pl

# NOTE : in case the workflow crashes you do not need to start again
# from te beginning, after fixing the problem you can start from the first
# failed step by providing an argument, e.g. :
#     qsub -cwd -b y masterscript.pl SNPEFF
# The script performs the following steps :
#   SNPEFF PARSESNPEFF BLAST PARSEBLAST FOLDX PARSEFOLDX GENE3D
#       AGADIR PARSEAGADIR TMHMM SIFT_PROVEAN FINAL

#The following definitions MUST be provided by the user
#The user MUST specify the correct pathway for each tool

#path of the directory containing the masterscript.pl
$scriptdir = '';  # Example: '/switchlab/group/ramdur/snpeffect-master'
#path to snpEff software
$SnpEffjar = '';  # Example: '/switchlab/group/guybot/snpEff/snpEff.jar'
#path to blastp
$BLAST = '';  # Example: '/switchlab/group/guybot/ncbi-blast-2.10.0+/bin/blastp'
#path to BLAST DB of PDB structures
$BLASTDIR = '';  # Example: '/switchlab/group/ramdur/Alphafold_Human/DB'
#path to FoldX software
$FoldX = ''; # Example: '/switchlab/group/tools/FoldX_2015/FoldX'
#path to directory containing PDB structures
$PDBDIR = ''; # Example: '/switchlab/group/ramdur/Alphafold_Human/RepairPDBs'
#path to python
$python = ''; # Example: '/usr/bin/python'
#path to hmmsearch
$hmmsearch = ''; # Example: '/switchlab/group/guybot/HMMER3/bin/hmmsearch'
#path to gene3d_hmmsearch
$gene3d = ''; # Example: '/switchlab/group/guybot/gene3d_hmmsearch'
#path to tmhmm
$tmhmm = ''; # Example: '/switchlab/group/guybot/tmhmm-2.0c/bin/tmhmm'
#path to sift directory
$siftdir = ''; # Example: '/switchlab/group/guybot/sift6.2.1'
#path to blast ncbi-blast 2.4.0+
$BLAST4SIFTdir = ''; # Example: '/switchlab/group/guybot/ncbi-blast-2.4.0+'
#path to Uniref90 sequences
$UniRef90_BLASTDB = ''; # Example: '/switchlab/group/guybot/UniRef/uniref90.fa'
#path to PROVEAN software
$PROVEAN = ''; # Example: '/switchlab/group/guybot/provean-1.1.5/bin/provean.sh'

# The following definitions might need to be adapted :
### THIS IS ESPECIALLY TRUE FOR $SnpEffgenome!!
$SwissProtstandard = "$scriptdir/NM_AC_ID.tab $scriptdir/NM_alternatives.tab";
#$SwissProtstandard = ''; # put this to analyze all variants, not just
  # one transcript per gene, corresponding to the UniProt standard
$SnpEffgenome = 'hg38'; #SnpEff genome installed
  # see SnpEff software for how to install other genomes
#Different elements within the BLAST DB. The user needs to make the blast DB.
$BLASTDB = "$BLASTDIR/PDB";
$PDBseqfile = "$BLASTDIR/PDBsequences.fa";
$PDBposfile = "$BLASTDIR/PDBpositions.fa";

#Other definitions that the user can change
$MAXLEN = 10000; # maximum allowed length for proteins (mainly for AGADIR)
$Nparts4BLAST = 5; # in how many parts to split the file
  # reference_sequences_nonreduncant.fa before sending BLAST job to GRID
$Nparts4FoldX = 50; # in how many parts to split the file variants4foldX.tab
  # before sending FoldX job to GRID
$Nparts4Gene3D = 20; # in how many parts to split the file
  # reference_sequences_nonredundant.fa before sending Gene3D pipeline to GRID
$Nparts4AGADIR = 10; # in how many parts to split the file variants.tab
  # before sending AGADIR pipeline to GRID
$Nparts4predictors = 50; # in how many parts to split the file
  # reference_sequences_nonredundant.fa before sending SIFT/PROVEAN job to GRID
$MEM = '2G'; # node memory allocation for most jobs
$JAVAMEM = '20G'; # node memory allocation for Java (for SnpEff)
$BLASTMEM = '4G'; # node memory allocation for BLAST
$FOLDXMEM = '20G'; # node memory allocation for FoldX
$rotabase = "$scriptdir/rotabase.txt"; #path to rotabase.txt
$agadir = "$scriptdir/agadirwrapper"; #path to agadirwrapper

$step = $ARGV[0];
if ($step ne '') {
  goto $step;
}

SNPEFF: # add extra annotation to VCF file and generate fastA file
open SCRIPT, '>SnpEffscript.sh';
  print SCRIPT "java -jar $SnpEffjar -i vcf -o vcf -fastaProt out.fa $SnpEffgenome in.vcf > out.vcf\n";
close SCRIPT;
if (system "qsub -cwd -sync y -l mem_limit=$JAVAMEM SnpEffscript.sh") {
  die "problem with step SNPEFF\n";
}
$WARNING = `grep WARNING_REF_DOES_NOT_MATCH_GENOME out.vcf | head`;
if ($WARNING) {
  die "$WARNING...\nout.vcf contains lines with WARNING_REF_DOES_NOT_MATCH_GENOME\nAre you sure the reads were mapped to $SnpEffgenome ?\n";
}
system "mv SnpEffscript.sh.e* SnpEff_errors.txt";
unlink 'SnpEffscript.sh';
unlink 'snpEff_summary.html'; # IF NEEDED OUTCOMMENT
# workflow output : out.vcf out.fa snpEff_genes.txt SnpEff_errors.txt

PARSESNPEFF: # parse out.fa to find SNP's that cause aa substitution
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseSnpEfffasta.pl $MAXLEN $SwissProtstandard") {
  die "problem with step PARSESNPEFF parseSnpEfffasta.pl\n";
}
system 'sort -n -k 1 sequence_identities_unsorted.tab > sequence_identities.tab';
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseVCF.pl") {
  die "problem with step PARSESNPEFF parseVCF.pl\n";
}
unlink 'sequence_identities_unsorted.tab';
unlink 'out.fa'; # IF NEEDED OUTCOMMENT
unlink 'out.vcf'; # IF NEEDED OUTCOMMENT
# workflow output :
#   variants.tab reference_sequences.fa mutated_sequences.fa
#   reference_sequences_nonredundant.fa sequence_identities.tab
#   SNPAAchange.vcf variants_not_investigated.txt

BLAST: # BLAST reference_sequences_nonredundant.fa againts switchLab sequences
$Nlines = `wc -l reference_sequences_nonredundant.fa`;
$Nlines = $Nlines * 1; # remove file name from wc output
$Nseq = $Nlines / 2;
if ($Nparts4BLAST > $Nseq) { $Nparts4BLAST = 1 }
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/splitfasta.pl $Nparts4BLAST") {
  die "problem with step BLAST splitfasta.pl\n";
}
for ($i = 1 ; $i <= $Nparts4BLAST ; $i++) {
  open SCRIPT, ">BLASTscript.${i}.sh";
    print SCRIPT "$BLAST -query reference_sequences_nonredundant.${i}.fa -db $BLASTDB -seg no -ungapped -matrix BLOSUM80 -evalue 0.01 -dbsize 10000 -comp_based_stats F -outfmt 7 > out.${i}.fmt7_blastp\n";
    print SCRIPT "touch flag.BLAST.${i}\n";
  close SCRIPT;
  system "qsub -cwd -l mem_limit=$BLASTMEM BLASTscript.${i}.sh";
}
$DONE = 0;
until ($DONE) {
  $DONE = 1;
  for ($i = 1 ; $i <= $Nparts4BLAST ; $i++) {
    if (not -e "flag.BLAST.${i}") { $DONE = 0 }
  }
}
unlink 'out.fmt7_blastp';
for ($i = 1 ; $i <= $Nparts4BLAST ; $i++) {
  system "cat out.${i}.fmt7_blastp >> out.fmt7_blastp";
}
if (-z 'out.fmt7_blastp') {
  die "problem with step BLAST\n";
}
system 'rm out.*.fmt7_blastp reference_sequences_nonredundant.*.fa BLASTscript.*.sh flag.BLAST.*';
# workflow output : out.fmt7_blastp

PARSEBLAST: # parse out.fmt7_blastp to find aa substitutions that are located
  # in the PDB structure. We assume that 1 SNP can affect multiple chains.
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/preparseBLASToutput.pl") {
  die "problem with step PARSEBLAST preparseBLASToutput.pl\n";
}
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseBLASToutput.pl $PDBseqfile $PDBposfile") {
  die "problem with step PARSEBLAST parseBLASToutput.pl\n";
}
unlink 'BLASTfile_preprocessed';
system 'sort -n -k 1 variants4foldX_unsorted.tab > variants4foldX.tab';
system 'sort -n -k 2 variants_without_structure_info_unsorted.txt > variants_without_structure_info.txt';
unlink 'variants4foldX_unsorted.tab', 'variants_without_structure_info_unsorted.txt';
unlink 'out.fmt7_blastp'; # IF NEEDED OUTCOMMENT
# workflow output : variants4foldX.tab variants_without_structure_info.txt

FOLDX: # run FoldX on structure - mutation pairs mentioned in
  # variants4foldX.tab and extract useful part of output
system 'cut -f 8,9 variants4foldX.tab | sort -u > protvariants4foldX.tab';
$Nlines = `wc -l protvariants4foldX.tab`;
$Nlines = $Nlines * 1; # remove file name from wc output
if ($Nparts4FoldX > $Nlines) { $Nparts4FoldX = 1 }
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/splitvariants4FoldX.pl $rotabase $Nparts4FoldX") {
  die "problem with step FOLDX splitvariants4FoldX.pl\n";
}
symlink "$scriptdir/FoldXoptions", 'FoldXoptions';
symlink "$scriptdir/FoldXcommands", 'FoldXcommands';
symlink "$scriptdir/FoldXcommands4AnalyseComplex",
  'FoldXcommands4AnalyseComplex';
for ($i = 1 ; $i <= $Nparts4FoldX ; $i++) {
  chdir "FOLDXrunfolder_${i}";
  system "qsub -cwd -b y -l mem_limit=$FOLDXMEM $scriptdir/submit2FoldX.pl $FoldX $PDBDIR";
  chdir '..';
}
$DONE = 0;
until ($DONE) {
  $DONE = 1;
  for ($i = 1 ; $i <= $Nparts4FoldX ; $i++) {
    if (not -e "FOLDXrunfolder_${i}/flag.FoldX") { $DONE = 0 }
  }
}
$error = 0;
for ($i = 1 ; $i <= $Nparts4FoldX ; $i++) {
  $e = `wc -l FOLDXrunfolder_${i}/submit2FoldX.pl.e*`;
  $e = $e * 1; # remove file name from wc output
  if ($e > 0) {
    warn "problem with step FOLDX in folder FOLDXrunfolder_${i}\n";
    $error = 1;
  }
}
if ($error) { exit(1) }
unlink 'BuildModel_lastline.tab', 'SequenceDetail_mutlines.tab',
  'AnalyseComplex_meandiff.tab', 'interactingchains.tab';
for ($i = 1 ; $i <= $Nparts4FoldX ; $i++) {
  system "cat FOLDXrunfolder_${i}/BuildModel_lastline.tab >> BuildModel_lastline.tab";
  system "cat FOLDXrunfolder_${i}/SequenceDetail_mutlines.tab >> SequenceDetail_mutlines.tab";
  system "cat FOLDXrunfolder_${i}/AnalyseComplex_meandiff.tab >> AnalyseComplex_meandiff.tab";
  system "cat FOLDXrunfolder_${i}/interactingchains.tab >> interactingchains.tab";
  system "rm -r FOLDXrunfolder_${i}";
}
unlink 'protvariants4foldX.tab';
unlink 'FoldXoptions', 'FoldXcommands', 'FoldXcommands4AnalyseComplex';
# workflow output : BuildModel_lastline.tab SequenceDetail_mutlines.tab
#   AnalyseComplex_meandiff.tab interactingchains.tab
# NOTE : the script submit2FoldX.pl deletes most of the FoldX output
#   because it is very bulky. The user should run FoldX "manually"
#   on the mutants that have been selected for final inspection.

PARSEFOLDX: # put the pieces together to make the final FoldX report
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/makeFoldXreport.pl") {
  die "problem with step PARSEFOLDX makeFoldXreport.pl\n";
}
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/makeFoldXSequenceDetailreport.pl $PDBDIR") {
  die "problem with step PARSEFOLDX makeFoldXSequenceDetailreport.pl\n";
}
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/makeFoldXAnalyseComplexreport.pl") {
  die "problem with step PARSEFOLDX makeFoldXAnalyseComplexreport.pl\n";
}
unlink 'BuildModel_lastline.tab', 'SequenceDetail_mutlines.tab',
  'AnalyseComplex_meandiff.tab', 'interactingchains.tab',
  'variants4foldX.tab', 'SNPAAchange.vcf';
# workflow output : FoldXreport.tab FoldXreport_SequenceDetail.tab
#   FoldXreport_AnalyseComplex.tab  withPDB.vcf

GENE3D: # run the Gene3D pipeline on reference_sequences_nonredundant.fa
  # and make mutated_domains.tab
$Nlines = `wc -l reference_sequences_nonredundant.fa`;
$Nlines = $Nlines * 1; # remove file name from wc output
$Nseq = $Nlines / 2;
if ($Nparts4Gene3D > $Nseq) { $Nparts4Gene3D = 1 }
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/splitreferences4Gene3D.pl $gene3d $Nparts4Gene3D") {
  die "problem with step GENE3D splitreferences4Gene3D.pl\n";
}
for ($i = 1 ; $i <= $Nparts4Gene3D ; $i++) {
  chdir "Gene3Drunfolder_${i}";
  system "qsub -cwd -b y -l mem_limit=$MEM $scriptdir/submit2Gene3D.pl $python $hmmsearch $gene3d";
  chdir '..';
}
$DONE = 0;
until ($DONE) {
  $DONE = 1;
  for ($i = 1 ; $i <= $Nparts4Gene3D ; $i++) {
    if (not -e "Gene3Drunfolder_${i}/flag.Gene3D") { $DONE = 0 }
  }
}
$error = 0;
for ($i = 1 ; $i <= $Nparts4Gene3D ; $i++) {
  $e = `wc -l Gene3Drunfolder_${i}/submit2Gene3D.pl.e*`;
  $e = $e * 1; # remove file name from wc output
  if ($e > 0) {
    warn "problem with GENE3D step submit2Gene3D.pl in folder Gene3Drunfolder_${i}\n";
    $error = 1;
  }
}
if ($error) { exit(1) }
unlink 'seqs.domains', 'no_domain_info.txt';
for ($i = 1 ; $i <= $Nparts4Gene3D ; $i++) {
  system "cat Gene3Drunfolder_${i}/seqs.domains >> seqs.domains";
  system "cat Gene3Drunfolder_${i}/no_domain_info.txt >> no_domain_info.txt";
  system "rm -r Gene3Drunfolder_${i}";
}
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseGene3Doutput.pl") {
  die "problem with GENE3D step parseGene3Doutput.pl\n";
}
unlink 'no_domain_info.txt';
unlink 'seqs.domains'; # IF NEEDED OUTCOMMENT
# workflow output : mutated_domains.tab  variants_without_domain_info.txt

AGADIR: # run TANGO and WALTZ on reference and mutated sequences
  # and parse on-the-fly the per residue scores
$Nlines = `wc -l variants.tab`;
$Nlines = $Nlines * 1; # remove file name from wc output
if ($Nparts4AGADIR > $Nlines) { $Nparts4AGADIR = 1 }
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/splitfiles4AGADIR.pl $scriptdir/Options.txt $Nparts4AGADIR") {
  die "problem with step AGADIR splitfiles4AGADIR.pl\n";
}
for ($i = 1 ; $i <= $Nparts4AGADIR ; $i++) {
  chdir "AGADIRrunfolder_${i}";
  system "qsub -cwd -b y -l mem_limit=$MEM $scriptdir/submit2AGADIR.pl $agadir";
  chdir '..';
}
$DONE = 0;
until ($DONE) {
  $DONE = 1;
  for ($i = 1 ; $i <= $Nparts4AGADIR ; $i++) {
    if (not -e "AGADIRrunfolder_${i}/flag.AGADIR") { $DONE = 0 }
  }
}
$error = 0;
for ($i = 1 ; $i <= $Nparts4AGADIR ; $i++) {
  $e = `wc -l AGADIRrunfolder_${i}/submit2AGADIR.pl.e*`;
  $e = $e * 1; # remove file name from wc output
  if ($e > 0) {
    warn "problem with AGADIR step submit2AGADIR.pl in folder AGADIRrunfolder_${i}\n";
    $error = 1;
  }
}
if ($error) { exit(1) }
unlink 'reference_TANGO.tab', 'reference_WALTZ.tab',
  'reference_AGADIR_summary.tab', 'mutated_TANGO.tab', 'mutated_WALTZ.tab',
  'mutated_AGADIR_summary.tab', 'domain_differences_AGADIR.tab';
for ($i = 1 ; $i <= $Nparts4AGADIR ; $i++) {
  system "cat AGADIRrunfolder_${i}/reference_TANGO.tab >> reference_TANGO.tab";
  system "cat AGADIRrunfolder_${i}/reference_WALTZ.tab >> reference_WALTZ.tab";
  system "cat AGADIRrunfolder_${i}/reference_AGADIR_summary.tab >> reference_AGADIR_summary.tab";
  system "cat AGADIRrunfolder_${i}/mutated_TANGO.tab >> mutated_TANGO.tab";
  system "cat AGADIRrunfolder_${i}/mutated_WALTZ.tab >> mutated_WALTZ.tab";
  system "cat AGADIRrunfolder_${i}/mutated_AGADIR_summary.tab >> mutated_AGADIR_summary.tab";
  system "cat AGADIRrunfolder_${i}/domain_differences_AGADIR.tab >> domain_differences_AGADIR.tab";
  system "rm -r AGADIRrunfolder_${i}";
}
unlink 'reference_sequences.fa';
# workflow output :
#   reference_TANGO.tab reference_WALTZ.tab reference_AGADIR_summary.tab
#   mutated_TANGO.tab mutated_WALTZ.tab mutated_AGADIR_summary.tab
#   domain_differences_AGADIR.tab

PARSEAGADIR: # compare the TANGO and WALTZ output for reference and mutant
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseAGADIR.pl") {
  die "problem with PARSEAGADIR step parseAGADIR.pl\n";
}
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/parseAGADIR4Gene3D.pl") {
  die "problem with PARSEAGADIR step parseAGADIR4Gene3D.pl\n";
}
unlink 'reference_AGADIR_summary.tab', 'mutated_AGADIR_summary.tab';
# workflow output :
#   region_differences_AGADIR.tab
#   domain_regions_differences_AGADIR.tab

TMHMM: # run TMHMM on sequences in reference_sequences_nonredundant.fa
  # and make check list with TMs in protein, domain and variants site
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/findTMs.pl $tmhmm") {
  die "problem with step TMHMM\n";
}
# workflow output : TMHMMoutput.tab TMcheck.tab

SIFT_PROVEAN: # run SIFT and PROVEAN on each SNP
$Nlines = `wc -l reference_sequences_nonredundant.fa`;
$Nlines = $Nlines * 1; # remove file name from wc output
$Nseq = $Nlines / 2;
if ($Nparts4predictors > $Nseq) { $Nparts4predictors = 1 }
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/splitreferences4predictors.pl $Nparts4predictors") {
  die "problem with step SIFT_PROVEAN splitreferences4predictors.pl\n";
}
for ($i = 1 ; $i <= $Nparts4predictors ; $i++) {
  chdir "runfolder4predictors_${i}";
  system "qsub -cwd -b y -l mem_limit=$MEM $scriptdir/runpredictors.pl $siftdir $BLAST4SIFTdir $UniRef90_BLASTDB $PROVEAN";
  chdir '..';
}
$DONE = 0;
until ($DONE) {
  $DONE = 1;
  for ($i = 1 ; $i <= $Nparts4predictors ; $i++) {
    if (not -e "runfolder4predictors_${i}/flag.SIFT") { $DONE = 0 }
  }
}
for ($i = 1 ; $i <= $Nparts4predictors ; $i++) {
  $e = `wc -l runfolder4predictors_${i}/runpredictors.pl.e*`;
  $e = $e * 1; # remove file name from wc output
  if ($e > 0) {
    warn "problem with step SIFT_PROVEAN in folder runfolder4predictors_${i}\n";
    $error = 1;
  }
}
if ($error) { exit(1) }
unlink 'SIFT_PROVEAN.tab', 'SIFT_warnings.txt';
for ($i = 1 ; $i <= $Nparts4predictors ; $i++) {
  system "cat runfolder4predictors_${i}/SIFT_PROVEAN.tab >> SIFT_PROVEAN.tab";
  system "cat runfolder4predictors_${i}/SIFT_warnings.txt >> SIFT_warnings.txt";
  system "cat runfolder4predictors_${i}/SIFT_errors.txt >> SIFT_errors.txt";
  system "rm -r runfolder4predictors_${i}";
}
# workflow output : SIFT_PROVEAN.tab SIFT_warnings.txt SIFT_errors.txt

FINAL: # make final report by merging Gene3D, TANGO/WALTZ, TMHMM output
  # and mutation predictors output
if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/makeSEQANALreport.pl $SwissProtstandard") {
  die "problem with step FINAL makeSEQANALreport.pl\n";
}
if ($SwissProtstandard) {
  if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/countgenes_withSwissProtstandard.pl") {
    die "problem with step FINAL countgenes_withSwissProtstandard.pl\n";
  }
} else {
  if (system "qsub -cwd -b y -sync y -l mem_limit=$MEM $scriptdir/countgenes.pl") {
    die "problem with step FINAL countgenes.pl\n";
  }
}
unlink 'snpEff_genes.txt'; # IF NEEDED OUTCOMMENT
unlink 'variants.tab', 'reference_TANGO.tab', 'mutated_TANGO.tab',
  'reference_WALTZ.tab', 'mutated_WALTZ.tab',
  'domain_differences_AGADIR.tab', 'region_differences_AGADIR.tab',
  'domain_regions_differences_AGADIR.tab', 'mutated_domains.tab',
  'TMcheck.tab', 'SIFT_PROVEAN.tab', 'withPDB.vcf';
system 'rm *.o* *.e*';
open MARK, '>PIPELINE_FINISHED' ; close MARK;
# workflow output :
#   SEQANALreport.tab SNPpipelinereport.vcf finalreport.txt
#   SNPpipelinereport.vcf 
