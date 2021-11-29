This repository accompanies the publication "SNPeffect 5.0: Large-scale structural phenotyping of protein coding variants extracted from next-generation sequencing data using AlphaFold models".

# Prerequisites

## GRID system
The current scripts can be executed on a cluster using the Sun GRID system. In case, you like to use a different cluster system, you need to change all `qsub` commands in the `masterscript.pl`.

## Tools

There are a number of bioinformatics tools which are used in this pipeline and need to be installed prior to use on the cluster system. Below, you can find an overview of tools, its versions, the download location, and license.

| tool | version | download location | license information |
| --- | --- | ---| --- | 
| SnpEff  5.0 | 2020-08-09 | http://pcingola.github.io/SnpEff/download/| SnpEff is open source, released as "LGPLv3". |
| ncbi-blast-2.10.0+/ | 2.10.0 | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ |https://en.wikipedia.org/wiki/Public_domain | 
| FoldX | 3.0 Beta 6 | http://foldxsuite.crg.eu/ | http://foldxsuite.crg.eu/academic-license-info | 
| AlphaFold Protein structure database | download 15-08-2021 | https://alphafold.ebi.ac.uk/download | CC-BY 4.0 |
| HMMER | 3.2.1 (June 2018) | http://hmmer.org/ | Freely distributed under the BSD open source license |
| Gene3D | n/a | http://gene3d.biochem.ucl.ac.uk/about#summary | not specified |
| Tango | n/a | http://tango.crg.es/ | https://switchlab.netlify.app/contact/ |
| Waltz | n/a | https://waltz.switchlab.org/ | https://switchlab.netlify.app/contact/ |
| tmhmm | 2.0c | https://services.healthtech.dtu.dk/cgi-bin/sw_request | dedicated license from DTU |
| PolyPhen | 2.2.2 | http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads | free for academic instruction and research use only |
| PDB/DSSP structural databases snapshot (38G) | 2.2.3 | http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads | free for academic instruction and research use only |
| perl | > 5.14.2 | https://www.perl.org/get.html | GPL or Artistic License |
| UniRef100 | 14-Dec-2011 | https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_releases/release-2011_01/uniref/| Creative Commons Attribution (CC BY 4.0) License to all copyrightable parts of our databases.|
| sift | 6.2.1 | https://s3.amazonaws.com/sift-public/nsSNV/sift6.2.1.tar.gz | more info https://sift.bii.a-star.edu.sg/www/SIFT_help.html |
| ncbi-blast | 2.4.0+ | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ | https://en.wikipedia.org/wiki/Public_domain |
| UniRef90 | 11-Jan-2011 release | https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_releases/release-2011_01/uniref/ | Creative Commons Attribution (CC BY 4.0) License to all copyrightable parts of our databases.|
| PROVEAN | v.1.1.5 (May 7, 2014)| http://provean.jcvi.org/downloads/README | http://provean.jcvi.org/downloads/LICENSE | 
| CD-HIT | 4.8.1 | https://github.com/weizhongli/cdhit | GPL v2 | 
| NCBI nr (non-redundant) protein database | Aug 2011 | https://www.ncbi.nlm.nih.gov/home/download/ | https://en.wikipedia.org/wiki/Public_domain |

# Perform analysis

To perform an analysis you must make a folder/directory, put in that folder
the VCF file, rename it in.vcf then with that folder as working directory
execute the command :
  qsub -cwd -b y <scriptdir>/masterscript.pl
It is important to use the same annotated genome as was used to generate
the VCF file. If you used something else than hg19 you will need to edit
the file masterscript.pl.

The script uses a database with PDB files from the AlphaFold database. If this
database is updated or you need an alternative, it is necessary to
recreate the BLAST database. Execute the commands:
  extractseqfromPDB.pl
  makeblastdb -dbtype prot -in PDBsequences.fa -out PDB

The script also uses a NP_* NM_* RefSeq correspondence table and a table with
UniProt standard (preferred) transcripts. See respectively the scripts
makeNP_NMtable.pl and makeNM_ACtable.pl for how to make them.

NOTE :
1. You can restrict the output and the CPU time used by analyzing only one
   alternative transcript per mutated gene, using the UniProt standards. The
   current version of masterscript.pl proposes this by default. Note however
   that you loose quite some data since some genes are not in the list and
   some mutations are not in the standard transcript.
2. You can perform a PolyPhen analysis. PolyPhen takes as input a SNP and
   predicts the structural effect on a standard protein product of the gene.
   This is however only possible for the hg19 genome and is quite CPU
   intensive. The current version of masterscript.pl does not do this by
   default.
3. The lists with standard transcripts in the UniProt and PolyPhen databases
   do not coincide strictly. Therefore it is not recommended to do both
   1 and 2.

# OUTPUT files

The files FoldXreport*.tab contain the final report of predictions made
using FoldX and PDB files with the structure of the protein or a closely
related protein. Each file contains the following columns :
- variant number (1 col). There are missing numbers since variants that
  cause an SNP outside a known structure are lacking.
- position SNP in reference genome (1 col)
- gene and transcript identifiers from the reference genome (3 cols)
- position SNP in protein (1 col)
- reference PDB sequence ID (1 col)
- FoldX mutation instruction (1 col)

The file FoldXreport.tab contains one row for each variant and has
besides the common information:
- the labels of the polypeptide chains in the model (1 col)
- the labels of those polypeptide chains that are reported to interact
  energenetically (1 col)
- the change in interaction energy caused by the mutation (1 col)
- the output of the FoldX BuildModel command, which predicts the effect
  of the mutation on the 3D structure (23 cols).

The file FoldXreport_compact.tab is a compacted version of FoldXreport.tab.
Lines that have the same PDB sequence ID and same FoldX mutation instruction
(and hence the same FoldX output) have been fused, the contents of the
noncommon fields have been replaced by comma-separated lists.

The file FoldXreport_SequenceDetail.tab contains for each variant as many
rows as there are polypeptide chains affected by the mutation. The FoldX
mutation instruction is lacking and you have instead :
- the amino acid (1 col)
- the label of the polypeptide chain and the position of the amino acid in
  this chain (2 cols)
- the output of the FoldX SequenceDetail command, which describes the
  environment of the amino acid (33 cols)

The file FoldXreport_AnalyseComplex.tab contains for those proteins for
which the mutation can affect the interaction between chains:
- the labels of each pair of chains (2 cols)
- a prediction of the effect of the mutation (5 cols)

========================================================================

The file SEQANALreport.tab contains the output of a series of software tools
that can be run on the sequence without knowledge of the structure. It contains
a row for each SNP that causes an amino acid change. It contains the following
columns :
- variant number (1 col)
- position SNP in reference genome and nucleotide change (2 cols)
- gene and transcript identifiers from the reference genome (3 cols)
- reference amino acid, position in protein and mutated amino acid (3 cols)
- information about CATH protein domain in which the mutation occurred,
  as found by the Gene3D pipeline (4 cols)
- information about transmembrane domains, as found by TMHMM (3 cols)
- information based on TANGO and WALTZ output (36 cols)
- SIFT prediction of effect mutation (4 cols)
- PROVEAN prediction of effect mutation (3 cols)
- (optionally) the accession number of the UniProt standard sequence,
  preceded by a ~ if the variant is for another transcript of the same
  gene (1 col)

The file SEQANALreport_withPolyPhen.tab contains only variants for
which PolyPhen found a standard transcript. It contains the same rows
as SEQANALreport.tab, with as extra :
- PolyPhen prediction of effect mutation (7 cols)
- UniProt ID of protein sequence used by PolyPhen and position mutation
  in this sequence (2 cols)

SEQANALreport_compact.tab contains a subset of the data in SEQANALreport.tab.
It contains one row for each chromosome location and amino acid change.

========================================================================

The file SNPpipelinereport.vcf (or SNPpipelinereport_withPolyPhen.vcf
when PolyPhen has been run) is a VCF file based on the output of
SnpEff, with only the records selected for further analysis and
extra items added to the INFO field (these items provide a selection
of the output of the analysis tools used by the pipeline).

========================================================================
  
The scripts in this folder, as run using the pipeline defined in
masterscript.pl, produce a lot of output files. Some of the intermediate
output is deleted, for the sake of saving disk space; you can, if neeeded,
outcomment unlink commands in masterscript.pl.

The masterscript routinely preserves the following intermediate and
supplementary results :

in.vcf : the origninal input VCF file with the variants
SnpEff_errors.txt : records in in.vcf that could not be parsed
  by SnpEff because of error
reference_sequences_nonredundant.fa : the protein sequences from the
  reference genome in which an aa changing SNF was observed
mutated_sequences.fa : the mutated protein sequences
sequence_identities.tab : a correspondence table, the first column are
  sequences in reference_sequences_nonredundant.fa, all columns together
  point to sequences in mutated_sequences.fa
variants_not_investigated.txt : reason why variants were not investigated
  (too long for AGADIR or no standard protein in UniProt)
variants_without_structure_info.txt : reasons why some variants could not
  be submitted to FoldX
variants_without_domain_info.txt : reasons why domain with mutation
  info could not be reported
variants_without_PolyPhen_standard.txt : reasons why a UniProt standard
  (preferred) protein based on PolyPhen could not be reported
TMHMMoutput.tab : TMHMM output
SIFT_warnings.txt : warnings issued by SIFT, mainly about amino acids
  in the reference sequence that seem poorly conserved
SIFT_errors.txt : errors issued by SIFT, reason why SIFT prediction could
  not be reported
PolyPhen.predictions.tab : PolyPhen output
finalreport.txt : a final count of how many mutated genes were reported
PIPELINE_FINISHED : a timestamp file to indicate when the pipeline finished

