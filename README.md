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

See README_OUTPUT for a description of the content of the major output files.

The script uses a database with PDB files created by switchLab. If this
database is updated or you need an alternative, it is necessary to
recreate the BLAST database. Execute the commands :
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


