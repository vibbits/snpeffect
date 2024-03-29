This repository accompanies the publication "SNPeffect 5.0: Large-scale structural phenotyping of protein coding variants extracted from next-generation sequencing data using AlphaFold models".

# Prerequisites

## GRID system
The current scripts can be executed on a cluster or HPC system using the Sun GRID system. In case, you like to use a different cluster system, you need to change all `qsub` commands in the `masterscript.pl`.

## Provisioning of installed tool for command line use on HPCs (Linux)
In order to kick-start the use of the pipeline, we provided a folder at [https://console.cloud.google.com/storage/browser/snpeffect-5-data](https://storage.googleapis.com/snpeffect-5-data/snpeffect-5.listing.txt) with around 3/4 of the software tools (where the license allows sharing) to be copied to a writeable folder on an HPC system. 
Once you have downloaded the folders to a Linux-based storage present within an HPC system, you will need to get the following software tools from the respective licensees and install them in the same storage location. FoldX from CRG and tmhmm from DTU. 

Once the tools have been downloaded, you need to specify the absolute path of each tool in the `masterscript.pl`.

## Tools

There are a number of bioinformatics tools which are used in this pipeline and need to be installed prior to use on the cluster system. Below, you can find an overview of tools, its versions, the download location, and license. 
For AGADIR (a wrapper containing the protein aggregation predictors TANGO and WALTZ), you can unzip the file `agadirwrapper.zip`.

| tool | version | download location | download location Singularity image | license information |
| --- | --- | ---| --- | --- |
| SnpEff  5.0 | 2020-08-09 | http://pcingola.github.io/SnpEff/download/ | https://depot.galaxyproject.org/singularity/snpeff%3A5.0--0 | SnpEff is open source, released as "LGPLv3". |
| ncbi-blast-2.10.0+/ | 2.10.0 | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ |https://en.wikipedia.org/wiki/Public_domain | 
| FoldX | 3.0 Beta 6 | http://foldxsuite.crg.eu/ | needs to be installed seperately | http://foldxsuite.crg.eu/academic-license-info | 
| AlphaFold Protein structure database | download 15-08-2021 | see creation of blastDB in this repository | https://alphafold.ebi.ac.uk/download | CC-BY 4.0 |
| HMMER | 3.2.1 (June 2018) | http://hmmer.org/ | https://depot.galaxyproject.org/singularity/hmmer%3A3.2.1--he1b5a44_2 | Freely distributed under the BSD open source license |
| Gene3D | n/a | http://gene3d.biochem.ucl.ac.uk/about#summary | preconfigured version see download via https://console.cloud.google.com/storage/browser/snpeffect-5-data| not specified |
| tmhmm | 2.0c | https://services.healthtech.dtu.dk/cgi-bin/sw_request | needs to be installed seperately | dedicated license from DTU |
| perl | > 5.14.2 | https://www.perl.org/get.html | needs to be present on the compute system | GPL or Artistic License |
| sift | 6.2.1 | https://s3.amazonaws.com/sift-public/nsSNV/sift6.2.1.tar.gz |  | more info https://sift.bii.a-star.edu.sg/www/SIFT_help.html |
| ncbi-blast | 2.4.0+ | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ | https://depot.galaxyproject.org/singularity/blast%3A2.10.1--pl526he19e7b1_0 | https://en.wikipedia.org/wiki/Public_domain |
| UniRef90 | 11-Jan-2011 release | https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_releases/release-2011_01/uniref/ |  | Creative Commons Attribution (CC BY 4.0) License to all copyrightable parts of our databases.|
| PROVEAN | v.1.1.5 (May 7, 2014)| http://provean.jcvi.org/downloads/README | https://depot.galaxyproject.org/singularity/provean%3A1.1.5--h87f3376_0 | http://provean.jcvi.org/downloads/LICENSE | 
| CD-HIT | 4.8.1 | https://github.com/weizhongli/cdhit | https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--h2e03b76_4  | GPL v2 | 
| NCBI nr (non-redundant) protein database | Aug 2011 | https://www.ncbi.nlm.nih.gov/home/download/ | to be downloaded via location (7.2 GB) | https://en.wikipedia.org/wiki/Public_domain |

Once the tools have been downloaded or the respective modules have been loaded on an HPC cluster, you need to specify the path of each tool in the `masterscript.pl`.

In addition you will nead to download and install the human genome database for SnpEff (hg19 and/or hg38). You can pre-install databases manually using the `SnpEff download` command (once SnpEff is installed). E.g. to download the human genome database hg38:
```
  java -jar <snpEff>/snpEff.jar download hg38
```

It is highly recommended to repair the PDB structures (energy minimization of the side chains) before you do any modelling with FoldX. You can repair the structures using the RepairPDB function in FoldX. It can be run from the command line:
```
  <foldxdir>/FoldX --command=RepairPDB --pdb=RP.pdb
```
To see more parameters you can visit http://foldxsuite.crg.eu/command/RepairPDB.

Finally, you need to create a BLAST database of the PDB files. In a new empty folder execute the commands:
```
  perl extractseqfromPDB.pl <PDBsdir>
  <blastdir>/makeblastdb -dbtype prot -in PDBsequences.fa -out PDB
```
If this database is updated or you want to use an alternative PDB database, you can follow the same steps to recreate the BLAST database.

# Pipeline testing

After installation the user can run the input VCF file from the SHP-77 carcinoma cell line as a test case to verify correct installment of the software. The input file, intermediate and output files of the test case can be found in the different folders of this repository (named after each step of the pipeline). 

# Perform analysis

To perform an analysis you must make a folder/directory, put in that folder
the VCF file, rename it in.vcf then with that folder as working directory (specified in the `masterscript.pl`)
execute the command :
```
  qsub -cwd -b y <scriptdir>/masterscript.pl
```
It is important to use the same annotated genome as was used to generate
the VCF file. As default the pipeline uses hg38. If you used something else, you will need to edit
the file `masterscript.pl`.

The script also uses a NP_* NM_* RefSeq correspondence table and a table with
UniProt standard (preferred) transcripts. See respectively the scripts
`makeNP_NMtable.pl` and `makeNM_ACtable.pl` in case you want to update them.

NOTE :
You can restrict the output and the CPU time used by analyzing only one
alternative transcript per mutated gene, using the UniProt standards. The
current version of `masterscript.pl` proposes this by default. Note however
that you loose some data since some genes are not in the list and
some mutations are not in the standard transcript.

# OUTPUT files
  
## FoldX reports
  
The files FoldXreport*.tab contain the final report of predictions made
using FoldX and PDB files with the structure of the protein or a closely
related protein. Each file contains the following columns :
- variant number (1 col). There are missing numbers since variants that
  cause an SNP outside a known structure are lacking.
- position SNP in reference genome (1 col)
- gene and transcript identifiers from the reference genome (3 cols)
- position SNP in protein (1 col)
- reference PDB sequence ID (1 col)

The file FoldXreport.tab contains one row for each variant and has
besides the common information:
- FoldX mutation instruction (1 col)
- the labels of the polypeptide chains in the model (1 col)
- the labels of those polypeptide chains that are reported to interact
  energenetically (1 col)
- the change in interaction energy caused by the mutation (1 col)
- the output of the FoldX BuildModel command, which predicts the effect
  of the mutation on the 3D structure (23 cols).
- part of the output of the FoldX SequenceDetail command, which describes the
  environment of the amino acid (e.g. main and side chain burial). (4 cols)
- the pLDDT score from AlphaFold structures (1 col)

The file FoldXreport_AnalyseComplex.tab contains for those proteins for
which the mutation can affect the interaction between chains (not useful for AlphaFold models):
- the labels of each pair of chains (2 cols)
- a prediction of the effect of the mutation (5 cols)

## SEQANAL reports
  
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

## SNPpipeline report

The file SNPpipelinereport.vcf is a VCF file based on the output of
SnpEff, with only the records selected for further analysis and
extra items added to the INFO field (these items provide a selection
of the output of the analysis tools used by the pipeline).

## Other files

The scripts in this folder, as run using the pipeline defined in
`masterscript.pl`, produce a lot of output files. Some of the intermediate
output is deleted, for the sake of saving disk space; you can, if neeeded,
outcomment unlink commands in `masterscript.pl`.

The masterscript routinely preserves the following intermediate and
supplementary results :

- in.vcf : the origninal input VCF file with the variants
- SnpEff_errors.txt : records in in.vcf that could not be parsed
  by SnpEff because of error
- reference_sequences_nonredundant.fa : the protein sequences from the
  reference genome in which an aa changing SNF was observed
- mutated_sequences.fa : the mutated protein sequences
- sequence_identities.tab : a correspondence table, the first column are
  sequences in reference_sequences_nonredundant.fa, all columns together
  point to sequences in mutated_sequences.fa
- variants_not_investigated.txt : reason why variants were not investigated
  (too long for AGADIR or no standard protein in UniProt)
- variants_without_structure_info.txt : reasons why some variants could not
  be submitted to FoldX
- variants_without_domain_info.txt : reasons why domain with mutation
  info could not be reported
- TMHMMoutput.tab : TMHMM output
- SIFT_warnings.txt : warnings issued by SIFT, mainly about amino acids
  in the reference sequence that seem poorly conserved
- SIFT_errors.txt : errors issued by SIFT, reason why SIFT prediction could
  not be reported
- finalreport.txt : a final count of how many mutated genes were reported
- PIPELINE_FINISHED : a timestamp file to indicate when the pipeline finished

