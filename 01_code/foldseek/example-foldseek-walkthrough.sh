##Using foldseek to compare AlphaFold structures
---

##Context

#In this notebook we will walkthrough how to use foldseek to compare AlphaFold protein structural predictions. The current release of the AlphaFold DB contains >200 million protein predictions from ~1 million phylogenetically diverse species. Here, we will programmitcally access full 'foldomes' from the AlphaFold DB using NCBI taxonomic IDs. We will then go through the foldseek pipeline, first creating databases for the foldomes we want to compare, and then running foldseek and clustering the results.

---
##Programmatically downloading AlphaFold foldomes

#AlphaFold data can be accessed using Google’s command line tool gsutil. AlphaFold proteins are named using the convention proteomes/proteome-tax_id-[TAX ID]-[SHARD ID].tar. We will replace the [TAX ID] portion with the desired NCBI taxonomy IDs, in this case the IDs for Xenopus [8355], mouse [10090], and zebrafish [7955].

Note: since AlphaFold foldomes can be split across multiple directories, we will use * at the end of the file path to download all relevant data.

Also, note: since this is notebook running R, the command line calls are wrapped in the system function. To run in bash, simply use the code inside the parantheses and quotes.

---

#Create directories to save AlphaFold data and foldseek files.

mkdir ~/Desktop/foldseek-walkthrough/
mkdir ~/Desktop/foldseek-walkthrough/foldseek/
mkdir ~/Desktop/foldseek-walkthrough/foldseek/dbs/
mkdir ~/Desktop/foldseek-walkthrough/foldseek/out/
mkdir ~/Desktop/foldseek-walkthrough/data/

#Change to foldseek directory

cd ~/Desktop/foldseek-walkthrough/

#Use `gsutil` to download the three foldomes.
gsutil -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-8355-* data/
gsutil -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-10090-* data/
gsutil -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-7955-* data/

#Untar the foldome directories. foldseek expects just one directory containing all of the relevant .pdb or .cif files (one per protein structure).

cd data/
ls *.tar |xargs -n1 tar -xzf
cd ../

---

##Running foldseek
#Create database with createdb.

foldseek createdb data/ foldseek/dbs/all_foldomesDB

#Compare the protein structures using search.
foldseek search foldseek/dbs/all_foldomesDB foldseek/dbs/all_foldomesDB foldseek/out/all_by_all foldseek/out/tmp -a

#Calculate tm scores (aln2tmscore) and compile into a .tsv (createtsv; useful for downstream analyses).
foldseek aln2tmscore foldseek/dbs/all_foldomesDB foldseek/dbs/all_foldomesDB foldseek/out/all_by_all foldseek/out/all_by_all_tmscore
foldseek createtsv foldseek/dbs/all_foldomesDB foldseek/dbs/all_foldomesDB foldseek/out/all_by_all_tmscore foldseek/out/all_by_all_tmscore.tsv

#Cluster foldseek results. foldseek has three built in clustering methods:
#greedy (--cluster-mode 0)
#blast (--cluster-mode 1)
#cdhit (--cluster-mode 2)
#Here we will use greedy (tends to generate better/more protein clusters than the others, at least among the species tested so far).

#Cluster
foldseek clust foldseek/dbs/all_foldomesDB foldseek/out/all_by_all foldseek/out/clu --cluster-mode 0 --similarity-type 2

#Compile into .tsv
foldseek createtsv foldseek/dbs/all_foldomesDB foldseek/out/all_by_all foldseek/out/clu foldseek/out/clu_greedy.tsv



