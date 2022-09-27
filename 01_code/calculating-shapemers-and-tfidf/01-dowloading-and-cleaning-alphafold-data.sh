#Ryan York
#092222
#NOTE: this code needs to be edited to make the directory paths/structure relative to the Github repo
#Raw data and intermediate files also need to taken into account and uploaded somewhere appropriate

###############################
#####Downloading proteomes#####
###############################
#Download taxonomy ids from s3
aws s3 cp s3://arcadia-protein-evolution/protein_lists/eukprot_species_lists/the-comparative-set_tax_ids.txt ~/environment/
aws s3 cp s3://arcadia-protein-evolution/protein_lists/eukprot_species_lists/archaeplastidae_ncbi_tax_ids.txt ~/environment/
  
#Create directories
mkdir data/proteomes
mkdir data/proteomes/the_comparative_set/
mkdir data/proteomes/archaeplastidae/

##Download and clean labels, etc. in R
#Start R
R

#Load ncbi taxon ids (collected from NCBI using species names at https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi)
tax = read.delim('~/environment/protein_lists/the-comparative-set_tax_ids.txt')

#Get ids
ids = as.numeric(na.omit(tax$taxid))

#Loop through and download via command line
for(i in 1:length(ids)){
  print(paste(i, 'out of', length(ids)))
  system(paste('bash gsutil -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-', 
               ids[i], 
               '-* ~/environment/data/proteomes/the_comparative_set/', sep = ''))}

#Load ncbi taxon ids (collected from NCBI using species names at https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi)
tax = read.delim('~/environment/protein_lists/archaeplastidae_ncbi_tax_ids.txt')

#Get ids
ids = as.numeric(na.omit(tax$taxid))

#Loop through and download via command line
for(i in 1:length(ids)){
  print(paste(i, 'out of', length(ids)))
  system(paste('bash gsutil -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-', 
               ids[i], 
               '-* ~/environment/data/proteomes/archaeplastidae/', sep = ''))}
               
               
#Loop through and untar
cd ~/environment/data/proteomes/the_comparative_set/
for f in *.tar; do destpath="$(echo "$f"| sed 's/.tar//')"; mkdir -p "$destpath"; tar xf "$f" -C "$destpath"; done

cd ~/environment/data/proteomes/archaeplastidae/
for f in *.tar; do destpath="$(echo "$f"| sed 's/.tar//')"; mkdir -p "$destpath"; tar xf "$f" -C "$destpath"; done

#Transfer to s3
aws s3 sync ~/environment/data/proteomes/ s3://arcadia-protein-evolution/proteomes/eukprot_proteomes/
