library(reticulate)
library(DescTools)
library(lsa)

############################################
#####Load protein tfidf representations#####
############################################
#Initiate reticulate
np <- import("numpy")
pd <- import("pandas")

#Load TFIDF PCs
proteins = np$load("/Users/ryanyork/Desktop/protein_structural_evolution/eukprot_analyses/02_analysis_files/tfidf/eukprot_comparative_set_shapemer_tfidf_top5000_components_PCA_components.npy")

#Load protein names (connect with rownames of TFIDF)
n = np$load("/Users/ryanyork/Desktop/protein_structural_evolution/eukprot_analyses/02_analysis_files/tfidf/eukprot_comparative_set_shapemer_tfidf_top5000_components_proteins.npy")

#Add names to TFIDF matrix
rownames(proteins) = unlist(lapply(strsplit(n, "-"), function(v){v[2]}))

####################################
#####Load and clean uniprot ids#####
####################################
#Set working directory
setwd('~/Desktop/protein_structural_evolution/orthogroups/TCS-OG-Annotations-083122/og_uniprot_accessions/')

#List files()
files = list.files()

#Progress bar
pb <- txtProgressBar(min = 0,      
                     max = length(files), 
                     style = 3,    
                     width = 50,
                     char = ".")

#Loop through and load
accessions = list()
for(i in 1:length(files)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  dat = read.delim(files[i], sep = '\t')
  dat = dat[!is.na(dat$Entry),]
  if(nrow(dat)>1){
    x = data.frame(orthogroup = unlist(lapply(strsplit(files[i], "_"), function(v){v[1]})),
                   uniprot = dat$Entry, 
                   taxon = dat$Taxon)
    accessions[[unlist(lapply(strsplit(files[i], "_"), function(v){v[1]}))]] = x
  }
}

#Save
saveRDS(accessions, '~/Desktop/protein_structural_evolution/eukprot_analyses/02_analysis_files/eukprot_orthogroup_uniprot_ids_091222.RDS')

##########################################################
#####Intersect tfidf representations with orthogroups#####
##########################################################
#Load accessions
accessions = readRDS('~/Desktop/protein_structural_evolution/eukprot_analyses/02_analysis_files/eukprot_orthogroup_uniprot_ids_091222.RDS')

#Combine accessions
accessions = do.call(rbind, accessions)

#Intersect with tfidf
x = intersect(rownames(proteins), accessions$uniprot)
p2 = as.data.frame(proteins[match(x, rownames(proteins)),])
a2 = as.data.frame(accessions[match(x, accessions$uniprot),])

#Split proteins on orthogroup
p2_split = split(p2, a2$orthogroup)

#Progress bar
pb <- txtProgressBar(min = 0,      
                     max = length(p2_split), 
                     style = 3,    
                     width = 50,
                     char = ".")

#Compare intra vs inter-orthogroup cosine distances via bootstrap
reps = 1000
res = list()
for(i in 121:length(p2_split)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  #Get all other orthogroups
  rest = p2[!a2$orthogroup == names(p2_split)[i],]
  
  #Calculate cosine distances of focal orthogroup
  obs = lsa::cosine(as.matrix(p2_split[[i]]))
  obs = mean(unlist(as.data.frame(obs)))

  
  #Calculate mean correlations for n number of random proteins
  perm = c()
  for(j in 1:1000){
    x = rest[sample(1:nrow(rest), nrow(p2_split[[i]]), replace = FALSE),]
    z = lsa::cosine(as.matrix(x))
    perm = c(perm, mean(unlist(as.data.frame(z))))
  }
  
  #Add to list
  l = list(obs, perm)
  names(l) = c('observed', 'permutations')
  res[[names(p2_split)[i]]] = l
}

#Calculate p values
ps = list()
for(i in 1:length(res)){
  ps[[names(res)[i]]] = sum(res[[i]]$observed>res[[i]]$permutations)/length(res[[i]]$permutations)
}
  
plot(unlist(ps), 
     unlist(lapply(p2_split[1:length(ps)], function(x) nrow(x))))

cor(unlist(ps), 
    unlist(lapply(p2_split[1:length(ps)], function(x) nrow(x))))





