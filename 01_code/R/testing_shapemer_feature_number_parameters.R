options(stringsAsFactors=F);
library(superml)
library(umap)
library(igraph)
library(phangorn)
library(reticulate)
library(DescTools)
library(lsa)
source("~/Documents/Research/Rfiles/TREBLE_walkthrough_functions_051920.R")

##################################################################################################################
#####Get structural kmers via python (following code from github.com/TurtleTools/alphafold-structural-space/)#####
##################################################################################################################
# #Set autoreload
# %load_ext autoreload
# %autoreload 2
# 
# #Source functions
# from src import make_data
# from pathlib import Path
# 
# #Set path to folder with foldomes
# DATA_FOLDER = Path("/Users/ryanyork/Desktop/alphafold_proteomes/proteomes/")
# 
# #Calculate shapemers for all proteins
# make_data.get_AF_shapemers(DATA_FOLDER)
# 
# #Load the corpus file(s)
# corpus_files = DATA_FOLDER.glob("*_ids_corpus_resolution_4_6*.txt")
# keys_corpus = (line.strip().split("\t") for line in itertools.chain.from_iterable((open(file) for file in corpus_files)))
# keys, corpus = itertools.tee(keys_corpus)
# keys = [k[0] for k in keys]
# corpus = (k[1] for k in corpus)
# 
# #Calculate TFIDF matrix
# print(f"Getting TFIDF matrix for {len(keys)} proteins...")
# vectorizer = TfidfVectorizer(min_df=2)
# tfidf_matrix = vectorizer.fit_transform(corpus)
# 
# #Save TFIDF matrix
# with open(DATA_FOLDER / "shapemer_embedding_TFIDF_matrix_6species_071322.pkl", "wb") as f:
#   pickle.dump((tfidf_matrix), f)

#############################
#####Load results into R#####
#############################
#Load raw shapemers
shapemers = read.delim('~/Desktop/protein_structural_evolution/shapemer_walkthrough/00_data/AF_ids_corpus_resolution_4_6_threshold_50.txt', 
                       sep = ' ',
                       header = FALSE)

#Remove rows without protein names
shapemers = shapemers[lapply(strsplit(shapemers[,1], "\t"), function(x) length(x))>1,]

#Add rownames
rownames(shapemers) = unlist(lapply(strsplit(shapemers[,1], "\t"), function(v){v[1]}))

#Fix column 1
shapemers[,1] = unlist(lapply(strsplit(shapemers[,1], "\t"), function(v){v[2]}))

#Convert to strings and split
s = list()
for(i in 1:nrow(shapemers)){
  s[[rownames(shapemers)[i]]] = paste(shapemers[i,], collapse = " ")
}

#Save
saveRDS(s, '~/Desktop/protein_structural_evolution/shapemer_walkthrough/00_data/human_shapemer_corpus_all_proteins_081122.RDS')

########################################################################################
#####Sweep through n features and select optimal based on downstream protein space)#####
########################################################################################
#Load corpus
s = readRDS('~/Desktop/protein_structural_evolution/shapemer_walkthrough/00_data/human_shapemer_corpus_all_proteins_081122.RDS')

#Set up feature sizes to sweep
toSweep = seq(5, 1000, 5)

#Set up random proteins to test
x = s[sample(seq(1, length(s), 1), 1000)]

#Sweep
res = list()
for(i in 1:length(toSweep)){
  
  print(paste('n features =', toSweep[i], i, 'out of', length(toSweep)))
  
  #Initialise the class
  tfv = TfIdfVectorizer$new(min_df = 0.1,
                            max_features = toSweep[i],
                            remove_stopwords = FALSE)
  
  #Generate the matrix
  #tf_mat = tfv$fit_transform(s[1:2000])
  tf_mat = tfv$fit_transform(x)
  
  #Save
  res[[as.character(toSweep[i])]] = tf_mat
}

#Save
saveRDS(res, '~/Desktop/protein_structural_evolution/shapemer_walkthrough/00_data/human_tfidf_min_df_sweeping_results_081522.RDS')

#Calculate trees for each tfidf
trees = list()
for(i in 1:length(res)){
  print(paste(i, 'out of', length(res)))
  d = dist(res[[i]])
  tr = nj(d)
  trees[[i]] = tr
}

#Compare trees made from tfidfs
dists = c()
for(i in 1:(length(trees)-1)){
  dists = c(dists, RF.dist(trees[[i]], trees[[i+1]]))
}

#Plot
plot(dists, 
     pch = 20,
     col = alpha('gray50', 0.5),
     ylab = 'Robinson-Foulds distance',
     xlab = 'n features',
     xaxt = 'n',
     bty = 'n',
     cex = 1.5,
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylim = c(0,2000))
axis(1, seq(1, length(dists), 20), names(res)[seq(1, length(dists), 20)], cex.axis = 1.5)
loess_fit <- loess(dists ~ seq(1, length(dists), 1))
lines(seq(1, length(dists), 1), 
      predict(loess_fit), 
      col = "black", 
      lwd = 1.5,
      lty = 'dashed')

