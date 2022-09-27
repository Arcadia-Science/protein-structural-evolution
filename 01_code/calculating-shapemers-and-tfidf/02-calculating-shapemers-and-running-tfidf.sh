#Ryan York
#092222
#NOTE: this code needs to be edited to make the directory paths/structure relative to the Github repo
#Raw data and intermediate files also need to taken into account and uploaded somewhere appropriate

#######################################################
#####Calculate shapemers and generate TFIDF matrix#####
#######################################################
#Point to working directory with supporting functions
import os
os.chdir('00_utils/alphafold-structural-space/')

#Import packages
from src import make_data
from pathlib import Path
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import StandardScaler
import itertools
from sklearn.decomposition import NMF
import openTSNE
import scipy.sparse
import pandas as pd
import numpy
import pickle

#Point to data directory (here we are using the Alphafold predictions of the human proteome). 
#Note, this directory is expected to be name in the convention used in the alphafold releases (e.g. UP000....).

DATA_FOLDER = Path('/home/ubuntu/environment/data/proteomes/the_comparative_set/')

#Calculate shapemers for all proteins
make_data.get_AF_shapemers(DATA_FOLDER)

#Load the data
corpus_files = DATA_FOLDER.glob("*_ids_corpus_resolution_4_6*.txt")
keys_corpus = (line.strip().split("\t") for line in itertools.chain.from_iterable((open(file) for file in corpus_files)))
keys, corpus = itertools.tee(keys_corpus)
keys = [k[0] for k in keys]
corpus = (k[1] for k in corpus)

#Calculate TFIDF matrix (all)
print(f"Getting TFIDF matrix for {len(keys)} proteins...")
vectorizer = TfidfVectorizer(min_df=2)
tfidf_matrix = vectorizer.fit_transform(corpus)

#Calculate mean
m = tfidf_matrix.mean(axis = 0)

#Save
numpy.save("/home/ubuntu/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_component_means.npy", m)

#Save
with open(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_all_component.pkl", "wb") as f:
    pickle.dump((tfidf_matrix), f)

#Copy to s3
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_all_component.pkl s3://arcadia-protein-evolution/shapemers/eukprot/

##Calculate TFIDF matrix (selecting n features)
##250
print(f"Getting TFIDF matrix for {len(keys)} proteins...")
vectorizer = TfidfVectorizer(min_df=2, max_features = 250)
tfidf_matrix = vectorizer.fit_transform(corpus)

#Save
mat = tfidf_matrix.toarray()
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top250_components.npy", mat)
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top250_components_proteins.npy", keys)

#Copy to s3
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top250_components.npy s3://arcadia-protein-evolution/shapemers/eukprot/
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top250_components_proteins.npy s3://arcadia-protein-evolution/shapemers/eukprot/

##5000
print(f"Getting TFIDF matrix for {len(keys)} proteins...")
vectorizer = TfidfVectorizer(min_df=2, max_features = 5000)
tfidf_matrix = vectorizer.fit_transform(corpus)

#Save
mat = tfidf_matrix.toarray()
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top5000_components.npy", mat)
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top5000_components_proteins.npy", keys)

#Copy to s3
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top5000_components.npy s3://arcadia-protein-evolution/shapemers/eukprot/
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top5000_components_proteins.npy s3://arcadia-protein-evolution/shapemers/eukprot/

##10000
print(f"Getting TFIDF matrix for {len(keys)} proteins...")
vectorizer = TfidfVectorizer(min_df=2, max_features = 10000)
tfidf_matrix = vectorizer.fit_transform(corpus)

#Save
mat = tfidf_matrix.toarray()
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top10000_components.npy", mat)
numpy.save(DATA_FOLDER / "eukprot_comparative_set_shapemer_tfidf_top10000_components_proteins.npy", keys)

#Copy to s3
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top10000_components.npy s3://arcadia-protein-evolution/shapemers/eukprot/
aws s3 cp ~/environment/data/proteomes/the_comparative_set/eukprot_comparative_set_shapemer_tfidf_top10000_components_proteins.npy s3://arcadia-protein-evolution/shapemers/eukprot/

