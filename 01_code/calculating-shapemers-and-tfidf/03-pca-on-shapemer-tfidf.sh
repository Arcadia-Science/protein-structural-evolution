#Ryan York
#092222
#NOTE: this code needs to be edited to make the directory paths/structure relative to the Github repo
#Raw data and intermediate files also need to taken into account and uploaded somewhere appropriate

############################################
#####Run PCA on a shapemer TFIDF matrix#####
############################################
#Import packages
import sklearn
import numpy
from sklearn.decomposition import PCA

#Set working directory
import os
os.chdir('/home/ubuntu/environment/analysis_files/eukprot/tfidf/')

##Top 10000
#load
dat = numpy.load('10000_components/eukprot_comparative_set_shapemer_tfidf_top10000_components.npy')

#PCA
pca = PCA()
pca.fit_transform(dat.T)

#Compute the loadings
loadings = pd.DataFrame(pca.components_.T)

#Variance explained
print(pca.explained_variance_ratio_)

#Save
numpy.save('10000_components/eukprot_comparative_set_shapemer_tfidf_top10000_components_PCA_var.npy', pca.explained_variance_ratio_)
numpy.save('10000_components/eukprot_comparative_set_shapemer_tfidf_top10000_components_PCA_components.npy', loadings)

##Top 5000
#load
dat = numpy.load('5000_components/eukprot_comparative_set_shapemer_tfidf_top5000_components.npy')

#PCA
pca = PCA(n_components = 100)
pca.fit_transform(dat.T)

#Compute the loadings
loadings = pca.components_.T

#Variance explained
print(pca.explained_variance_ratio_)

#Save
numpy.save('5000_components/eukprot_comparative_set_shapemer_tfidf_top5000_components_PCA_var.npy', pca.explained_variance_ratio_)
numpy.save('5000_components/eukprot_comparative_set_shapemer_tfidf_top5000_components_PCA_components.npy', loadings)
