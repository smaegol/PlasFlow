#!/usr/bin/env python
import argparse
import os, sys
from os.path import basename
import numpy as np
import scipy as sp
import pandas as pd
from array import array
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfTransformer
import tensorflow as tf
from rpy2.robjects import pandas2ri
import matplotlib.pyplot as plt



r = robjects.r

#parse command line arguments
parser = argparse.ArgumentParser(description='Classify data using Tensorflow model')

parser.add_argument('--input', dest='inputfile', action='store', help='input fasta file with sequences to classify (required)',required=True)
parser.add_argument('--output', dest='outputfile', action='store', help='Output of classification (required)',required=True)
parser.add_argument('--threshold', dest='threshold', action='store', type=float, help='threshold for probability filtering (default=0.9)',default=0.9)
parser.add_argument('--kmer', dest='kmer', action='store', help='kmer length used for training (default=6)',choices=[3,4,5,6,7,8],default=6,type=int)
parser.add_argument('--modeldir', dest='modeldir', action='store', help='custom modeldir')
parser.add_argument('--labels', dest='labels', action='store', help='custom labels file')
parser.add_argument('--hidden',dest='hidden',action='store',help='Hidden units in custom providd model',type=int)
#parser.add_argument('--outputproba', dest='outputfileproba', action='store', help='Output of classification (probabilities)',required=False)
args = parser.parse_args()


# import Biostrings package for kmer quantification
biostrings = importr('Biostrings')

#import labels description
if (args.labels):
	labels_df = pd.read_csv(args.labels,sep="\t")
else:
	labels_df = pd.read_csv('/home/smaegol/storage//analyses/plasmids/tensorflow/split/test/labels_df2.tsv',sep="\t")


no_classes = labels_df.shape[0]

if (args.modeldir):
	modeldir=args.modeldir
elif (args.kmer==3):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer3_30_neurons_relu/"
elif (args.kmer==4):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer4_30_neurons_relu/"
elif (args.kmer==5):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer5_30_neurons_relu/"
elif (args.kmer==6):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer6_30_neurons_relu/"
elif (args.kmer==7):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer7_30_neurons_relu/"
elif (args.kmer==8):
	modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer8_30_neurons_relu/"

	

#get input file with contigs to process
inputfile = '/home/smaegol/storage/analyses/annotation/ap2ap3sza_metavelvet/AP2AP3SZA_contigs.fasta'
inputfile = args.inputfile

print("Importing sequences")
#read data to classify and quanitify kmers
test_data = r.readDNAStringSet(inputfile)
accessions = r.names(test_data)
print("Calculating kmer frequencies using kmer",args.kmer)
kmer_count = r.oligonucleotideFrequency(test_data,args.kmer)
no_features = kmer_count.ncol

#create pandas frame with info about contigs (id, name, length)
pd_accessions = pandas2ri.ri2py(accessions)
pd_accessions = pd.DataFrame(pd_accessions)
pd_accessions.index.name = 'contig_name'
pd_accessions.reset_index(inplace=True)
pd_accessions.columns = ['contig_id','contig_name']
lengths = r.width(test_data)
pd_lengths = pandas2ri.ri2py(lengths)
pd_lengths = pd.DataFrame(pd_lengths)
pd_lengths.index.name = 'contig_id'
pd_lengths.reset_index(inplace=True)
pd_lengths.columns = ['contig_id','contig_length']
pd_contigs_info = pd.merge(pd_accessions,pd_lengths,on=['contig_id'])
pd_contigs_info


print("Transforming kmer frequencies")
#Tfidf transform data
from sklearn.feature_extraction.text import TfidfTransformer
transformer = TfidfTransformer(smooth_idf=False)
test_tfidf = transformer.fit_transform(kmer_count)
test_tfidf_nd = test_tfidf.toarray()

#import trained tensorflow model
print("Importing model for classification")
#modeldir="/home/smaegol/storage/analyses/plasmids/tensorflow/final_files/fin2/kmer5_30_neurons_relu/"
import tensorflow as tf
sess = tf.InteractiveSession()
feature_columns = [tf.contrib.layers.real_valued_column("", dimension=no_features)]
classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                            hidden_units=[args.hidden],
                                            n_classes=no_classes,
                                            model_dir=modeldir)
#tf.logging.set_verbosity(tf.logging.INFO)



print("Predicting labels")
#predict labels
new_test = classifier.predict(test_tfidf_nd)
#predict probabilities
new_test_proba = classifier.predict_proba(test_tfidf_nd)

#convert classification results to pandas DataFrame
pd_n_proba = pd.DataFrame(new_test_proba)
pd_n = pd.DataFrame(new_test)
#add columns with contig_id
pd_n.index.name = 'contig_id'
pd_n.reset_index(inplace=True)
pd_n_proba.columns=labels_df['label']
pd_n_proba.index.name = 'contig_id'
pd_n_proba.reset_index(inplace=True)
pd_n.columns = ['contig_id','id']

#add labels to classification results
results_merged = pd.merge(pd_n,labels_df,on=['id'])
results_merged_proba = pd.merge(results_merged,pd_n_proba,on=['contig_id'])
results_merged_proba_with_names = pd.merge(pd_contigs_info,results_merged_proba,on=['contig_id'])


print("Filtering by probability threshold",args.threshold)
for index,row in results_merged_proba_with_names.iterrows():
    label_name = row.label
    if row[label_name]<args.threshold:
        temp = results_merged_proba_with_names.set_value(index,'label','unclassified')


results_merged_proba_with_names.to_csv(args.outputfile, sep='\t')

results_merged_proba_with_names.label.value_counts().plot(kind="barh")

plt.tight_layout()
plt.savefig(args.outputfile+".png",dpi=300,)
#np.savetxt(args.outputfile,new_test,delimiter="\t",fmt="%2.3f")
#np.savetxt(args.outputfileproba,new_test_proba,delimiter="\t",fmt="%2.3f")

