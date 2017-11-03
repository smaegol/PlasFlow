#!/usr/bin/env python

#######################################################################################
###                                                                                 ###
###     Copyright (C) 2017  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

import os
import sys
import argparse


# parse command line arguments
parser = argparse.ArgumentParser(
    description='PlasFlow - predicting plasmid sequences in metagenomic data using genome signatures. PlasFlow is based on the TensorFlow artificial neural network classifier')

parser.add_argument('--input', dest='inputfile', action='store',
                    help='Input fasta file with sequences to classify (required)', required=True)
parser.add_argument('--output', dest='outputfile', action='store',
                    help='Output file with classification results (required)', required=True)
parser.add_argument('--threshold', dest='threshold', action='store', type=float,
                    help='Threshold for probability filtering (default=0.7)', default=0.7)
parser.add_argument('--labels', dest='labels',
                    action='store', help='Custom labels file')
parser.add_argument('--models', dest='models',
                    action='store', help='Custom models localization')


args = parser.parse_args()

import numpy as np
import pandas as pd
import rpy2
#rpy2.rinterface.set_initoptions(('rpy2', '--no-save', '--no-restore', '--quiet'))
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import re


# srcipt path is required to find the location of models used for classification (script_path/models)
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

if(args.models):
    models_path = args.models
else:
    models_path = script_path + '/models'

r = robjects.r

# import Biostrings package for kmer quantification
biostrings = importr('Biostrings')

# import labels description
if (args.labels):
    labels_df = pd.read_csv(args.labels, sep="\t")
else:
    labels_df = pd.read_csv(
        models_path + '/class_labels_df.tsv', sep="\t")


# number of classes in the labels file - should equal the number of classes in trained model, otherwise an error will be thrown on the later step
no_classes = labels_df.shape[0]


# get input file with sequences to process
inputfile = args.inputfile

print("Importing sequences")
# read data to classify and quanitify kmers
input_data = r.readDNAStringSet(inputfile)
no_sequences = r.length(input_data)
print("Imported ", no_sequences, " sequences")
# get accessions of files
accessions = r.sub("(\S*)\s.*", "\\1", r.names(input_data), perl=True)

# create pandas frame with info about contigs (id, name, length)
pd_accessions = pandas2ri.ri2py(accessions)
pd_accessions = pd.DataFrame(pd_accessions)
pd_accessions.index.name = 'contig_name'
pd_accessions.reset_index(inplace=True)
pd_accessions.columns = ['contig_id', 'contig_name']
lengths = r.width(input_data)
pd_lengths = pandas2ri.ri2py(lengths)
pd_lengths = pd.DataFrame(pd_lengths)
pd_lengths.index.name = 'contig_id'
pd_lengths.reset_index(inplace=True)
pd_lengths.columns = ['contig_id', 'contig_length']
pd_contigs_info = pd.merge(pd_accessions, pd_lengths, on=['contig_id'])
pd_contigs_info


# Create class defining single classifier

class tf_classif:
    """Main classifier."""

    def __init__(self, kmer, hidden):
        """Initialize the class using kmer length and hidden neurons comfiguration."""
        self.kmer = kmer

        if (hidden == "30"):
            self.hidden = [30]
        elif (hidden == "20_20"):
            self.hidden = [20, 20]
        # set locations of models
        if (kmer == 5):
            if(hidden == "30"):
                self.modeldir = models_path+ "/kmer5_split_30_neurons_relu/"
            elif (hidden == "20_20"):
                self.modeldir = models_path+ "/kmer5_split_20_20_neurons_relu/"
            else:
                print("Wrong hidden layers specification. Exiting...")
                sys.exit()
        elif (kmer == 6):
            if(hidden == "30"):
                self.modeldir = models_path+ "/kmer6_split_30_neurons_relu/"
            elif (hidden == "20_20"):
                self.modeldir = models_path+ "/kmer6_split_20_20_neurons_relu/"
            else:
                print("Wrong hidden layers specification. Exiting...")
                sys.exit()
        elif (kmer == 7):
            if(hidden == "30"):
                self.modeldir = models_path+ "/kmer7_split_30_neurons_relu/"
            elif (hidden == "20_20"):
                self.modeldir = models_path+ "/kmer7_split_20_20_neurons_relu/"
            else:
                print("Wrong hidden layers specification. Exiting...")
                sys.exit()
        else:
            print("Wrong kmer number. Exiting...")
            sys.exit()

    def calculate_freq(self, data):
        """Calculate kmer frequencies and perform td-idf transformation."""
        kmer = self.kmer
        print("Calculating kmer frequencies using kmer", kmer)
        kmer_count = r.oligonucleotideFrequency(data, kmer)
        self.no_features = kmer_count.ncol

        print("Transforming kmer frequencies")
        # Tfidf transform data
        from sklearn.feature_extraction.text import TfidfTransformer
        transformer = TfidfTransformer(smooth_idf=False)
        test_tfidf = transformer.fit_transform(kmer_count)
        test_tfidf_nd = test_tfidf.toarray()
        self.testing_data = test_tfidf_nd

    def predict_proba_tf(self, data):
        """Perform actual prediction (with probabilities)."""
        kmer = self.kmer
        if not hasattr(self, 'testing_data'):
            self.calculate_freq(data)

        # import trained tensorflow model
        import tensorflow as tf
        feature_columns = [tf.contrib.layers.real_valued_column(
            "", dimension=self.no_features)]
        classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=self.hidden,
                                                    n_classes=no_classes,
                                                    model_dir=self.modeldir)

        print("Predicting labels using kmer", kmer, " frequencies")

        # predict probabilities
        new_test_proba = classifier.predict_proba(self.testing_data)
        return new_test_proba

    def predict(self, data):
        """Perform actual prediction (Without probabilities)."""
        if not hasattr(self, 'testing_data'):
            self.calculate_freq(data)

        # import trained tensorflow model
        import tensorflow as tf
        feature_columns = [tf.contrib.layers.real_valued_column(
            "", dimension=self.no_features)]
        classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=self.hidden,
                                                    n_classes=no_classes,
                                                    model_dir=self.modeldir)

        # predict classes
        new_test = classifier.predict(self.testing_data)
        return new_test


# class for voting classifier
class TF_Vote_Classifier:
    """Voting classifier class."""

    def __init__(self, clfs, weights=None):
        """Initialize the voting classifier class."""
        self.clfs = clfs
        self.weights = weights

    def predict_proba(self, X):
        """Return average probabilities."""
        self.probas_ = [clf.predict_proba_tf(X) for clf in self.clfs]
        print("Voting...")
        avg = np.average(self.probas_, axis=0, weights=self.weights)

        return avg

    def predict(self, X):
        """Perform actual prediction."""
        self.classes_ = np.asarray([clf.predict(X) for clf in self.clfs])

        if self.weights:
            avg = self.predict_proba_tf(X)

            maj = np.apply_along_axis(lambda x: max(
                enumerate(x), key=operator.itemgetter(1))[0], axis=1, arr=avg)

        else:
            maj = np.asarray([np.argmax(np.bincount(self.classes_[:, c]))
                              for c in range(self.classes_.shape[1])])

        return maj

    def return_individual_probas(self, data):
        """Return probabilities for individual classifiers."""
        if hasattr(self, "probas_"):
            return self.probas_
        else:
            return 0

    def return_individual_classes(self, data):
        """Return classes outputted by each classifier."""
        if hasattr(self, "classes_"):
            return self.classes_
        else:
            return 0


# classifiers used in PlasFlow (2 hidden layers with 20 neurons in each for 5 and 6-mers and 1 hidden layer with 30 neurons for 7-mers)

kmer5_20_20 = tf_classif(5, "20_20")
kmer6_20_20 = tf_classif(6, "20_20")
kmer7_30 = tf_classif(7, "30")

# voting classifier
vote_class = TF_Vote_Classifier(clfs=[kmer5_20_20, kmer6_20_20, kmer7_30])

vote_proba = vote_class.predict_proba(input_data)
vote = vote_class.predict(input_data)


# results pandas dataframe:
pd_n = pd.DataFrame(vote)
# add columns with contig_id
pd_n.index.name = 'contig_id'
pd_n.reset_index(inplace=True)
pd_n.columns = ['contig_id', 'id']

pd_n_proba = pd.DataFrame(vote_proba)
pd_n_proba.columns = labels_df['label']
pd_n_proba.index.name = 'contig_id'
pd_n_proba.reset_index(inplace=True)


# add labels to classification results
results_merged = pd.merge(pd_n, labels_df, on=['id'])
results_merged_proba = pd.merge(results_merged, pd_n_proba, on=['contig_id'])
results_merged_proba_with_names = pd.merge(
    pd_contigs_info, results_merged_proba, on=['contig_id'])


print("Filtering by probability threshold", args.threshold)

for index, row in results_merged_proba_with_names.iterrows():
    label_name = row.label
    taxname = label_name.split(".", 1)[1]

    if row[label_name] < args.threshold:
        plasmids = row[[col for col in results_merged_proba_with_names.columns if re.match(
            r'^plasmid.*', col)]]
        plasmidssum = plasmids.sum()
        chromosomes = row[[col for col in results_merged_proba_with_names.columns if re.match(
            r'^chromosom.*', col)]]
        chromosomessum = chromosomes.sum()
        my_regex = r".*" + re.escape(taxname) + r""
        taxnames = row[[
            col for col in results_merged_proba_with_names.columns if re.match(my_regex, col)]]
        taxnamessum = taxnames.sum()
        if plasmidssum > args.threshold:
            temp = results_merged_proba_with_names.set_value(
                index, 'label', 'plasmid.unclassified')
        elif chromosomessum > args.threshold:
            temp = results_merged_proba_with_names.set_value(
                index, 'label', 'chromosome.unclassified')
        elif taxnamessum > args.threshold:
            temp = results_merged_proba_with_names.set_value(
                index, 'label', 'unclassified.' + taxname)
        else:
            temp = results_merged_proba_with_names.set_value(
                index, 'label', 'unclassified.unclassified')


results_merged_proba_with_names.to_csv(args.outputfile, sep='\t')

taxons = {}
plasmids = {}
taxon_column = []
plasmid_column = []

for index, row in results_merged_proba_with_names.iterrows():
    label_name = row.label
    taxname = label_name.split(".", 1)[1]
    plasmid = label_name.split(".", 1)[0]
    if taxname in taxons.keys():
        taxons[taxname] = taxons[taxname] + 1
    else:
        taxons[taxname] = 1
    if plasmid in plasmids.keys():
        plasmids[plasmid] = plasmids[plasmid] + 1
    else:
        plasmids[plasmid] = 1
    taxon_column.append(taxname)
    plasmid_column.append(plasmid)


plasmids_pd = pd.DataFrame.from_dict(plasmids, orient="index")
taxons_pd = pd.DataFrame.from_dict(taxons, orient="index")

plasmids_pd = plasmids_pd.transpose()
taxons_pd = taxons_pd.transpose()

print("\nResulting plasmid sequences prediction:")
print(plasmids_pd)
print("\nResulting taxonomical assignment:")
print(taxons_pd)

from Bio import SeqIO

print("\nOutputting fasta files with classified sequences")

sequences_dict = SeqIO.index(args.inputfile, "fasta")
plasmid_sequences = []
chromosome_sequences = []
unclassified_sequences = []
for index, row in results_merged_proba_with_names.iterrows():
    label_name = row.label
    contig_name = row.contig_name
    processed_sequence = sequences_dict[contig_name]
    processed_sequence.id = processed_sequence.id + " " + label_name
    if re.match(r'^chromosome.*', label_name):
        chromosome_sequences.append(processed_sequence)
    elif re.match(r'^plasmid.*', label_name):
        plasmid_sequences.append(processed_sequence)
    else:
        unclassified_sequences.append(processed_sequence)


SeqIO.write(chromosome_sequences, args.outputfile +
            "_chromosomes.fasta", "fasta")
SeqIO.write(plasmid_sequences, args.outputfile + "_plasmids.fasta", "fasta")
SeqIO.write(unclassified_sequences, args.outputfile +
            "_unclassified.fasta", "fasta")
