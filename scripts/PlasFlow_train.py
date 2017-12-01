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


##################
#
# Training TensorFlow DNN on kmer counts from sequences
#
##################


import argparse

# parse command line arguments
parser = argparse.ArgumentParser(
    description='Train TensorFlow neural network on kmer counts data.')

parser.add_argument('--input', dest='inputfile', action='store',
                    help='input file with kmer counts (required)', required=True)
parser.add_argument('--hidden1', dest='hidden_units1', action='store',
                    help='number of neurons in 1st hidden layer (default=30)', default=30, type=int)
parser.add_argument('--hidden2', dest='hidden_units2', action='store',
                    help='number of neurons in 2nd hidden layer (optional)', type=int)
parser.add_argument('--hidden3', dest='hidden_units3', action='store',
                    help='number of neurons in 3rd hidden layer (optional, --hidden2 required)', type=int)
parser.add_argument('--modeldir', dest='modeldir',
                    action='store', help='path to the dir for storing model')
parser.add_argument('--steps', dest='training_steps', action='store',
                    help='number of training steps (default=1000)', default=1000, type=int)
parser.add_argument('--activation', dest='activation_fun', action='store', help='Activation funcion (default=sigmoid)',
                    choices=['softmax', 'sigmoid', 'relu', 'relu6', 'crelu', 'elu'], default="sigmoid", type=str)

args = parser.parse_args()

import os
import sys
from os.path import basename
import numpy as np
from array import array
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfTransformer
import tensorflow as tf

# check if 3rd layer was set in commandline
if args.hidden_units3 is not None:
    if args.hidden_units2 is not None:
        hidden_units = [args.hidden_units1, args.hidden_units2,
                        args.hidden_units3]  # set neurons for 3 layers
    else:
        print('Error!!\n\n	--hidden2 option was not set!! Terminating...\n')
        sys.exit()
elif args.hidden_units2 is not None:
    # set neurons for 2 layers
    hidden_units = [args.hidden_units1, args.hidden_units2]
else:
    hidden_units = [args.hidden_units1]  # set neurons for 1 layer


# create modeldir path based on input file or command line argument
modeldir = basename(os.path.splitext(args.inputfile)[0]) + "_tf"
if (args.modeldir):
    modeldir = args.modeldir
print('Location of model and checkpoints: ', modeldir)

# read input from file using numpy
print('reading input from', args.inputfile)
training_data = np.recfromcsv(args.inputfile, delimiter='\t', dtype=np.float64)

# retrieve class information and divide into training and testing datasets
features = np.delete(training_data.view(np.float64).reshape(
    training_data.size, -1), training_data.dtype.names.index('plasmid'), axis=1)
number_of_features = features.shape[1]
print('number of features if input file', number_of_features)
print('creating training and testing dataset')
#stratify option is used to maintain the structure of classes between training and testing datasets
training_features, testing_features, training_classes, testing_classes = train_test_split(
    features, training_data['plasmid'], stratify=training_data['plasmid'], random_state=42)

# get number of classes in training data
number_of_classes = np.unique(training_classes).shape[0]
print('number of classes in training data:', number_of_classes)

# perform tfidf (term-frequency times inverse document-frequency) transformation using scikit-learn
print('tfidf transforming data')
transformer = TfidfTransformer(smooth_idf=False)
training_tfidf = transformer.fit_transform(training_features)
testing_tfidf = transformer.fit_transform(testing_features)

# convert transformed values to numpy array - acceptable by TensorFlow
training_tfidf_nd = training_tfidf.toarray()
testing_tfidf_nd = testing_tfidf.toarray()

# convert class labels to int64 (acceptable by TensorFlow)
training_classes = training_classes.astype(np.int64)
testing_classes = testing_classes.astype(np.int64)


#Save training and testing datasets (and classes), tf-idf transformed, to separate files (for further testing)
np.savetxt(args.inputfile + "_" + "training_tfidf_nd_" +
           args.modeldir + ".tsv", training_tfidf_nd, '%s', '\t')
np.savetxt(args.inputfile + "_" + "testing_tfidf_nd_" +
           args.modeldir + ".tsv", testing_tfidf_nd, '%s', '\t')

np.savetxt(args.inputfile + "_" + "training_classes_" +
           args.modeldir + ".tsv", training_classes, '%s', '\t')
np.savetxt(args.inputfile + "_" + "testing_classes_" +
           args.modeldir + ".tsv", testing_classes, '%s', '\t')


# define validation metrics saved in checkpoints - for visualization of learning in TensorBoard
validation_metrics = {"accuracy": tf.contrib.metrics.streaming_accuracy,
                      "precision": tf.contrib.metrics.streaming_precision,
                      "recall": tf.contrib.metrics.streaming_recall,
                      "MAE": tf.contrib.metrics.streaming_mean_absolute_error,
                      "MSE": tf.contrib.metrics.streaming_mean_squared_error}

# define validation monitor - run every 50 steps on testing data
validation_monitor = tf.contrib.learn.monitors.ValidationMonitor(
    testing_tfidf_nd,
    testing_classes,
    every_n_steps=50,
    metrics=validation_metrics)

# define number of feature columns
feature_columns = [tf.contrib.layers.real_valued_column(
    "", dimension=number_of_features)]

# define classifier usinf tf.contrib.learn API
# hidden units describes size of hidden layers (defined in commandline)
# in config checkpoints set to be saved every 1s
if (args.activation_fun == 'sigmoid'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.sigmoid,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))
elif (args.activation_fun == 'relu'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.nn.relu,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))
elif (args.activation_fun == 'relu6'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.nn.relu6,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))
elif (args.activation_fun == 'crelu'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.nn.crelu,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))
elif (args.activation_fun == 'elu'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.nn.elu,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))
elif (args.activation_fun == 'softmax'):
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_units,
                                                n_classes=number_of_classes,
                                                activation_fn=tf.nn.softmax,
                                                model_dir=modeldir,
                                                config=tf.contrib.learn.RunConfig(save_checkpoints_secs=1))

else:
    print("Unknown activation function\n")
    sys.exit()

# define verbosity of classifier - for debugging - INFO
tf.logging.set_verbosity(tf.logging.INFO)

# start model fitting
#training_steps - defined in commandline
# validation monitor allow for visualization in TensorBoard
classifier.fit(x=training_tfidf_nd, y=training_classes,
               steps=args.training_steps, monitors=[validation_monitor])
