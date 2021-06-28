
import os.path
import tensorflow as tf
import pandas as pd
import numpy as np

from os.path import join

from .exceptions import UnavailableLayerSpecError, UnavailableKmerSpecError
from .constants import MODELS_PATH
from .utils import get_kmer_counts


class Tf_Classif:
    """Main classifier.

    Create class defining single classifier
    """

    def __init__(self, kmer, hidden):
        """Initialize the class using kmer length and hidden neurons comfiguration."""
        if kmer not in (5, 6, 7):
            raise UnavailableKmerSpecError('k-mer length ' + str(kmer) + ' not available')
        if hidden not in ('30', '20_20'):
            raise UnavailableLayerSpecError('Layer ' + hidden + ' not available')

        self.kmer = kmer
        self.modeldir = join(MODELS_PATH, 'kmer' + str(kmer) + '_split_' + hidden + '_neurons_relu')
        self.hidden = {
            '30': [30],
            '20_20': [20, 20]
        }[hidden]

    def calculate_freq(self, seqs, input_data_path, no_chkpt=True):
        """Calculate kmer frequencies and perform td-idf transformation."""
        kmer = self.kmer
        file_name = str(input_data_path) + "_kmer_" + str(kmer) + '_freqs.npy'
        # Try to load previously saved frequncies (TF-IDF transformed)
        if (not no_chkpt) and os.path.isfile(file_name):
            test_tfidf_nd = np.load(file_name)
            self.num_features = test_tfidf_nd.shape[1]
            self.testing_data = test_tfidf_nd
            print("Succesfully read previously calculated kmer frequencies for kmer", kmer)
        # if previous calculations are not available - calculate frequencies
        else:
            print("Calculating kmer frequencies using kmer", kmer)
            tbl = {}
            for seq in seqs:
                tbl[seq.id] = get_kmer_counts(seq.seq, self.kmer)
            kmer_count = pd.DataFrame.from_dict(tbl, orient='index').fillna(0).values

            self.num_features = kmer_count.shape[1]

            print("Transforming kmer frequencies")
            # Tfidf transform data
            from sklearn.feature_extraction.text import TfidfTransformer
            transformer = TfidfTransformer(smooth_idf=False)
            test_tfidf = transformer.fit_transform(kmer_count)
            test_tfidf_nd = test_tfidf.toarray()
            self.testing_data = test_tfidf_nd
            print("Finished transforming, saving transformed values")
            np.save(file_name, test_tfidf_nd)

    def predict_proba_tf(self, seqs, labels, data):
        """Perform actual prediction (with probabilities)."""
        kmer = self.kmer
        if not hasattr(self, 'testing_data'):
            self.calculate_freq(seqs, data)

        # import trained tensorflow model
        feature_columns = [tf.contrib.layers.real_valued_column("", dimension=self.num_features)]
        classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=self.hidden,
                                                    n_classes=labels.shape[0],
                                                    model_dir=self.modeldir)

        print("Predicting labels using kmer", kmer, " frequencies")

        # predict probabilities
        new_test_proba = classifier.predict_proba(self.testing_data)
        return new_test_proba

    def predict(self, seqs, labels, data):
        """Perform actual prediction (Without probabilities)."""
        if not hasattr(self, 'testing_data'):
            self.calculate_freq(seqs, data)

        # import trained tensorflow model
        feature_columns = [tf.contrib.layers.real_valued_column("", dimension=self.num_features)]
        classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=self.hidden,
                                                    n_classes=labels.shape[0],
                                                    model_dir=self.modeldir)

        # predict classes
        new_test = classifier.predict(self.testing_data)
        return new_test
