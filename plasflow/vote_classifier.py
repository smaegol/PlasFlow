
import numpy as np
import operator


class TF_Vote_Classifier:
    """Voting classifier class.

    class for voting classifier
    based on http://sebastianraschka.com/Articles/2014_ensemble_classifier.html
    """

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
