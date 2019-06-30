
import argparse
import pandas as pd
from os.path import join

from .plasflow import main
from .constants import MODELS_PATH


# parse command line arguments
parser = argparse.ArgumentParser(
    description='PlasFlow v.1.1 - predicting plasmid sequences in metagenomic data using genome signatures. PlasFlow is based on the TensorFlow artificial neural network classifier')

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
parser.add_argument('--batch_size', dest='batch_size',
                    action='store', default=25000, help='Batch size for large datasets')
parser.add_argument('--no_chkpt', help='Do not use checkpoints', action='store_true')
parser.add_argument('--min_len', dest='min_length',
                    action='store', default=1000, help='Minimum size for scaffolds')
args = parser.parse_args()


def cli_main():
    labels_df = pd.read_csv(join(MODELS_PATH, 'class_labels_df.tsv'), sep="\t")
    if args.labels:
        labels_df = pd.read_csv(args.labels, sep="\t")
    main(labels_df, args.inputfile, args.outputfile, args.threshold, min_length=args.min_length)


if __name__ == '__main__':
    cli_main()
