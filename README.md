# PlasFlow 1.0
A set of scripts used for prediction of plasmid sequences in metagenomic contigs


## Requirements:
- Python 3.5
- Python packages:
  - Scikit-learn 0.18.1 
  - Numpy
  - Pandas
  - TensorFlow
  - rpy2
- R 3.25 
- R packages:
  - Biostrings


## Installation

PlasFlow can be easily installed as an Anaconda package from my Anaconda channel using:

    conda install plasflow -c smaegol

To exclude the possibility of dependencies conflicts its encouraged to create spearate conda environment for Plasflow
using command: 

    conda create --name plasflow python=3.5

to activate created environment type:

    source activate plasflow

When you decide to finish tour work with PlasFlow, you can simply deactivate current anaconda environment with command:

    source deactivate

Of course, PlasFlow repo can be cloned using 

    git clone https://github.com/smaegol/PlasFlow

but in that case all dependencies have to be installed manually.

## Getting started

PlasFlow is designed to take a metagenomic assembly and identify contigs which may come from plasmids. It outputs several files, from which the most important is a tabular file containing all predictions (specified with `--output` option). 

Options available in PlasFlow include:
* `--input` - specifies input  fasta file with assembly contigs to classify [required]
* `--output` - should be tsv file with the tabular output of classification [required]
* `--threshold` - manually specified threshold for probability filtering (default = 0.7)

To invoke PlasFlow on `test.fasta` dataset simply:

    PlasFlow.py --input test.fasta --output test.fasta.plasflow_predictions.tsv --threshold 0.7

## Output

The most important output of PlasFlow is a tabular file containing all predictions (specified with `--output` option), consiting of several columns including:
<table><tr><td>contig_id</td><td>contig_name</td><td>contig_length</td><td>id</td><td>label</td><td>...</td></<tr></table>
where:

* `contig_id`is an internal id of sequence used for the classification
* `contig_name` is a name of contig used in the classification 
* `contig_length` shows the length of a classified sequence
* `id` is an internal id of a produced label (classification)
* `label` is the actual classification
* `...` represents additional columns showing probabilities of assignment to each possible class

Sequences can be classified to 26 classes including: chromosome.Acidobacteria, chromosome.Actinobacteria, chromosome.Bacteroidetes, chromosome.Chlamydiae, chromosome.Chlorobi, chromosome.Chloroflexi, chromosome.Cyanobacteria, chromosome.DeinococcusThermus, chromosome.Firmicutes, chromosome.Fusobacteria, chromosome.Nitrospirae, chromosome.other, chromosome.Planctomycetes, chromosome.Proteobacteria, chromosome.Spirochaetes, chromosome.Tenericutes, chromosome.Thermotogae, chromosome.Verrucomicrobia, plasmid.Actinobacteria, plasmid.Bacteroidetes, plasmid.Chlamydiae, plasmid.Cyanobacteria, plasmid.DeinococcusThermus, plasmid.Firmicutes, plasmid.Fusobacteria, plasmid.other, plasmid.Proteobacteria, plasmid.Spirochaetes.

If the probability of assignment to given class is lower than threshold (default = 0.7) then the sequence is treated as unclassified.

Additionaly, PlasFlow produces fasta files containing input sequences binned to plasmids, chromosomes and unclassified. 

