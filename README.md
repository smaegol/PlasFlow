# PlasFlow 1.0

PlasFlow is a set of scripts used for prediction of plasmid sequences in metagenomic contigs. It relies on the neural network models trained on full genome and plasmid sequences and is able to differentiate between plasmids and chromosomes with accuracy reaching 96%. It outperforms other available solutions for plasmids recovery from metagenomes and incorporates the thresholding which allows for exclusion of incertain predictions.


# Table of contents

  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Conda-based](https://github.com/smaegol/PlasFlow#conda-based---recommended)
    - [Pip installer](https://github.com/smaegol/PlasFlow#pip-installer)
    - [Manual installation](https://github.com/smaegol/PlasFlow#manual-installation)
  - [Getting started](https://github.com/smaegol/PlasFlow#getting-started)
  - [Output](https://github.com/smaegol/PlasFlow#output)
  - [Test dataset](https://github.com/smaegol/PlasFlow#test-dataset)
  - [Detailed information](https://github.com/smaegol/PlasFlow#detailed-information)

## Requirements:

- Python 3.5
- Python packages:

  - Scikit-learn 0.18.1
  - Numpy
  - Pandas
  - [TensorFlow 0.10.0](https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.10.0rc0-cp35-cp35m-linux_x86_64.whl)
  - rpy2
  - scipy

- R 3.25

- R packages:

  - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

## Installation

##### Conda-based - recommended

Conda is recommended option for installation as it properly manage all dependencies and allows installation without messing with other packages installed.
Conda can be used both as the [Anaconda](https://www.anaconda.com/download/), and [Miniconda](https://conda.io/miniconda.html) (which is easier to install and maintain).

To exclude the possibility of dependencies conflicts its encouraged to create spearate conda environment for Plasflow using command:

```
conda create --name plasflow python=3.5
```

Python 3.5 is required becuase of TensorFlow requirements.


to activate created environment type:

```
source activate plasflow
```

PlasFlow can be easily installed as an Anaconda package from my Anaconda channel using:

```
conda install plasflow -c smaegol
```

When you decide to finish your work with PlasFlow, you can simply deactivate current anaconda environment with command:

```
source deactivate
```

##### Pip installer

There is a possibility of pip based installation. However, some requirements have to be met:
1. Python 3.5 is required (due to TensorFlow requirements)
2. TensorFlow has to be installed manually:

```
pip install https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.10.0rc0-cp35-cp35m-linux_x86_64.whl
```

then install PlasFlow with

```
pip install plasflow
```

However, models used for prediction have to be downloaded separately.


##### Manual installation

Of course, PlasFlow repo can be cloned using

```
git clone https://github.com/smaegol/PlasFlow
```

but in that case all dependencies have to be installed manually.
TensorFlow can be installed as specified above:

```
pip install https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.10.0rc0-cp35-cp35m-linux_x86_64.whl
```

python dependencies can be installed using pip:

```
pip install numpy pandas scipy rpy2 scikit-learn
```

to install R Biostrings go to https://bioconductor.org/packages/release/bioc/html/Biostrings.html and follow instructions therein


## Getting started

PlasFlow is designed to take a metagenomic assembly and identify contigs which may come from plasmids. It outputs several files, from which the most important is a tabular file containing all predictions (specified with `--output` option).

Options available in PlasFlow include:

- `--input` - specifies input fasta file with assembly contigs to classify [required]
- `--output` - a name of the tsv file with the tabular output of classification [required]
- `--threshold` - manually specified threshold for probability filtering (default = 0.7)

To invoke PlasFlow on `test.fasta` dataset (available in the test folder) simply:

```
PlasFlow.py --input test.fasta --output test.fasta.plasflow_predictions.tsv --threshold 0.7
```

## Output

The most important output of PlasFlow is a tabular file containing all predictions (specified with `--output` option), consiting of several columns including:

contig_id | contig_name | contig_length | id | label | ...
--------- | ----------- | ------------- | -- | ----- | ---
|

where:

- `contig_id`is an internal id of sequence used for the classification
- `contig_name` is a name of contig used in the classification
- `contig_length` shows the length of a classified sequence
- `id` is an internal id of a produced label (classification)
- `label` is the actual classification
- `...` represents additional columns showing probabilities of assignment to each possible class

Sequences can be classified to 26 classes including: chromosome.Acidobacteria, chromosome.Actinobacteria, chromosome.Bacteroidetes, chromosome.Chlamydiae, chromosome.Chlorobi, chromosome.Chloroflexi, chromosome.Cyanobacteria, chromosome.DeinococcusThermus, chromosome.Firmicutes, chromosome.Fusobacteria, chromosome.Nitrospirae, chromosome.other, chromosome.Planctomycetes, chromosome.Proteobacteria, chromosome.Spirochaetes, chromosome.Tenericutes, chromosome.Thermotogae, chromosome.Verrucomicrobia, plasmid.Actinobacteria, plasmid.Bacteroidetes, plasmid.Chlamydiae, plasmid.Cyanobacteria, plasmid.DeinococcusThermus, plasmid.Firmicutes, plasmid.Fusobacteria, plasmid.other, plasmid.Proteobacteria, plasmid.Spirochaetes.

If the probability of assignment to given class is lower than threshold (default = 0.7) then the sequence is treated as unclassified.

Additionaly, PlasFlow produces fasta files containing input sequences binned to plasmids, chromosomes and unclassified.

## Test dataset

Test dataset is located in the `test` folder (file `Citrobacter_freundii_strain_CAV1321_scaffolds.fasta`). It is the SPAdes 3.9.1 assembly of Citrobacter freundii strain CAV1321 genome (NCBI assembly ID: GCA_001022155.1), which contains 1 chromosome and 9 plasmids. In the `test_output` subfolder results of classification can be found

## Detailed information

Detailed information concerning the alogrithm and assumptions on which the PlasFlow is based can be found in the publication "*PlasFlow - Predicting Plasmid Sequences in Metagenomic Data Using Genome Signatures*" (*Nucleic Acids Research*, submitted). The flowchart illustrating major steps of training and prediction is shown below

![PlasFlow Flowchart](https://github.com/smaegol/PlasFlow/blob/master/flowchart.png)
