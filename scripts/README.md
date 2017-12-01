Scripts used during the development of [PlasFlow](https://github.com/smaegol/PlasFlow) and manuscript preparation.

# Table of contents

- [Sequence preprocessing](#sequence-preprocessing)
- [Neural network training](#neural-network-training)
- [Other prediction software testing](#other-prediction-software-testing)

  - [cBar](#cbar)
  - [Recycler](#recycler)

## Sequence preprocessing

All preprocessing steps producing kmer counts used for the neural network training can be run using command:

```
  bash PlasFlow_preprocessing.sh
```

Alternatively, all steps can be ran individually:

- Use `get_refseq_genomes_complete.sh` to download complete genomes (Archaea and Bacteria)
- Run `get_refseq_catalog_filter.sh` to get accessions numbers of each downloaded sequence and filtered RefSeq catalog files
- Run `Rscript create_annotation_of_data.R` to create the annotation table based on the RefSeq catalog (filtered to contain only downloaded sequences).

  - Uses rentrez to collect data from NCBI - may be slow

- Run `RScript generate_fragments_for_training.R`

  - May require large amounts of RAM - especially for longer kmers (hexamers, heptamers)

## Neural network training

- Run training using `PlasFlow_train.py` script - requires TensorFlow 0.10.0 (see [requirements for PlasFlow](https://github.com/smaegol/PlasFlow#requirements))

  - As an input `kmern_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv` files should be used (where **n** is a kmer length)
  - Typical invocation of a script:

  ```
  python PlasFlow_train.py --input kmern_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv --hidden1 20 --hidden2 10 --activation relu --modeldir kmern_split_20_10_neurons_relu --steps 50000
  ```

  will run training on specified file, using 20 neurons in first hidden layer and 10 neurons in second hidden layer, using `relu` activation function and 50000 steps of training. Model will be saved in the `kmern_split_20_10_neurons_relu` folder

## Other prediction software testing

Scripts used for plasmid prediction using other available software are also provided. The results of comparison to PlasFlow can be found in the publication manuscript (coming soon).

### cBar

cBar([1](https://www.ncbi.nlm.nih.gov/pubmed/20538725)) can be downloaded from its [website](http://csbl.bmb.uga.edu/~ffzhou/cBar)

To invoke cBar prediction on test dataset use `calculate_cbar_predictions.sh`. It is required to provide cBar path at the top of the script and that `get_seqs_by_ids_list_cBar.pl` and `get_seqs_by_ids_list_cBar.pl` scripts are available and their paths are also provided.

Run the script using:

```
bash calculate_cbar_predictions.sh -i test.fasta
```

This will produce 3 output files:

- `test.fasta_cbar.tsv` - containing raw cBar predictions
- `test.fasta_cbar.tsv_plasmids.fasta` - fasta file containing sequences predicted to come from plasmids
- `test.fasta_cbar.tsv_chromosomes.fasta` - fasta file containing sequence predicted to come from chromosomes

### Recycler

Recycler([2](https://www.ncbi.nlm.nih.gov/pubmed/28003256)) can be downloaded from [github](https://github.com/Shamir-Lab/Recycler).

To invoke Recycler prediction on any dataset it is required to have assembly graph (fastg) file from SPAdes([3](https://www.ncbi.nlm.nih.gov/pubmed/22506599)) assembly. Here we provide the script which performs the analysis starting with the paired-end sequencing files (`calculate_recycler_predictions.sh`).

To run the scripts it is required to have installed:

- SPADes (can be obtained from its [website](http://bioinf.spbau.ru/spades))
- BWA([4](https://www.ncbi.nlm.nih.gov/pubmed/19451168)) (can be obtained form its [website](http://bio-bwa.sourceforge.net/))
- Samtools([5](https://www.ncbi.nlm.nih.gov/pubmed/19505943)) (can be obtained from the [website](http://samtools.sourceforge.net/))
- Recycler

All the paths have to be specified in the top section of the script.

Run the script using:

```
bash calculate_recycler_predictions.sh -1 input_R1.fastq -2 input_R2.fastq
```

SPAdes assembly will be run using provided input files and located in the output dir (`input_R1_spades`). Then the mapping of reads to contigs is done and Recycler is run. The output of Recycler (fasta with circular contigs) will be available in the `input_R1_spades/assembly_graph.cycs.fasta` file.

## REFERENCES

1. Zhou, F., and Xu, Y. (2010). cBar: a computer program to distinguish plasmid-derived from chromosome-derived sequence fragments in metagenomics data. Bioinformatics 26, 2051–2052.
2. Rozov, R., Brown Kav, A., Bogumil, D., Shterzer, N., Halperin, E., Mizrahi, I., and Shamir, R. (2017). Recycler: an algorithm for detecting plasmids from de novo assembly graphs. Bioinformatics 33, 475–482.
3. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D., et al. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19, 455–477.
4. Li, H., and Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754–1760.
5. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., and Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079.
