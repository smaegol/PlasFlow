#!/usr/bin/env python

#######################################################################################
###                                                                                 ###
###     PlasFlow 1.1                                                                ###
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
import numpy as np
import pandas as pd
import re
from Bio import SeqIO

from os.path import join

from .classifier import Tf_Classif
from .vote_classifier import TF_Vote_Classifier
from .constants import MODELS_PATH
from .utils import re_select_columns


def predict(inputfile):
    vote_class = TF_Vote_Classifier(clfs=[
        Tf_Classif(5, "20_20"),
        Tf_Classif(6, "20_20"),
        Tf_Classif(7, "30"),
    ])
    vote_proba = vote_class.predict_proba(inputfile)
    vote = vote_class.predict(inputfile)
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

    return pd_n, pd_n_proba


def parse_seqs(inputfile):
    """Return a list of the seqs and a pandas dataframe with metadata."""
    print("Importing sequences")
    seqs = list(SeqIO.parse(inputfile, 'fasta'))
    print("Imported ", len(seqs), " sequences")
    contig_info = pd.DataFrame([
        {'contig_id': i, 'contig_name': rec.id, 'contig_length': len(rec.seq)}
        for i, rec in enumerate(seqs)
    ])
    return seqs, contig_info


def write_seqs(inputfile, outputfile, results_merged):
    sequences_dict = SeqIO.index(inputfile, "fasta")
    plasmid_sequences, chromosome_sequences, unclassified_sequences = [], [], []
    for index, row in results_merged.iterrows():
        processed_sequence = sequences_dict[row.contig_name]
        processed_sequence.id = processed_sequence.id + " " + row.label
        if re.match(r'^chromosome.*', row.label):
            chromosome_sequences.append(processed_sequence)
        elif re.match(r'^plasmid.*', row.label):
            plasmid_sequences.append(processed_sequence)
        else:
            unclassified_sequences.append(processed_sequence)

    SeqIO.write(chromosome_sequences, outputfile + "_chromosomes.fasta", "fasta")
    SeqIO.write(plasmid_sequences, outputfile + "_plasmids.fasta", "fasta")
    SeqIO.write(unclassified_sequences, outputfile + "_unclassified.fasta", "fasta")


def filter_results(results_merged, threshold):
    for index, row in results_merged.iterrows():
        taxname = row.label.split(".", 1)[1]
        if row[row.label] > threshold:
            continue

        plasmids = re_select_columns(row, results_merged, r'^plasmid.*')
        chromosomes = re_select_columns(row, results_merged, r'^chromosom.*')
        taxnames = re_select_columns(row, results_merged, r".*" + re.escape(taxname) + r"")
        if plasmids.sum() > threshold:
            results_merged.at[index, 'label'] = 'plasmid.unclassified'
        elif chromosomes.sum() > threshold:
            results_merged.at[index, 'label'] = 'chromosome.unclassified'
        elif taxnames.sum() > threshold:
            results_merged.at[index, 'label'] = 'unclassified.' + taxname
        else:
            results_merged.at[index, 'label'] = 'unclassified.unclassified'
    return results_merged


def main(labels, inputfile, outputfile, threshold):
    seqs, contig_info = parse_seqs(inputfile)
    vote, vote_proba = predict(inputfile)

    results_merged = pd.merge(vote, labels, on=['id'])
    results_merged = pd.merge(results_merged, vote_proba, on=['contig_id'])
    results_merged = pd.merge(contig_info, results_merged, on=['contig_id'])

    print("Filtering by probability threshold", threshold)
    results_merged = filter_results(results_merged, threshold)
    results_merged.to_csv(outputfile, sep='\t')

    taxons, plasmids = {}, {}
    for index, row in results_merged.iterrows():
        plasmid, taxname = row.label.split(".", 1)[0], row.label.split(".", 1)[1]
        taxons[taxname] = taxons.get(taxname, 0) + 1
        plasmids[plasmid] = plasmids.get(taxname, 0) + 1
    plasmids = pd.DataFrame.from_dict(plasmids, orient="index").transpose()
    print(f"\nResulting plasmid sequences prediction:\n{plasmids}")
    taxons = pd.DataFrame.from_dict(taxons, orient="index").transpose()
    print(f"\nResulting taxonomical assignment:\n{taxons}")

    print("\nOutputting fasta files with classified sequences")
    write_seqs(inputfile, outputfile, results_merged)
