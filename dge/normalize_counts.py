#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='This script normalizes a count matrix to either TPM or FPKM')
parser.add_argument('-i', '--dataset', type=str, help='input filename - count matrix')
parser.add_argument('-d', '--delimiter', type=str, default='\t', help='field delimiter character')
parser.add_argument('-t', '--norm_type', type=str, default='tpm', choices=['tpm', 'fpkm'], help='type of normalization')
parser.add_argument('-c', '--lengths_column', type=int, default=2, help=('Which column contains gene length information'))
parser.add_argument('-s', '--samples_first_column', type=int, default=3, help=('Which is the first sample column?'))

args=parser.parse_args()

args.lengths_column = args.lengths_column - 1
args.samples_first_column = args.samples_first_column - 1
filebase, extension = args.dataset.split('.')[0], args.dataset.split('.')[1]

df = pd.read_csv(args.dataset, sep=args.delimiter)
df_samples = df.iloc[:, args.samples_first_column:]
length_column = df.iloc[:, args.lengths_column]
length_column = length_column / 1000

if args.norm_type == 'fpkm':
    
    scaling_factor = df_samples.sum() / 1000000
    df_samples = df_samples.div(scaling_factor, axis=1)
    df_samples = df_samples.div(length_column, axis=0)
    
    
elif args.norm_type == 'tpm':
    
    df_samples = df_samples.div(length_column, axis=0)
    scaling_factor = df_samples.sum() / 1000000
    df_samples = df_samples.div(scaling_factor, axis=1)


df_out = pd.concat([df.iloc[:,:args.samples_first_column], df_samples], axis=1)

df_out.to_csv(filebase + '_normalized.' + extension, sep=args.delimiter, index=False)
