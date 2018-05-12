#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import pandas
import seaborn as sns


def main(f1, f2, outprefix):
    tpm1 = pandas.read_table(f1, header=None, names=['id', 'tpm1'], sep=' ')
    print(tpm1.head())
    tpm2 = pandas.read_table(f2, header=None, names=['id', 'tpm'], sep=' ')
    tpm1['tpm2'] = tpm2['tpm']
    ax = sns.jointplot('tpm1', 'tpm2', data=tpm1, kind='reg')
    ax.savefig('{}.png'.format(outprefix))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
