import sys
import os
import pandas as pd

def main(filepath, fileout=None):
    if not fileout:
        fname, ext = os.path.splitext(filepath)
        fileout = fname+'_txid' + '.' + ext
    df = pd.read_table(filepath)
    df['target_id'] = df['target_id'].str.split('|').str.get(0)
    df.to_csv(fileout, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[2])
    else:
        print('Run: python fix_kallisto_output.py <abundance.tsv> <fileout.tsv>')

