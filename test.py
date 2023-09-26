from rdkit import rdBase
from gram_utils import *
from smiles_utils import *
import numpy as np

def main():
    rdBase.DisableLog('rdApp.error')

    filename = 'PI1M_clean.smi'
    smiles = []

    with open(filename) as f:
        for i, line in enumerate(f, start=1):
            orig = canonicalize(line)
            smiles.append(orig)
            if len(smiles) > 1000:
                break

    smi = canonicalize(smiles[777])
    c_gene = CFGtoGene(encode(smi), max_len=300)
    decoded_smi = canonicalize(decode(GenetoCFG(c_gene)))

    print('polymer smiles:', smi)
    print('c_gene:', c_gene)
    print('decoded polymer smiles:', decoded_smi)
    print('matched?', smi == decoded_smi)

if __name__ == "__main__":
    main()
