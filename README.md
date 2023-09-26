# PolymerGE

## Description

Based on the paper "Population-based De Novo Molecule Generation, Using Grammatical Evolution"  ([Read Paper](https://doi.org/10.1246/cl.180665)). This version  adapts and extends the original Zinc grammar to represent and generate polymers, enabling specialized molecule design using grammatical evolution for polymers. 

## Usage

### Basic Usage

1. **Clone the Repository**

    ```bash
    git clone https://github.com/yourusername/PolymerGE.git
    ```

2. **Navigate to the Directory**

    ```bash
    cd PolymerGE
    ```

3. **Run `test.py`**

    ```bash
    python test.py
    ```

    Performing a test about encoding/decoding of polymer smiles in PI1M database through grammatical evolution. 

## Dependencies

- RDKit
- NumPy

## More Information

### Testing and Accuracy

We applied extensive testing with a dataset of 100,000 PI1M polymers. The encoding and decoding accuracy stands at 98.3%.

### Limitations

- Some deviations when working with large molecules.
- Smiles with large rings of size 10 or above cannot be computed, similar to the limitations found in the original ChemGE.



