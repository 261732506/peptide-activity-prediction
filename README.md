# Dual-Functional Peptide Prediction Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://doi.org/)

An integrated computational pipeline for predicting and prioritizing dual-functional antimicrobial and cell-penetrating peptide candidates using protein language models and machine learning.

## 📄 Publication

**Title**: An integrated computational pipeline for prioritizing dual-functional antimicrobial and cell-penetrating peptide candidates

**Journal**: PLOS ONE (submitted)

**Authors**: [Your Name], [Institution]

**Citation**: [To be added upon publication]

---

## 🎯 Overview

This repository provides a complete, reproducible implementation of our dual-functional peptide prediction and design pipeline combining:

- **ESM-2 protein language model** embeddings (2560 dimensions)
- **Physicochemical features** (15 dimensions): charge, hydrophobicity, Boman index, etc.
- **Sequence motif features** (15 dimensions): function-specific patterns

The pipeline achieves **84-88% accuracy** for simultaneous AMP (antimicrobial peptide) and CPP (cell-penetrating peptide) prediction, comparable to established single-function tools.

---

## ✨ Key Features

- ✅ **Multi-functional prediction**: Simultaneous prediction of AMP, CPP, AOP, AHP
- ✅ **State-of-the-art embeddings**: ESM-2-3B protein language model
- ✅ **Two design strategies**: Point mutation editing and modular assembly
- ✅ **Rigorous validation**: Cross-validation, literature validation, comparison with 6 existing tools
- ✅ **Fully reproducible**: All data, models, and code publicly available
- ✅ **Well-documented**: Comprehensive documentation and usage examples

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/261732506/peptide-activity-prediction.git
cd peptide-activity-prediction

# Create conda environment
conda create -n peptide-pred python=3.9
conda activate peptide-pred

# Install dependencies
pip install -r requirements.txt

# Install fair-esm for ESM-2
pip install fair-esm
```

### Basic Usage

```python
import pickle
import numpy as np

# Load trained models
with open('models/improved_predictors_optimized.pkl', 'rb') as f:
    predictors = pickle.load(f)

# The models dictionary contains Random Forest classifiers for each function
# predictors = {'AMP': RandomForestClassifier, 'CPP': RandomForestClassifier, ...}

# To use the models, you need to extract features from sequences
# See src/peptide_prediction.py for the ESMPeptideAnalyzer class
# which extracts ESM-2 embeddings + physicochemical + motif features (2590-dim)

sequence = "RKKRRQRRR"  # TAT peptide
print(f"Sequence: {sequence}")
# Full prediction requires feature extraction - see src/peptide_prediction.py
```

**Note**: The `ESMPeptideAnalyzer` class in `src/peptide_prediction.py` provides the complete feature extraction pipeline. The basic workflow is:
1. Extract ESM-2 embeddings (2560-dim)
2. Calculate physicochemical properties (15-dim)
3. Extract motif patterns (15-dim)
4. Concatenate to 2590-dimensional feature vector
5. Predict using trained Random Forest models

### Design Dual-Functional Candidates

```python
# Example: Modular assembly strategy (as described in paper)
from itertools import product

# Define modules (examples from Table 5)
cpp_modules = ['RRRRRRRRR', 'RKKRRQRRR', 'RRWRRWRR']
amp_modules = ['KLAKLAK', 'KRWWKWIRW', 'GIGKFLHSAKKF']
linkers = ['', 'G', 'GG', 'GGS', 'GGGGS']

# Generate candidates
candidates = []
for cpp, amp, linker in product(cpp_modules, amp_modules, linkers):
    candidate_seq = cpp + linker + amp
    if 15 <= len(candidate_seq) <= 25:
        candidates.append(candidate_seq)
        # Score with predictors (requires feature extraction)
        # features = extract_features(candidate_seq)
        # cpp_prob = predictors['CPP'].predict_proba([features])[0][1]
        # amp_prob = predictors['AMP'].predict_proba([features])[0][1]
        # joint_prob = cpp_prob * amp_prob

print(f"Generated {len(candidates)} candidate sequences")
# Example output: ['RRRRRRRRRKLAKL AK', 'RRRRRRRRRGKLAKLAK', ...]
```

**See `scripts/modular_assembly.py` and `scripts/point_mutation.py` for complete implementations.**

---

## 📊 Dataset

### Training Data

21,825 unique peptide sequences (10-30 amino acids) from four public databases:

| Database | Function | Sequences | URL |
|----------|----------|-----------|-----|
| DBAASP | AMP | 18,218 | https://dbaasp.org/ |
| CPPsite 2.0 | CPP | 1,328 | http://crdd.osdd.net/raghava/cppsite/ |
| BIOPEP-UWM | AOP | 935 | http://www.uwm.edu.pl/biochemia/biopep |
| AHTPDB | AHP | 1,693 | http://crdd.osdd.net/raghava/ahtpdb/ |

**File**: `data/processed/training_data.csv`

### Validation Data

Independent literature-curated set:
- **CPP positives**: TAT, R9, Penetratin, Pep-1, etc. (n=15)
- **AMP positives**: LL-37, Magainin-2, Melittin, etc. (n=15)
- **Negatives**: Poly-amino acids, epitope tags (n=10 per function)

**File**: `data/validation/literature_validation.csv`

---

## 🧬 Features

### ESM-2 Embeddings (2560d)

```python
from esm_features import extract_esm2_features

# Extract ESM-2 embeddings
embeddings = extract_esm2_features(
    sequences=['RKKRRQRRR'],
    model='esm2_t36_3B_UR50D',
    layer=36,
    pooling='mean'
)
```

### Physicochemical Features (15d)

1. Sequence length
2-7. Amino acid composition (hydrophobic, positive, negative, aromatic, polar, small)
8-10. Special residues (proline, cysteine, glycine)
11. Net charge (pH 7.4)
12. Hydrophobic moment
13. Boman index
14. Instability index
15. Estimated pI

### Sequence Motif Features (15d)

- AMP motifs: KLAK, RWR, RRWW, etc.
- CPP motifs: RRR, RKKR, etc.
- AOP motifs: HH, YY, HY
- AHP motifs: IPP, VPP, LPP
- Repeat counts: KK, RR, WW, PP

---

## 📈 Performance

### Cross-Validation (5-fold)

| Classifier | Accuracy | Precision | Recall | F1-Score |
|-----------|----------|-----------|--------|----------|
| AMP (τ=0.80) | 87.6 ± 1.0% | 85.6% | 89.6% | 87.6% |
| CPP (τ=0.50) | 84.4 ± 1.0% | 83.4% | 85.6% | 84.4% |
| AOP (τ=0.50) | 82.3 ± 1.5% | 81.2% | 84.5% | 82.8% |
| AHP (τ=0.50) | 83.1 ± 1.2% | 81.9% | 85.2% | 83.5% |

### Comparison with Existing Tools

| Task | Our Model | AMPlify | iAMP-2L | DBAASP | CellPPD | MLCPP |
|------|-----------|---------|---------|--------|---------|-------|
| AMP | 88.0% | 89.0% | 85.0% | 87.0% | - | - |
| CPP | 84.0% | - | - | - | 82.0% | 83.0% |

### Design Strategy Comparison

| Strategy | Candidates | Max Joint Prob | Improvement |
|----------|-----------|----------------|-------------|
| Point Mutation | 304 | 0.308 | Baseline |
| Modular Assembly | 193 | 0.364 | **+18.3%** |

---

## 🔬 Key Results

### Top Dual-Functional Candidate

**Sequence**: `RRRRRRRRRGGGGSKRWWKWIRW`

**Structure**: R9 (CPP module) + GGGGS (linker) + KRWWKWIRW (AMP motif)

**Predictions**:
- CPP probability: 0.639
- AMP probability: 0.570
- Joint probability: 0.364

**Properties**:
- Length: 23 amino acids
- Net charge: +10 (pH 7.4)
- Hydrophobic content: 17.4%
- Aromatic content: 13.0%

---

## 📁 Repository Structure

```
peptide-activity-prediction/
├── data/
│   ├── raw/                    # Raw sequences from databases
│   ├── processed/              # Preprocessed training data
│   │   └── training_data.csv   # 21,825 sequences
│   └── validation/             # Validation sets
├── models/
│   └── improved_predictors_optimized.pkl  # Trained Random Forest models
├── features/
│   └── extracted_features.npy  # Pre-computed 2590d feature matrix
├── scripts/
│   ├── peptide_editing_system.py  # Main prediction system
│   ├── extract_esm2_features.py   # ESM-2 feature extraction
│   ├── train_models.py            # Model training
│   ├── cross_validation.py        # Cross-validation
│   ├── literature_validation.py   # Literature validation
│   ├── comparison_with_tools.py   # Benchmark comparison
│   └── ablation_study.py          # Feature ablation
├── design/
│   ├── modular_assembly.py     # Modular design strategy
│   └── point_mutation.py       # Point mutation editing
├── results/
│   ├── figures/                # All publication figures (300 DPI)
│   └── tables/                 # All publication tables
├── docs/
│   ├── README.md              # This file
│   ├── USAGE.md               # Detailed usage guide
│   └── API.md                 # API documentation
├── requirements.txt           # Python dependencies
├── environment.yml            # Conda environment
├── LICENSE                    # MIT License
└── CITATION.cff              # Citation metadata
```

---

## 🛠️ Dependencies

### Core Requirements

- Python ≥ 3.9
- PyTorch ≥ 1.12
- fair-esm ≥ 2.0.0
- scikit-learn ≥ 1.3.0
- NumPy ≥ 1.23.5
- pandas ≥ 1.5.3

### Optional (for visualization)

- matplotlib ≥ 3.5.3
- seaborn ≥ 0.12.2

See `requirements.txt` for complete list.

---

## 📖 Documentation

### Detailed Guides

- **[Usage Guide](docs/USAGE.md)**: Step-by-step tutorials
- **[API Reference](docs/API.md)**: Complete API documentation
- **[Methods](docs/METHODS.md)**: Detailed methodology
- **[FAQ](docs/FAQ.md)**: Frequently asked questions

### Jupyter Notebooks

- `notebooks/01_Feature_Extraction.ipynb`: Extract features from sequences
- `notebooks/02_Model_Training.ipynb`: Train custom models
- `notebooks/03_Prediction.ipynb`: Predict peptide functions
- `notebooks/04_Design.ipynb`: Design dual-functional candidates
- `notebooks/05_Visualization.ipynb`: Visualize results

---

## 🔍 Validation

### Cross-Validation

```bash
python scripts/cross_validation.py --n-folds 5
```

### Literature Validation

```bash
python scripts/literature_validation.py \
    --validation-set data/validation/literature_validation.csv
```

### Comparison with Existing Tools

```bash
python scripts/comparison_with_tools.py
```

### Ablation Study

```bash
python scripts/ablation_study.py
```

---

## 💡 Example Use Cases

### 1. Predict Function of Novel Peptide

```python
# Load predictor
predictor = DualFunctionalPredictor('models/improved_predictors_optimized.pkl')

# Your novel sequence
novel_seq = "KWKLFKKIEKVGQN"

# Predict all functions
results = predictor.predict_all_functions(novel_seq)

for func, prob in results.items():
    print(f"{func}: {prob:.3f}")
```

### 2. Batch Prediction

```python
import pandas as df

# Load sequences
sequences = pd.read_csv('my_sequences.csv')['sequence'].tolist()

# Batch predict
batch_results = predictor.predict_batch(sequences, batch_size=32)

# Save results
batch_results.to_csv('predictions.csv', index=False)
```

### 3. Design Custom Dual-Functional Peptide

```python
# Define custom modules
my_cpp_modules = ['RKKRRQRRR', 'RRRRRRRRR']
my_amp_modules = ['KLAKLAK', 'KRWWKWIRW']
my_linkers = ['GGS', 'GGGGS']

# Generate candidates
designer = ModularAssembly(predictor)
candidates = designer.generate_candidates(
    cpp_modules=my_cpp_modules,
    amp_modules=my_amp_modules,
    linkers=my_linkers,
    length_min=15,
    length_max=25
)

# Rank and export
top = designer.rank_by_joint_probability(candidates, n=20)
top.to_csv('my_dual_candidates.csv', index=False)
```

---

## ⚠️ Limitations

Please note the following important limitations:

1. **No experimental validation**: All results are computational predictions requiring wet-lab verification
2. **Binary classification**: Functions treated as binary (present/absent), not quantitative (e.g., MIC values)
3. **Sequence-based only**: No explicit 3D structure or dynamics modeling
4. **Training data biases**: Limited by coverage and annotations in public databases
5. **Context-independent**: Does not model salt sensitivity, serum stability, or target specificity

**Recommendations**:
- Use predictions to *prioritize* candidates for experimental validation
- Validate top candidates through MIC assays, cellular uptake experiments, and cytotoxicity assays
- Consider structure prediction (AlphaFold2) and molecular dynamics for mechanistic insights

---

## 📊 Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{YourName2025,
  title={An integrated computational pipeline for prioritizing dual-functional antimicrobial and cell-penetrating peptide candidates},
  author={[Your Name]},
  journal={PLOS ONE},
  year={2025},
  doi={[To be added]},
  url={https://github.com/261732506/peptide-activity-prediction}
}
```

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution

- [ ] Add more design strategies (genetic algorithms, reinforcement learning)
- [ ] Implement quantitative prediction (MIC, EC50)
- [ ] Add structure-aware features (AlphaFold2 integration)
- [ ] Expand to additional functions (anticancer, antiviral, etc.)
- [ ] Web interface for easy access
- [ ] Docker container for reproducibility

---

## 🆘 Support

- **Issues**: Please open an [issue](https://github.com/261732506/peptide-activity-prediction/issues)
- **Email**: [your.email@institution.edu]
- **Documentation**: [https://dual-functional-peptides.readthedocs.io](https://dual-functional-peptides.readthedocs.io)

---

## 🙏 Acknowledgments

- ESM-2 model: [fair-esm](https://github.com/facebookresearch/esm) by Meta AI Research
- Database maintainers: DBAASP, CPPsite 2.0, BIOPEP-UWM, AHTPDB
- Open-source community: scikit-learn, PyTorch, NumPy, pandas

---

## 📅 Version History

### v1.0.0 (2025-02-01)
- Initial release with PLOS ONE submission
- Complete pipeline for dual-functional peptide prediction
- Two design strategies implemented
- Comprehensive validation framework

---

## 🔗 Related Resources

- **Databases**:
  - [DBAASP](https://dbaasp.org/) - Antimicrobial peptides
  - [CPPsite 2.0](http://crdd.osdd.net/raghava/cppsite/) - Cell-penetrating peptides
  - [APD3](http://aps.unmc.edu/AP/) - Antimicrobial Peptide Database

- **Prediction Tools**:
  - [AMPlify](https://github.com/bcgsc/AMPlify) - AMP prediction
  - [CellPPD](http://crdd.osdd.net/raghava/cellppd/) - CPP prediction
  - [MLCPP](http://www.thegleelab.org/MLCPP/) - Machine learning CPP prediction

- **Protein Language Models**:
  - [ESM-2](https://github.com/facebookresearch/esm) - Evolutionary Scale Modeling
  - [ProtTrans](https://github.com/agemagician/ProtTrans) - Protein transformers

---

**Last updated**: January 2025

**Maintained by**: [Your Name] ([your.email@institution.edu](mailto:your.email@institution.edu))

**Project Status**: Active Development

**Keywords**: antimicrobial peptides, cell-penetrating peptides, protein language models, ESM-2, machine learning, peptide design, dual-functional peptides, computational biology
