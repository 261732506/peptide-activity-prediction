# Cross-Validation Data Splits

This directory contains the exact data splits used for cross-validation experiments to ensure reproducibility.

## Files:
- `cv_folds.json` - 5-fold cross-validation indices
- `train_indices.txt` - Training set indices
- `test_indices.txt` - Test set indices
- `validation_indices.txt` - Validation set indices

Generated from the full training dataset (21,825 sequences) using stratified sampling to maintain class balance across folds.

## Usage:
```python
import json
import pandas as pd

# Load CV folds
with open('cv_folds.json', 'r') as f:
    cv_folds = json.load(f)

# Load training data
df = pd.read_csv('../processed/training_data.csv')

# Access fold 0 training data
fold_0_train = df.iloc[cv_folds['fold_0']['train']]
fold_0_val = df.iloc[cv_folds['fold_0']['val']]
```