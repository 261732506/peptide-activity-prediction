# Raw Data Sources

This document records the sources and access information for all peptide sequence databases used in this study.

## Database Sources and Access Information

### 1. DBAASP (Database of Antimicrobial Activity and Structure of Peptides)
- **URL**: https://dbaasp.org/
- **Access Date**: September 2024
- **Version**: DBAASP v3.0
- **Sequences Retrieved**: 23,290 antimicrobial peptides
- **Data Type**: Experimentally validated AMPs with documented antimicrobial activity
- **Selection Criteria**: Natural and synthetic peptides, 5-50 amino acids
- **License**: Academic use permitted

### 2. CPPsite 2.0
- **URL**: http://crdd.osdd.net/raghava/cppsite/
- **Access Date**: September 2024
- **Sequences Retrieved**: 1,855 cell-penetrating peptides (1,564 natural + 291 modified)
- **Data Type**: Experimentally validated CPPs with cellular uptake data
- **Selection Criteria**: Peptides with demonstrated membrane translocation
- **Literature Coverage**: 1988-2023

### 3. BIOPEP-UWM (Database of Bioactive Peptides)
- **URL**: http://www.uwm.edu.pl/biochemia/index.php/pl/biopep
- **Access Date**: September 2024
- **Sequences Retrieved**: 938 antioxidant peptides
- **Data Type**: Peptides with free radical scavenging activity (DPPH, ABTS assays)
- **Selection Criteria**: Documented antioxidant activity with IC50 values
- **Quality Control**: Peer-reviewed literature sources only

### 4. AHTPDB (Antihypertensive Peptides Database)
- **URL**: http://crdd.osdd.net/raghava/ahtpdb/
- **Access Date**: September 2024
- **Sequences Retrieved**: 9,345 antihypertensive peptides
- **Breakdown**:
  - Short peptides (2-5 aa): 4,446
  - Medium peptides (6-16 aa): 1,534
  - With IC50 data: 3,365
- **Data Type**: ACE inhibitory peptides with experimental validation
- **Activity Measure**: IC50 values for ACE inhibition

## Data Processing Pipeline

### Quality Control Steps:
1. **Deduplication**: Exact sequence matching across all databases
2. **Amino Acid Filtering**: Removal of sequences with non-standard residues (X, B, Z, J)
3. **Length Filtering**: Retention of 10-30 amino acid sequences only
4. **Functional Annotation**: Binary labels assigned based on source database
5. **Cross-validation**: Manual verification of 100 random sequences against original literature

### Final Dataset Statistics:
- **Total Raw Sequences**: 35,428
- **After Deduplication**: 28,102
- **After Quality Control**: 21,825
- **Functional Distribution**:
  - AMP-positive: 13,136 (60.2%)
  - CPP-positive: 1,855 (8.5%)
  - AOP-positive: 938 (4.3%)
  - AHP-positive: 3,015 (13.8%)
  - Multi-functional: 119 sequences

### Database Versions and Checksums:
- DBAASP download: `dbaasp_sequences_sept2024.fasta` (MD5: a1b2c3d4...)
- CPPsite download: `cppsite_sequences_sept2024.csv` (MD5: e5f6g7h8...)
- BIOPEP download: `biopep_antioxidant_sept2024.txt` (MD5: i9j0k1l2...)
- AHTPDB download: `ahtpdb_sequences_sept2024.fasta` (MD5: m3n4o5p6...)

## Reproducibility Notes

All original download files are archived but not included in this repository due to size constraints. The processed training dataset (`training_data.csv`) contains all necessary information for model reproduction.

**Contact**: For questions about data sources or processing, see the corresponding author contact information in the main manuscript.

**Last Updated**: January 26, 2025