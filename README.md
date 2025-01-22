# MQuad Complementing Scripts

This tutorial provides guidance on using two Python scripts that complement the MQuad model for selecting informative mutations and analyzing mitochondrial DNA (mtDNA) mutations. These scripts provide 1) the beta-binomial mode for MQuad; 2) the ability to select a predefined number of top informative variants based on deltaBIC scores.
The original MQuad model is available at: https://github.com/single-cell-genetics/MQuad/tree/c2d750c5f279dca4f200274fc126d89e619fdd58

This repository will **keep updating** to fulfill requirements from users. If you're interested in those updates, please **star** this repository to stay informed!

### Setting Up the Environment

To use these scripts, ensure the required Python environment is properly configured. Follow these steps:

1. **Install Dependencies**:
   Make sure you have Python 3.7 or later installed, and install the required packages using `pip`:
   ```bash
   pip install numpy pandas scipy matplotlib seaborn vireoSNP mquad
   ```

2. **Clone the Repository**:
   Clone the repository containing the scripts:
   ```bash
   git clone https://github.com/linxy29/MQuad2.git
   cd mquad-scripts
   ```

3. **Prepare Input Data**:
   Ensure you have the necessary input files, including AD and DP matrices and a VCF file containing variant information. These can be generated using tools like cellSNP.

4. **Test the Setup**:
   Run a quick test using the example dataset provided with MQuad. The dataset includes:
   - 500 background variants.
   - 9 variants highlighted in Ludwig et al., *Cell*, 2019 (Supp Fig. 2F and main Fig. 2F).
   - 1 additional informative variant not mentioned in the paper.

   **Test Command**:
   ```bash
   mquad --vcfData example/example.vcf.gz -o example_test -p 5
   ```
   or with batch mode:
   ```bash
   mquad --vcfData example/example.vcf.gz -o example_test -p 5 --batchFit 1 --batchSize 5
   ```

Once the environment is set up, you can use the provided scripts for extended analysis.

## Script 1: Selecting Informative Mutations (`test_betabin.py`)

### Purpose

This script identifies informative mutations using the MQuad model. The key parameters are dynamically configurable via command-line arguments.

### Usage

```bash
python test_betabin.py --input_folder <input_folder> --output_folder <output_folder> --nproc <nproc> --minDP <minDP>
```

### Arguments

- `--input_folder`: Path to the input folder containing `cellSNP` output files.
- `--output_folder`: Path to the folder where results will be saved.
- `--nproc`: Number of processors to use for parallel computation (default: 16).
- `--minDP`: Minimum depth of coverage for variants (default: 2).

### Example

```bash
python test_betabin.py --input_folder ./example_input --output_folder ./example_output --nproc 32 --minDP 5
```

### Output

The script saves the following:

- Delta BIC values.
- A list of selected informative variants.
- Optionally, heatmaps and allele frequency matrices.

## Script 2: Selecting Predefined Number of Variants (`run_numbcutoff.py`)

### Purpose

This script allows selecting a predefined number of top informative variants based on deltaBIC scores.

### Usage

```bash
python run_numbcutoff.py --input_folder <input_folder> --output_folder <output_folder> --num_variants <num_variants> --nproc <nproc> --minDP <minDP>
```

### Arguments

- `--input_folder`: Path to the input folder containing `cellSNP` output files.
- `--output_folder`: Path to the folder where results will be saved.
- `--num_variants`: Number of top variants to select (default: None).
- `--nproc`: Number of processors to use for parallel computation (default: 32).
- `--minDP`: Minimum depth of coverage for variants (default: 2).

### Example

```bash
python run_numbcutoff.py --input_folder ./example_input --output_folder ./example_output --num_variants 2000 --nproc 32 --minDP 5
```

### Output

The script saves:

- The sorted deltaBIC values.
- The top N informative variants as specified by the `--num_variants` argument.
- Heatmaps and matrices of allele frequencies (optional).

## General Notes

- **Heatmaps**: Heatmaps are generated to visualize allele frequencies of selected variants. These are saved in the output folder.
- **Output Matrices**: Passed AD and DP matrices are exported as `.mtx` files for further analysis.
- **Knee Point Detection**: If no cutoff or predefined number of variants is provided, the scripts automatically determine the optimal threshold for variant selection using the knee point method.

## Common Questions

### 1. What is the difference between the two scripts?
- `test_betabin.py` is designed for identifying informative mutations using the beta-binomial model in MQuad.
- `run_numbcutoff.py` allows you to specify a fixed number of top variants based on deltaBIC scores, offering flexibility for predefined variant selection.

### 2. Can I use these scripts independently of MQuad?
No, these scripts complement the MQuad model and rely on its outputs (e.g., deltaBIC scores and variant data).

### 3. How do I prepare input data for these scripts?
Ensure you have:
- AD and DP matrices (e.g., from cellSNP outputs).
- A VCF file containing variant information.

### 4. How do I troubleshoot common errors?
- **File Not Found**: Verify the input file paths and folder permissions.
- **Permission Denied**: Ensure you have write permissions for the output directory.
- **Incomplete Outputs**: Check for memory or CPU limitations and consider increasing resources.

### 5. Can I modify parameters for better performance?
Yes, both scripts allow customization of parameters like `--nproc`, `--minDP`, and `--num_variants`. Adjust these based on your dataset and computational resources.

### 6. What is the recommended cutoff for selecting variants?
If unsure, allow the script to automatically determine the cutoff using the knee point method. Alternatively, use domain knowledge to set a custom `--num_variants` or `--minDP`.



