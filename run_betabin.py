# Xinyi Lin, 202303
# Goal: Using the beta-binomial model to select informative mutations (MQuad2)
# Usage: python /home/linxy29/code/MQuad2/test_betabin.py --input_folder "path_to_input" --output_folder "path_to_output" --nproc 16 --minDP 2 --beta_mode True --minCell 10 --cutoff None

import os
import logging
import argparse
from scipy.io import mmread
from scipy.sparse import csc_matrix
from mquad.mquad import Mquad

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

def validate_files(files):
    """Validate the existence of required files."""
    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File not found: {file}")

def main(input_folder, output_folder, nproc, minDP, beta_mode, minCell, cutoff):
    """Main function to process data and run MQuad."""
    setup_logging()

    # Define file paths
    ad_file = os.path.join(input_folder, "cellSNP.tag.AD.mtx")
    dp_file = os.path.join(input_folder, "cellSNP.tag.DP.mtx")

    # Validate files
    validate_files([ad_file, dp_file])

    logging.info("Reading allele frequency data...")
    cell_dat = {
        "AD": csc_matrix(mmread(ad_file)),
        "DP": csc_matrix(mmread(dp_file)),
        "variants": [f"SNP{x + 1}" for x in range(csc_matrix(mmread(ad_file)).shape[0])],
    }
    logging.info("Loaded AD and DP files.")

    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    logging.info("Initializing MQuad...")
    mdphd = Mquad(AD=cell_dat['AD'], DP=cell_dat['DP'], variant_names=cell_dat['variants'])

    logging.info("Calculating deltaBIC...")
    df = mdphd.fit_deltaBIC(out_dir=output_folder, nproc=nproc, minDP=minDP, beta_mode=beta_mode)

    logging.info("Selecting informative variants...")
    best_ad, best_dp = mdphd.selectInformativeVariants(
        min_cells=minCell, out_dir=output_folder, tenx_cutoff=cutoff
    )

    logging.info("Process completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mitochondrial informative mutation selection for MAESTER data.")
    parser.add_argument("--input_folder", required=True, help="Path to the input folder containing cellSNP output files.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder.")
    parser.add_argument("--nproc", type=int, default=16, help="Number of processors to use.")
    parser.add_argument("--minDP", type=int, default=2, help="Minimum depth of coverage.")
    parser.add_argument("--beta_mode", type=bool, default=True, help="Use beta mode (True/False).")
    parser.add_argument("--minCell", type=int, default=10, help="Minimum number of cells.")
    parser.add_argument("--cutoff", type=float, default=None, help="Cutoff for tenx selection.")
    args = parser.parse_args()

    try:
        main(
            args.input_folder,
            args.output_folder,
            args.nproc,
            args.minDP,
            args.beta_mode,
            args.minCell,
            args.cutoff
        )
    except Exception as e:
        logging.error(f"An error occurred: {e}")
