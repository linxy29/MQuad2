# Xinyi Lin, 202303
# Goal: Using the beta-binomial model to select informative mutations (MQuad2) and provide pre-defined number of variants
# Usage: python /home/linxy29/code/MQuad2/run_numbcutoff.py --input_folder ./maester_assemble_trimmed_aligned_mt_addtag_cellSNP0 --output_folder ./maester_assemble_trimmed_aligned_mt_addtag_cellSNP0_betamquad_2000variant --num_variants 2000 --nproc 32 --minDP 2

import os
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
from mquad.mquad_utils import findKnee
from mquad.mquad import Mquad
from vireoSNP.utils.vcf_utils import read_sparse_GeneINFO, load_VCF

def select_informative_variants(self, min_cells=2, export_heatmap=True, export_mtx=True, out_dir=None, existing_df=None, tenx_cutoff=None, numb_cutoff=None):
    """
    Select informative variants based on deltaBIC and/or other thresholds.
    
    Args:
        self: Mquad instance.
        min_cells: Minimum number of cells.
        export_heatmap: Whether to export the heatmap.
        export_mtx: Whether to export AD and DP matrices.
        out_dir: Directory to save outputs.
        existing_df: Existing pre-fitted model data.
        tenx_cutoff: Cutoff value for tenx mode.
        numb_cutoff: Number of top variants to select.

    Returns:
        best_ad, best_dp: Sparse matrices of selected variants.
    """
    if existing_df is not None:
        print(f'[MQuad] Fitted model detected, using {existing_df}...')
        self.df = pd.read_csv(existing_df)
        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

    if out_dir and not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except PermissionError:
            print("[MQuad] Cannot create output directory, check permissions.")

    if numb_cutoff is not None and numb_cutoff > self.df.shape[0]:
        print('[MQuad] Number of variants to select exceeds total variants, disabling numb mode.')
        numb_cutoff = None

    if numb_cutoff is not None and numb_cutoff > 50:
        print('[MQuad] Number of variants to select exceeds 50, disabling heatmap export.')
        export_heatmap = False

    if tenx_cutoff is not None and numb_cutoff is not None:
        print('[MQuad] Both tenx and numb modes specified, disabling numb mode.')
        numb_cutoff = None

    if tenx_cutoff is not None:
        print(f'[MQuad] Tenx mode used with cutoff = {tenx_cutoff}')
        self.final_df = self.sorted_df[(self.sorted_df.deltaBIC >= tenx_cutoff) & (self.sorted_df.num_cells_minor_cpt >= min_cells)]
    elif numb_cutoff is not None:
        print(f'[MQuad] Numb mode used with cutoff = {numb_cutoff}')
        self.final_df = self.sorted_df.head(numb_cutoff)
        self.final_df = self.final_df[self.final_df.num_cells_minor_cpt >= min_cells]
    else:
        print('[MQuad] Finding knee point for deltaBIC cutoff...')
        x, y, knee_x, cutoff = findKnee(self.df.deltaBIC)
        plt.plot(x, y)
        plt.axvline(x=knee_x, color="black", linestyle='--', label="Cutoff")
        plt.legend()
        plt.ylabel("log10(ΔBIC)")
        plt.xlabel("Cumulative probability")
        plt.savefig(os.path.join(out_dir, 'deltaBIC_cdf.pdf'))
        print('deltaBIC cutoff =', cutoff)

        self.sorted_df['PASS_KP'] = self.sorted_df.deltaBIC >= cutoff
        self.sorted_df['PASS_MINCELLS'] = self.sorted_df.num_cells_minor_cpt >= min_cells
        self.final_df = self.sorted_df[self.sorted_df.PASS_KP & self.sorted_df.PASS_MINCELLS]

    idx = self.final_df.index
    best_ad = self.ad[idx]
    best_dp = self.dp[idx]

    print(f'Number of variants passing threshold: {len(best_ad)}')

    if self.variants is not None:
        best_vars = np.array(self.variants)[idx]
        with open(os.path.join(out_dir, 'passed_variant_names.txt'), "w") as var_file:
            var_file.write('\n'.join(best_vars))

    if export_heatmap:
        af = best_ad / best_dp
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.title("Allele frequency of top variants")
        sns.heatmap(af, cmap="YlGnBu", yticklabels=best_vars)
        plt.yticks(rotation=0)
        plt.savefig(os.path.join(out_dir, 'top_variants_heatmap.pdf'))

    if export_mtx:
        mmwrite(os.path.join(out_dir, 'passed_ad.mtx'), sparse.csr_matrix(best_ad))
        mmwrite(os.path.join(out_dir, 'passed_dp.mtx'), sparse.csr_matrix(best_dp))

    return best_ad, best_dp

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mitochondrial informative mutation selection for MAESTER data.')
    parser.add_argument('--input_folder', required=True, help='Path to the cellsnp output folder.')
    parser.add_argument('--output_folder', required=True, help='Path to the output folder.')
    parser.add_argument('--num_variants', type=int, help='Number of variants to be selected.', default=None)
    parser.add_argument('--nproc', type=int, default=32, help='Number of processors to use.')
    parser.add_argument('--minDP', type=int, default=2, help='Minimum depth of coverage.')
    args = parser.parse_args()

    print("Start to read data...")
    vcf_file = os.path.join(args.input_folder, "cellSNP.cells.vcf.gz")
    cell_vcf = load_VCF(vcf_file, biallelic_only=True)
    print(f"Loaded VCF file: {vcf_file}")

    cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
    cell_dat['variants'] = cell_vcf['variants']

    out_dir = args.output_folder

    print("Start to select variants...")
    mdphd = Mquad(AD=cell_dat['AD'], DP=cell_dat['DP'], variant_names=cell_dat['variants'])
    df = mdphd.fit_deltaBIC(out_dir=out_dir, nproc=args.nproc, minDP=args.minDP)
    df.sort_values(by=['deltaBIC'], ascending=False).to_csv(os.path.join(out_dir, 'sorted_debug_BIC_params.csv'), index=False)

    best_ad, best_dp = select_informative_variants(self=mdphd, out_dir=out_dir, numb_cutoff=args.num_variants)
    print("The MQuad process is done!")
