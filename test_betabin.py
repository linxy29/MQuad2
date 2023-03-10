# Xinyi Lin, 202303
# Goal: Using the betabinomial model to select informative mutations (MQuad2)
# Useage: python /home/linxy29/code/MQuad2/test_betabin.py "day7"

import vireoSNP
from vireoSNP.utils.vcf_utils import read_sparse_GeneINFO, load_VCF, write_VCF, parse_donor_GPb
import mquad
from scipy.io import mmread
import pickle
import argparse

parser = argparse.ArgumentParser(description='Mitochondrial informative mutation selection for MAESTER data.')
parser.add_argument('input_folder', help='Path of the cellsnp output folder.')
parser.add_argument('output_folder', help='Path of the output folder.')
args = parser.parse_args()

# sample = "gct86"
# vcf_file = "/home/linxy29/data/maester/oagct/" + sample + "/maester_cellSNP/cellSNP.cells.vcf.gz"
vcf_file = args.input_folder + "cellSNP.cells.vcf.gz"
cell_vcf = load_VCF(vcf_file, biallelic_only=True)
print("Loaded VCF file: %s" % vcf_file)
cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
cell_dat['variants'] = cell_vcf['variants']

## setting for mquad
out_dir = args.output_folder
nproc = 32
minDP = 2
beta_mode = True
minCell = 10
cutoff = None

from mquad.mquad import Mquad
mdphd = Mquad(AD = cell_dat['AD'], DP = cell_dat['DP'], variant_names = cell_dat['variants'])
df = mdphd.fit_deltaBIC(out_dir = out_dir, nproc = nproc, minDP = minDP, beta_mode = True)
best_ad, best_dp = mdphd.selectInformativeVariants(min_cells = minCell, out_dir = out_dir, tenx_cutoff=cutoff)
print("Done!")
