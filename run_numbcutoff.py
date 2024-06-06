# Xinyi Lin, 202303
# Goal: Using the betabinomial model to select informative mutations (MQuad2) and provide pre-defined number of variants
# Useage: python /home/linxy29/code/MQuad2/run_numbcutoff.py ./maester_assemble_trimmed_aligned_mt_addtag_cellSNP0 ./maester_assemble_trimmed_aligned_mt_addtag_cellSNP0_betamquad_2000variant 2000

import mquad
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix
import numpy as np
import mquad
from mquad.mquad import Mquad
import pandas as pd
import os
from os import path
import matplotlib.pyplot as plt
from mquad.mquad_utils import findKnee
from scipy import sparse
import seaborn as sns
import argparse
from vireoSNP.utils.vcf_utils import read_sparse_GeneINFO, load_VCF

def selectInformativeVariants(self, min_cells=2, export_heatmap=True, export_mtx=True, out_dir=None, existing_df=None, tenx_cutoff=None, numb_cutoff = None):
        #takes self.df, return best_ad and best_dp as array

        if existing_df is not None:
            #input /path/to/unsorted_debug_BIC_params.csv for existing df if model is already fit
            print('[MQuad] Fitted model detected, using' + existing_df + '...')
            self.df = pd.read_csv(existing_df)
            self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        if out_dir is not None:
            if path.exists(out_dir) is not True:
                try:
                    os.mkdir(out_dir)
                except:
                    print("[MQuad] Can't make directory, do you have permission?")
        else:
            print('[MQuad] Out directory already exists, overwriting content inside...')
        
        if numb_cutoff > self.df.shape[0]:
            print('[MQuad] Number of variants to select exceeds number of variants in dataset, turning off numb mode...')
            numb_cutoff = None
        
        if numb_cutoff > 50:
            print('[MQuad] Number of variants to select exceeds 50, turning off plot heatmap...')
            export_heatmap = False
        
        if tenx_cutoff is not None and numb_cutoff is not None:
            print('[MQuad] Both tenx and numb mode used, turning off numb mode...')
            numb_cutoff = None

        if tenx_cutoff is not None:
            print('[MQuad] Tenx mode used with cutoff = ' + str(tenx_cutoff))
            self.final_df = self.sorted_df[self.sorted_df.deltaBIC >= float(tenx_cutoff)]
            self.final_df = self.final_df[self.sorted_df.num_cells_minor_cpt >= min_cells]
        elif numb_cutoff is not None:
            print('[MQuad] Numb mode used with cutoff = ' + str(numb_cutoff))
            ## select numb_cutoff number of variants
            self.final_df = self.sorted_df.head(numb_cutoff)
            self.final_df = self.final_df[self.sorted_df.num_cells_minor_cpt >= min_cells]
        else:
            print('[MQuad] Finding knee point for deltaBIC cutoff...')
            #self.filt_df = self.sorted_df[self.sorted_df.deltaBIC >= 10]
            x,y,knee_x, cutoff = findKnee(self.df.deltaBIC)
            plt.plot(x, y)
            plt.axvline(x=knee_x, color="black", linestyle='--',label="cutoff")
            plt.legend()
            plt.ylabel("log10(\u0394BIC)")
            plt.xlabel("Cumulative probability")
            plt.savefig(out_dir + '/' + 'deltaBIC_cdf.pdf')

            print('deltaBIC cutoff = ', cutoff)
            #self.sorted_df['VALID'] = self.validateSNP(self.sorted_df.variant_name)
            self.sorted_df['PASS_KP'] = self.sorted_df.deltaBIC.apply(lambda x: True if x >= cutoff else False)
            self.sorted_df['PASS_MINCELLS'] = self.sorted_df.num_cells_minor_cpt.apply(lambda x: True if x >= min_cells else False)

            self.final_df = self.sorted_df[(self.sorted_df.PASS_KP == True) & (self.sorted_df.PASS_MINCELLS == True)]


        idx = self.final_df.index
        best_ad = self.ad[idx]
        best_dp = self.dp[idx]

        print('Number of variants passing threshold: '  + str(len(best_ad)))

        #fname = by + '_' + str(threshold) + '_'

        if self.variants is not None:
            best_vars = np.array(self.variants)[idx]
            #renamed_vars = []
            #for var in best_vars:
                #renamed_vars.append((var.split('_')[1] + var.split('_')[2] + '>' + var.split('_')[3]))

            with open(out_dir + '/' + 'passed_variant_names.txt', "w+") as var_file:
                #var_file.write('\n'.join(str(var) for var in renamed_vars))
                var_file.write('\n'.join(str(var) for var in best_vars))
                
        if export_heatmap:
            af = best_ad/best_dp
            #af = af.fillna(0)
            fig, ax = plt.subplots(figsize=(8,6))
            plt.title("Allele frequency of top variants")
            plt.style.use('seaborn-dark')
            pal = "YlGnBu"
            if self.variants is not None:
                sns.heatmap(af, cmap=pal, yticklabels=best_vars)
                #sns.heatmap(af, cmap=pal, yticklabels=renamed_vars)
                plt.yticks(rotation=0)
            else:
                sns.heatmap(af, cmap=pal)
                plt.yticks(rotation=0)
            plt.savefig(out_dir + '/' + 'top variants heatmap.pdf')

        #export ad dp mtx out for vireo
        if export_mtx is True:
            mmwrite(out_dir + '/' + 'passed_ad.mtx', sparse.csr_matrix(best_ad))
            mmwrite(out_dir + '/' + 'passed_dp.mtx', sparse.csr_matrix(best_dp))

        return best_ad, best_dp

parser = argparse.ArgumentParser(description='Mitochondrial informative mutation selection for MAESTER data.')
parser.add_argument('input_folder', help='Path of the cellsnp output folder.')
parser.add_argument('output_folder', help='Path of the output folder.')
parser.add_argument('num_variants', type=int, help='Number of variants to be selected.', default=None)
args = parser.parse_args()

print("Start to read data ...")
vcf_file = args.input_folder + "/cellSNP.cells.vcf.gz"
cell_vcf = load_VCF(vcf_file, biallelic_only=True)
print("Loaded VCF file: %s" % vcf_file)
cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
cell_dat['variants'] = cell_vcf['variants']

## setting for mquad
out_dir = args.output_folder
nproc = 32
minDP = 2

print("Start to select variants ...")
mdphd = Mquad(AD = cell_dat['AD'], DP = cell_dat['DP'], variant_names = cell_dat['variants'])
df = mdphd.fit_deltaBIC(out_dir = out_dir, nproc = nproc, minDP = minDP, beta_mode = True)
df.sort_values(by=['deltaBIC'], ascending=False).to_csv(out_dir + '/' + 'sorted_debug_BIC_params.csv', index=False)
best_ad, best_dp = selectInformativeVariants(self=mdphd, out_dir = out_dir, numb_cutoff = args.num_variants)
print("The MQuad process is done!")
