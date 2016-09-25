"""
Read one or more RVCF files and write new VCF-like files with stability selection probabilities for each taxon.

test on one file:
  python \
      src/hominid/lasso/stability_selection_features_lasso_cv_C.py \
      ~/lab/glowing-happiness/hominid/project/hmp/16S_laf_sadj/mesabi/lasso_lars_cv_C/results/hmp_16S_laf_sadj_anterior_nares/results_mesabi_lasso_lars_cv_C_hmp_16S_laf_sadj_0_anterior_nares.rvcf

run on all body sites (don't forget the double quotes around the glob):
  python \
      src/hominid/lasso/stability_selection_features_lasso_cv_C.py \
      "~/lab/glowing-happiness/hominid/project/hmp/16S_laf_sadj/mesabi/lasso_lars_cv_C/results/hmp_16S_laf_sadj_*/results_mesabi_lasso_lars_cv_C_hmp_16S_laf_sadj_0_*.rvcf"

"""
import argparse
import glob
import os
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LassoLarsCV, RandomizedLasso
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.utils import ConvergenceWarning

from hominid.hominid import LassoMPI


def stability_selection_features_lasso_cv_C(
        rvcf_input_file_path,
        taxon_table_input_file_path,
        output_file_path,
        lo_alpha_coef,
        transform,
        maf_lower_cutoff,
        snp_limit
):
    # the input file path may be a glob
    rvcf_input_file_path = os.path.expanduser(rvcf_input_file_path)
    rvcf_input_file_path_list = glob.glob(rvcf_input_file_path)
    print('input file(s):\n\t{}'.format('\n\t'.join(rvcf_input_file_path_list)))
    print('body_site\tall_snps\ttested_snps\tretained_snps')

    if transform == 'none' or transform == 'None':
        transform = None

    for rvcf_input_file_path in sorted(rvcf_input_file_path_list):

        # rsq_median will be NA if a SNP did not meet the minimum criteria for testing
        #   greater than 50 samples
        #   greater than 5 homozygous and heterozygous samples
        # count() does not count NAs

        # use LassoMPI to help with aligning SNP values and taxon abundances
        lasso_mpi = LassoMPI(
            input_vcf_fp=rvcf_input_file_path,
            input_taxon_table_fp=taxon_table_input_file_path,
            output_vcf_fp='stability_selection_features_lasso_cv_C_dummy_output.txt',
            transform=transform,
            maf_lower_cutoff=maf_lower_cutoff,
            snp_limit=-1,
            cv_count=100,
            permutation_method='none'
        )
        lasso_mpi.initialize_controller()

        with open(output_file_path, 'w') as output_file:
            snp_count = 0
            snp_lasso_task_generator = lasso_mpi.get_task()
            for snp_lasso_task in snp_lasso_task_generator:
                if snp_lasso_task is None:
                    print('all done!')
                    break
                if snp_lasso_task.snp_with_rsq_df.rsq_pibsp_median_95ci_lo[0] > 0.0:
                    snp_count += 1
                    # found a SNP result for feature selection
                    print('  rsq_median: {:5.4f}'.format(float(snp_lasso_task.snp_with_rsq_df.rsq_median)))
                    t0 = time.time()
                    feature_scores_df, alphas = select_features(
                        snp_lasso_task.aligned_snp_df,
                        snp_lasso_task.aligned_taxa_df,
                        lo_alpha_coef
                    )
                    t1 = time.time()
                    print('time for feature selection: {:4.3f}s'.format(t1-t0))

                    feature_scores_df = feature_scores_df.reindex(lasso_mpi.taxon_table_df.index)
                    snp_df = pd.concat([snp_lasso_task.snp_with_rsq_df, feature_scores_df.transpose()], axis=1)
                    if snp_count == 1:
                        print('printing header')
                        snp_df.to_csv(
                            output_file,
                            index=False,
                            sep='\t',
                            na_rep=np.nan,
                            header=True
                        )
                    else:
                        snp_df.to_csv(
                            output_file,
                            index=False,
                            sep='\t',
                            na_rep=np.nan,
                            header=False
                        )

                if snp_count == snp_limit:
                    print('SNP limit reached')
                    break


def select_features(aligned_snp_df, aligned_taxa_df, lo_alpha_coef):
    X = aligned_taxa_df.values
    y = aligned_snp_df.values.flatten()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        warnings.simplefilter('ignore', ConvergenceWarning)
        lars_cv = LassoLarsCV(cv=StratifiedShuffleSplit(y, n_iter=100, test_size=0.2)).fit(X, y)

    print('lars_cv.alphas_: {}'.format(lars_cv.alphas_))
    alphas = np.linspace(lars_cv.alphas_[0], lo_alpha_coef * lars_cv.alphas_[0], 10)
    print('alphas: {}'.format(alphas))
    clf = RandomizedLasso(
        alpha=alphas,
        sample_fraction=0.8,
        n_resampling=1000
        #random_state=13
    ).fit(X, y)

    feature_scores_df = pd.DataFrame(clf.scores_, index=aligned_taxa_df.columns)

    return feature_scores_df, lars_cv.alphas_


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('rvcf_input_file_path')
    argparser.add_argument('taxon_table_input_file_path')
    argparser.add_argument('output_file_path')
    argparser.add_argument(
        'lo_alpha_coef',
        type=float
    )
    argparser.add_argument('transform')
    argparser.add_argument(
        '--maf-lower-cutoff',
        type=float,
        default=0.2
    )
    argparser.add_argument(
        'snp_limit',
        type=int
    )
    args = argparser.parse_args()
    print(args)
    stability_selection_features_lasso_cv_C(**vars(args))


if __name__ == '__main__':
    main()
