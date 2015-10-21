"""
Read a rvcf file with stability selection scores for taxa.
Sort the dataframe by rsq_median.
Print results.

usage:
    python sort_results.py \
        ../example/stability_selection_example_output.vcf \
        ../example/lime_example_taxon_table_input.txt \
        arcsinsqrt \
        0.5 \
        10
"""
import argparse
import sys

import pandas as pd

from lime import read_taxon_file, align_snp_and_taxa


def sort_results(rvcf_input_file_path, taxon_table_file_path, transform,
                 r_sqr_median_cutoff, stability_cutoff, snp_count, no_tables):
    print('plotting {} SNPs from {}'.format(snp_count, rvcf_input_file_path))

    # read the rvcf file and sort by rsq_median
    df = pd.read_csv(rvcf_input_file_path, sep='\t')
    #print('df.shape: {}'.format(df.shape))

    sorted_rsq_best_medians_df = df.sort_values(by='rsq_median', ascending=False)

    x_df = sorted_rsq_best_medians_df[sorted_rsq_best_medians_df.rsq_median > r_sqr_median_cutoff]
    print('{} SNPs with r_sqr > {:5.3f}'.format(x_df.shape[0], r_sqr_median_cutoff))

    #print('df:\n{}'.format(sorted_rsq_best_medians_df.rsq_median.head()))
    #print('df.rsq_median type: {}'.format(sorted_rsq_best_medians_df.rsq_median.dtype))

    taxon_table_df = read_taxon_file(taxon_table_file_path, transform=transform)

    for row_i in range(sorted_rsq_best_medians_df.shape[0]):
        if row_i >= snp_count:
            break
        else:
            # get a 1-row dataframe
            snp_df = sorted_rsq_best_medians_df.iloc[[row_i]]
            #print('snp_df.shape {}'.format(snp_df.shape))
            #print('row {}: rsq_median: {}'.format(row_i, snp_df.iloc[0].rsq_median))
            aligned_snp_df, aligned_taxa_df = align_snp_and_taxa(
                snp_df,
                taxon_table_df
            )
            #print('aligned_snp_df.shape {}'.format(aligned_snp_df.shape))
            #print('aligned_taxa_df.shape {}'.format(aligned_taxa_df.shape))
            # get the taxon stability selection scores
            # use the taxon table df index to get column names for snp_df
            taxon_scores_df = snp_df.loc[:, taxon_table_df.index].transpose()
            sorted_taxon_scores_df = taxon_scores_df.sort_values(by=taxon_scores_df.columns[0], ascending=False)
            # print all sorted taxon scores to verify they are sorted high to low
            ##print(sorted_taxon_scores_df)

            p_df_list = []
            print('{} {} {:5.3f}'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID, snp_df.iloc[0].rsq_median))
            summary_line = '{}\t{}\t'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID)
            for i, (selected_taxon, selected_taxon_row) in enumerate(sorted_taxon_scores_df.iterrows()):
                # use selected_taxon_row.index[0] to index the first and only column
                selected_taxon_score = selected_taxon_row.iloc[0]
                #print('taxon: {} score: {}'.format(selected_taxon, selected_taxon_score))
                if selected_taxon_score < stability_cutoff:
                    #print('done with selected taxa')
                    break
                else:
                    # trim 'Root;' from the front of the taxon name
                    if selected_taxon.startswith('Root;'):
                        taxon_name = selected_taxon[5:]
                    else:
                        taxon_name = selected_taxon
                    print('  {:5.3f} {}'.format(selected_taxon_score, taxon_name))
                    summary_line += '{}, '.format(taxon_name)
                    # print a box plot
                    '''
                    r_pdf_file_path = \
                        os.path.join(
                            plot_output_dir_path,
                            'best_taxa_{}_{}_{}_boxplot_{}.pdf'.format(
                                row_i,
                                snp_df.iloc[0].GENE,
                                snp_df.iloc[0].ID,
                                i
                            )
                        )
                    '''
                    #print('writing file {}'.format(r_pdf_file_path))
                    gts = [
                        snp_df.iloc[0].REF + snp_df.iloc[0].REF,  # 0
                        snp_df.iloc[0].REF + snp_df.iloc[0].ALT,  # 1
                        snp_df.iloc[0].ALT + snp_df.iloc[0].ALT   # 2
                    ]
                    aligned_snp_value_list = aligned_snp_df.values.flatten().tolist()
                    p_df = pd.DataFrame({
                        'chromosome': [snp_df.iloc[0].CHROM] * aligned_snp_df.shape[1],
                        'snp_id': [snp_df.iloc[0].ID] * aligned_snp_df.shape[1],
                        'gene': [snp_df.iloc[0].GENE] * aligned_snp_df.shape[1],
                        'taxon': [selected_taxon] * aligned_snp_df.shape[1],
                        'abundance': aligned_taxa_df[selected_taxon].values.tolist(),
                        'variant_allele_count': [str(int(v)) for v in aligned_snp_value_list],
                        'genotype': [gts[int(v)] for v in aligned_snp_value_list]
                    })
                    p_df_list.append(p_df)
                    '''
                    r_df = rpy2.robjects.vectors.DataFrame({
                        'abundance': rpy2.robjects.FloatVector(aligned_taxa_df[selected_taxon].values.tolist()),
                        'variant_allele_count': rpy2.robjects.StrVector([str(int(v)) for v in aligned_snp_value_list]),
                        'genotype': rpy2.robjects.StrVector([gts[int(v)] for v in aligned_snp_value_list])
                    })
                    '''
                    if no_tables:
                        pass
                    else:
                        #print(taxon_name)
                        #print(snp_df.iloc[0].GENE)
                        #print(snp_df.iloc[0].ID)
                        p_df[['abundance', 'variant_allele_count', 'genotype']].to_csv(
                            sys.stdout,
                            sep='\t'
                        )
                        #print(p_df[['abundance', 'variant_allele_count', 'gt']])
            # write a summary line and
            #print(summary_line[:-2])
            # save a stacked bar plot
            if len(p_df_list) > 0:
                file_name = 'stacked_bar_plot_selected_taxa_{}_{}.pdf'.format(
                    snp_df.iloc[0].GENE,
                    snp_df.iloc[0].ID
                )
                #stacked_bar_plot_file_path = os.path.join(plot_output_dir_path, file_name)
                p_df = pd.concat(p_df_list, axis=0)
                # at this point the index for p_df looks like
                #   0...76.0...76.0...76
                # replace the index
                p_df.index = range(p_df.shape[0])
                #p_df.to_csv(file_path, sep='\t')
                '''
                r_all_df = rpy2.robjects.vectors.DataFrame({
                    'abundance': rpy2.robjects.FloatVector(p_df['abundance'].values.tolist()),
                    'variant_allele_count': rpy2.robjects.StrVector([str(int(v)) for v in p_df['variant_allele_count'].values]),
                    'taxon': rpy2.robjects.StrVector(p_df['taxon']),
                    'gene': rpy2.robjects.StrVector(p_df['gene']),
                    'genotype': rpy2.robjects.StrVector(p_df['gt'])
                })
                '''
                stacked_bar_title = '{}\n{}'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('rvcf_input_file_path')
    argparser.add_argument('taxon_table_file_path')
    argparser.add_argument('transform')
    argparser.add_argument(
        'r_sqr_median_cutoff',
        type=float
    )
    argparser.add_argument(
        'stability_cutoff',
        type=float
    )
    argparser.add_argument(
        'snp_count',
        type=int
    )
    argparser.add_argument(
        '--no-tables',
        action='store_true'
    )
    args = argparser.parse_args()
    print(args)
    sort_results(**vars(args))
