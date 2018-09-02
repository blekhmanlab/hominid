"""
Read a rvcf file with stability selection scores for taxa.
Sort the dataframe by rsq_median.
Save box and bar plots for the top N SNPs.

Note: R must have packages dplyr and ggplot2 installed.

I installed rpy2 with conda, which installs R in the virtual environment directory.
In my case that is ~/miniconda3/envs/hve/lib/R. This takes care of everything needed to
run this script.

conda install r-essentials

This installs R and packages such as ggplot2 and dplyr.

When using Canopy or Anaconda python distributions it may happen that
R and rpy2 disagree about where to find gfortran and cause this:

  version `GFORTRAN_1.4' not found (required by /usr/lib/liblapack.so.3)

The solution in Ubuntu and related distributions is to point at the
system libgfortran with LD_PRELOAD like this:

  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0 python box_bar_plot.py <command line arguments>

usage:
    (16S)
    python box_bar_plot_lasso_lars_cv_C_stability_selection_features.py \
        ~/lab/glowing-happiness/hominid/project/hmp/16S_laf_sadj/mesabi/lasso_lars_c_C/results/hmp_16S_laf_sadj_anterior_nares/selected_features_results_mesabi_lasso_lars_cv_C_hmp_16S_laf_sadj_0_anterior_nares.rvcf \
        taxon table file path
        transform
        10 <- maximum SNP count

    (MGS modules)
    python box_bar_plot_lasso_lars_cv_C_stability_selection_features.py \
        ~/lab/glowing-happiness/hominid/project/hmp/MGS_humann_mpm_sadj/mesabi/lasso_lars_c_C/results/MGS_humann_mpm_sadj_anterior_nares/selected_features_results_mesabi_lasso_lars_cv_C_MGS_humann_mpm_sadj_0_anterior_nares.rvcf \
        taxon table file path
        transform
        10

"""
import argparse
import os
import shutil

import pandas as pd

# the import of readline is a work-around for a
# compatibility issue between anaconda and rpy2
# see https://github.com/ContinuumIO/anaconda-issues/issues/152
##import readline
import rpy2.robjects
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
base = importr('base')
grdevices = importr('grDevices')

from hominid.hominid import read_taxon_file, align_snp_and_taxa


def get_taxon_abundance_box_plot():
    box_plot_fnc = """
        require("dplyr")
        require("ggplot2")

        taxon_abundance_box_plot <- function(data, plot_file_path, title, xlabel, ylabel) {

            temp <- data[order(data$variant_allele_count),] #sort by variant_allele_count
            temp$genotype <- factor(temp$genotype,levels=unique(temp$genotype)) #use reordered genotypes as levels

            pdf(plot_file_path)

            ap <- ggplot(data=temp,
                         aes(x=genotype,y=abundance)
                  )
            ap <- ap + geom_boxplot()
            ap <- ap + ggtitle(title)
            ap <- ap + labs(x=xlabel, y=ylabel)
            ap <- ap + geom_jitter(position=position_jitter(w=0.1))

            print(ap)
            dev.off()
        }
        """
    pck = SignatureTranslatedAnonymousPackage(box_plot_fnc, 'pck')
    return pck.taxon_abundance_box_plot


def direct_taxon_abundance_box_plot(data, plot_file_path, title, xlabel, ylabel):
    grdevices.pdf(file=plot_file_path)

    gp = ggplot2.ggplot(data)
    pp = gp \
        + ggplot2.aes_string(x='genotype', y='abundance') \
        + ggplot2.geom_boxplot() \
        + ggplot2.ggtitle(title) \
        + ggplot2.labs(x=xlabel, y=ylabel) \
        + ggplot2.geom_jitter(position=ggplot2.position_jitter(w=0.1)) \
        + ggplot2.geom_point()

    pp.plot()

    grdevices.dev_off()


def get_taxon_abundance_stacked_bar_plot():
    box_plot_fnc = """
        require("dplyr")
        require("ggplot2")

        taxon_abundance_stacked_bar_plot <- function(data, plot_file_path, title, xlabel, ylabel) {

            temp <- data[order(data$variant_allele_count),] #sort by variant_allele_count
            temp$genotype <- factor(temp$genotype,levels=unique(temp$genotype)) #use reordered genotypes as levels

            #creates a new data frame with median abundance from each combo
            result <- temp %>%
                group_by(genotype, gene, taxon) %>%
                    summarize(medianAbundance = median(abundance))
            #If you want the heights of the bars to represent values in the data,
            #use stat="identity" and map a value to the y aesthetic.

            pdf(plot_file_path, width=8, height=4)

            ap <- ggplot(data=result, aes(x=genotype,y=medianAbundance,fill=taxon)) +
                geom_bar(stat='identity') +
                ggtitle(title)
            ap <- ap + labs(x=xlabel, y=ylabel)
            ap <- ap + theme(legend.direction = 'vertical', legend.position = 'bottom')
            ap <- ap + guides(fill = guide_legend(reverse = TRUE))
            print(ap)
            dev.off()
        }
        """
    pck = SignatureTranslatedAnonymousPackage(box_plot_fnc, 'pck')
    return pck.taxon_abundance_stacked_bar_plot


#def direct_abundance_stacked_bar_plot(data, plot_file_path, title, xlabel, ylabel):
#    pass


def box_bar_lasso_lars_cv_C_stability_selection_features(
        rvcf_input_file_path, taxon_table_file_path, transform, plot_output_dir_path, stability_cutoff, snp_count):
    print('plotting {} SNPs from {}'.format(snp_count, rvcf_input_file_path))

    if os.path.exists(plot_output_dir_path):
        # delete it
        print('deleting old plots')
        shutil.rmtree(plot_output_dir_path)
    os.makedirs(plot_output_dir_path)

    # read the rvcf file and sort by rsq_median
    df = pd.read_csv(rvcf_input_file_path, sep='\t', dtype={'CHROM': str})

    sorted_rsq_best_medians_df = df.sort_values(by='rsq_median', ascending=False)

    taxon_table_df = read_taxon_file(taxon_table_file_path, transform=transform)

    # these are proxies for R functions
    taxon_abundance_box_plot = get_taxon_abundance_box_plot()
    taxon_abundance_stacked_bar_plot = get_taxon_abundance_stacked_bar_plot()

    for row_i in range(sorted_rsq_best_medians_df.shape[0]):
        if row_i >= snp_count:
            break
        else:
            # get a 1-row dataframe
            snp_df = sorted_rsq_best_medians_df.iloc[[row_i]]
            aligned_snp_df, aligned_taxa_df = align_snp_and_taxa(
                snp_df,
                taxon_table_df
            )
            # get the taxon stability selection scores
            # use the taxon table df index to get column names for snp_df
            taxon_scores_df = snp_df.loc[:, taxon_table_df.index].transpose()
            sorted_taxon_scores_df = taxon_scores_df.sort_values(
                by=taxon_scores_df.columns[0], ascending=False)
            # print all sorted taxon scores to verify they are sorted high to low
            ##print('sorted_taxon_scores_df:\n{}'.format(sorted_taxon_scores_df))
            p_df_list = []
            summary_line = '{}\t{}\t'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID)
            for i, (selected_taxon, selected_taxon_row) in enumerate(sorted_taxon_scores_df.iterrows()):
                # use selected_taxon_row.index[0] to index the first and only column
                print('selected_taxon_row:\n{}'.format(selected_taxon_row))
                selected_taxon_score = selected_taxon_row.iloc[0]
                if selected_taxon_score >= stability_cutoff:
                    # trim 'Root;' from the front of the taxon name
                    if selected_taxon.startswith('Root;'):
                        taxon_name = selected_taxon[5:]
                    else:
                        taxon_name = selected_taxon
                    summary_line += '{}, '.format(taxon_name)
                    # print a box plot
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
                        'gt': [gts[int(v)] for v in aligned_snp_value_list]
                    })
                    p_df_list.append(p_df)
                    r_df = rpy2.robjects.vectors.DataFrame({
                        'abundance': rpy2.robjects.FloatVector(aligned_taxa_df[selected_taxon].values.tolist()),
                        'variant_allele_count': rpy2.robjects.StrVector([str(int(v)) for v in aligned_snp_value_list]),
                        'genotype': rpy2.robjects.StrVector([gts[int(v)] for v in aligned_snp_value_list])
                    })
                    print(taxon_name)
                    print('r_df:\n'.format(r_df))
                    taxon_abundance_box_plot(
                        r_df,
                        r_pdf_file_path,
                        '{} (score: {:4.3f})'.format(snp_df.iloc[0].GENE, selected_taxon_score),
                        '{} {}'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID),
                        selected_taxon
                    )
                else:
                    # the selected_taxon_score is below the cutoff or is nan
                    break

            # write a summary line and
            print(summary_line[:-2])
            #summary_file.write(summary_line[:-2])
            #summary_file.write('\n')
            # save a stacked bar plot
            if len(p_df_list) > 0:
                file_name = 'stacked_bar_plot_selected_taxa_{}_{}.pdf'.format(
                    snp_df.iloc[0].GENE,
                    snp_df.iloc[0].ID
                )
                stacked_bar_plot_file_path = os.path.join(plot_output_dir_path, file_name)
                p_df = pd.concat(p_df_list, axis=0)
                # at this point the index for p_df looks like
                #   0...76.0...76.0...76
                # replace the index
                p_df.index = range(p_df.shape[0])
                r_all_df = rpy2.robjects.vectors.DataFrame({
                    'abundance': rpy2.robjects.FloatVector(p_df['abundance'].values.tolist()),
                    'variant_allele_count': rpy2.robjects.StrVector([str(int(v)) for v in p_df['variant_allele_count'].values]),
                    'taxon': rpy2.robjects.StrVector(p_df['taxon']),
                    'gene': rpy2.robjects.StrVector(p_df['gene']),
                    'genotype': rpy2.robjects.StrVector(p_df['gt'])
                })
                stacked_bar_title = '{}\n{}'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID)
                taxon_abundance_stacked_bar_plot(
                    r_all_df,
                    stacked_bar_plot_file_path,
                    stacked_bar_title,
                    '{} {}'.format(snp_df.iloc[0].GENE, snp_df.iloc[0].ID),
                    'median abundance'
                )


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('rvcf_input_file_path')
    argparser.add_argument('taxon_table_file_path')
    argparser.add_argument('transform')
    argparser.add_argument('plot_output_dir_path')
    argparser.add_argument(
        'stability_cutoff',
        type=float
    )
    argparser.add_argument(
        'snp_count',
        type=int
    )
    args = argparser.parse_args()
    print(args)
    box_bar_lasso_lars_cv_C_stability_selection_features(**vars(args))


if __name__ == '__main__':
    main()
