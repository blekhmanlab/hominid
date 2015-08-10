"""
lime.py

Python MPI program using LASSO regression to find associations between host genetics and microbiome.

This program is based on example 9 from https://github.com/jbornschein/mpi4py-examples.

This program requires these packages:
  SciPy
  scikit-learn
  scikit-bootstrap
  statsmodels
  mpi4py
  pandas

"""
import argparse
import collections
import os
import time
import traceback
import warnings

from mpi4py import MPI
import numpy as np
import pandas as pd
import scikits.bootstrap
import scipy.stats
import sklearn.cross_validation
import sklearn.linear_model
from statsmodels.robust.scale import mad

# Define MPI message tags
READY_ = 0
DONE_ = 1
EXIT_ = 2
START_ = 3
EXCEPTION_= 4

# Initializations and preliminaries
comm = MPI.COMM_WORLD # get MPI communicator object
size = comm.size # total number of processes
rank = comm.rank # rank of this process
status = MPI.Status() # get MPI status object


class SnpLassoTask(object):
    def __init__(self, aligned_snp_df, aligned_taxa_df, snp_with_rsq_df, cv_count, permutation_method):
        self.aligned_snp_df = aligned_snp_df
        self.aligned_taxa_df = aligned_taxa_df
        self.snp_with_rsq_df = snp_with_rsq_df
        self.cv_count = cv_count
        self.permutation_method = permutation_method
        self.cv_score_list = None

    def do(self):
        print('testing SNP {} {}'.format(self.snp_with_rsq_df.GENE[0], self.snp_with_rsq_df.ID[0]))

        if self.permutation_method == 'no_permutation':
            y_labels = self.aligned_snp_df.values.flatten()
        elif self.permutation_method == 'uniform_permutation':
            y_labels = self.uniform_snp_permutation()
        elif self.permutation_method == 'group_permutation':
            y_labels = self.group_snp_permutation()
        else:
            raise Exception('unknown permutation_method {}'.format(self.permutation_method))

        self.cv_score_list = self.score_cv(y_labels)
        if len(self.cv_score_list) > self.cv_count:
            print('{} {} {} {} has fewer scores than expected: {}'.format(
                self.snp_with_rsq_df.CHROM[0],
                self.snp_with_rsq_df.POS[0],
                self.snp_with_rsq_df.GENE[0],
                self.snp_with_rsq_df.ID[0],
                len(self.cv_score_list)
            ))
        validation_score_array = np.asarray(self.cv_score_list)

        (rsq_mean_pibs95ci_lo, rsq_mean_pibs95ci_hi) = bootstrap_ci_lo_hi(
            validation_score_array, alpha=0.05, method='pi'
        )
        (rsq_mean_pibs99ci_lo, rsq_mean_pibs99ci_hi) = bootstrap_ci_lo_hi(
            validation_score_array, alpha=0.01, method='pi'
        )

        rsq_median = np.median(validation_score_array)
        (rsq_median_pibs95ci_lo, rsq_median_pibs95ci_hi) = bootstrap_ci_lo_hi(
            validation_score_array,
            alpha=0.05,
            statistic=np.median,
            method='pi'
        )
        (rsq_median_pibs99ci_lo, rsq_median_pibs99ci_hi) = bootstrap_ci_lo_hi(
            validation_score_array,
            alpha=0.01,
            statistic=np.median,
            method='pi'
        )

        rsq_mean = np.mean(validation_score_array)

        self.snp_with_rsq_df.loc[0, 'rsq_mean'] = rsq_mean
        self.snp_with_rsq_df.loc[0, 'rsq_std'] = np.std(validation_score_array)
        self.snp_with_rsq_df.loc[0, 'rsq_sem'] = scipy.stats.sem(validation_score_array)
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_mean_95ci_lo'] = rsq_mean_pibs95ci_lo
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_mean_95ci_hi'] = rsq_mean_pibs95ci_hi
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_mean_99ci_lo'] = rsq_mean_pibs99ci_lo
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_mean_99ci_hi'] = rsq_mean_pibs99ci_hi
        self.snp_with_rsq_df.loc[0, 'rsq_median'] = rsq_median
        self.snp_with_rsq_df.loc[0, 'rsq_mad'] = mad(validation_score_array)
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_median_95ci_lo'] = rsq_median_pibs95ci_lo
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_median_95ci_hi'] = rsq_median_pibs95ci_hi
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_median_99ci_lo'] = rsq_median_pibs99ci_lo
        self.snp_with_rsq_df.loc[0, 'rsq_pibsp_median_99ci_hi'] = rsq_median_pibs99ci_hi
        self.snp_with_rsq_df.loc[0, 'cv_skewness'] = scipy.stats.skew(validation_score_array)
        self.snp_with_rsq_df.loc[0, 'cv_kurtosis'] = scipy.stats.kurtosis(validation_score_array)
        print('{} {} rsq_mean 95% (pi)  : {:6.4f} <-- {:6.4f} --> {:6.4f}'.format(
            self.snp_with_rsq_df.GENE[0], self.snp_with_rsq_df.ID[0],
            rsq_mean_pibs95ci_lo, rsq_mean, rsq_mean_pibs95ci_hi
        ))
        print('{} {} rsq_median 95% (pi): {:6.4f} <-- {:6.4f} --> {:6.4f}'.format(
            self.snp_with_rsq_df.GENE[0], self.snp_with_rsq_df.ID[0],
            rsq_median_pibs95ci_lo, rsq_median, rsq_median_pibs95ci_hi
        ))
        print('{} {} rsq_mean 99% (pi)  : {:6.4f} <-- {:6.4f} --> {:6.4f}'.format(
            self.snp_with_rsq_df.GENE[0], self.snp_with_rsq_df.ID[0],
            rsq_mean_pibs99ci_lo, rsq_mean, rsq_mean_pibs99ci_hi
        ))
        print('{} {} rsq_median 99% (pi): {:6.4f} <-- {:6.4f} --> {:6.4f}'.format(
            self.snp_with_rsq_df.GENE[0], self.snp_with_rsq_df.ID[0],
            rsq_median_pibs99ci_lo, rsq_median, rsq_median_pibs99ci_hi
        ))

    def uniform_snp_permutation(self):
        y_true = self.aligned_snp_df.values.flatten()
        permuted_y_true = np.copy(y_true)
        np.random.shuffle(permuted_y_true)
        return permuted_y_true

    def group_snp_permutation(self):
        taxa_groups_by_sex = self.aligned_taxa_df.groupby('Sex')
        if len(taxa_groups_by_sex) == 2:
            # everything is good
            pass
        else:
            raise Exception('unexpected number of groups: {}'.format(len(taxa_groups_by_sex)))

        permuted_group_list = []
        for name, group in taxa_groups_by_sex:
            group_snp_copy = np.copy(self.aligned_snp_df[group.index].values.flatten())
            np.random.shuffle(group_snp_copy)
            permuted_group_snp_sr = pd.DataFrame(group_snp_copy, index=group.index)
            permuted_group_list.append(permuted_group_snp_sr)
        permuted_aligned_snp = pd.concat(permuted_group_list, axis=0).transpose()
        reindexed_permuted_aligned_snp = permuted_aligned_snp.reindex_like(self.aligned_snp_df)
        return reindexed_permuted_aligned_snp.values.flatten()

    def score_cv(self, y_true):
        validation_score_list = []
        val_skf = sklearn.cross_validation.StratifiedShuffleSplit(
            y_true,
            n_iter=self.cv_count,
            test_size=0.2
        )
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sklearn.utils.ConvergenceWarning)
            for train, test in val_skf:
                skf = sklearn.cross_validation.StratifiedKFold(
                    y_true[train],
                    n_folds=5,
                    #shuffle=True  # already shuffling above
                )
                lasso_lars_cv = sklearn.linear_model.LassoLarsCV(cv=skf)
                model = lasso_lars_cv.fit(
                    self.aligned_taxa_df.values[train],
                    y_true[train]
                )
                score = model.score(
                    self.aligned_taxa_df.values[test],
                    y_true[test]
                )
                validation_score_list.append(score)

        return validation_score_list


def bootstrap_ci_lo_hi(validation_score_array, alpha=0.05, method='bca', statistic=np.mean):
    return scikits.bootstrap.ci(
        data=validation_score_array,
        statfunction=statistic,
        alpha=alpha,
        method=method
    )


class LassoMPI(object):
    def __init__(self, input_vcf_fp, input_taxon_table_fp, output_vcf_fp, permutation_method,
                 transform=None, snp_limit=-1, cv_count=100):
        self.input_vcf_fp = os.path.expanduser(input_vcf_fp)
        self.input_taxon_table_fp = os.path.expanduser(input_taxon_table_fp)
        self.output_vcf_fp = os.path.expanduser(output_vcf_fp)
        self.permutation_method = permutation_method

        self.transform = transform
        self.snp_limit = snp_limit
        self.cv_count = cv_count

        self.taxon_table_df = None
        self.output_line_count = 0
        self.output_file = None

        output_vcf_dir_path, output_vcf_name = os.path.split(self.output_vcf_fp)
        output_vcf_base_name, output_vcf_ext = os.path.splitext(output_vcf_name)
        self.output_cv_scores_fp = os.path.join(
            output_vcf_dir_path,
            output_vcf_base_name + '_cv_scores.txt'
        )
        self.output_cv_scores_file = None
        self.complete_snp_task_count = None

    def initialize_controller(self):
        self.taxon_table_df = read_taxon_file(
            self.input_taxon_table_fp,
            transform=self.transform
        )

    def initialize_worker(self):
        pass

    def go(self):
        if rank == 0:
            self.initialize_controller()
            self.complete_snp_task_count = 0
            with open(self.output_vcf_fp, 'w') as self.output_file, \
                 open(self.output_cv_scores_fp, 'w') as self.output_cv_scores_file:
                # this process is the controller
                # while there are still running workers wait for a work request
                num_workers = size - 1
                closed_workers = 0
                a_task_gen = self.get_task()
                # need a fake task that is not None to get started
                a_task = object()
                while closed_workers < num_workers:
                    msg = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                    source = status.Get_source()
                    tag = status.Get_tag()
                    print('[controller] recv message from worker {} with tag {}'.format(source, tag))
                    if tag == READY_ and a_task is None:
                        #if a_task is None there are no more SNPs to test
                        # and we should not call next(a_task_gen) again
                        print('[controller] sending exit message to worker {}'.format(source))
                        comm.send(None, dest=source, tag=EXIT_)
                    elif tag == READY_ and a_task is not None:
                        a_task = next(a_task_gen)
                        if a_task is None:
                            print('[controller] sending exit message to worker {}'.format(source))
                            comm.send(None, dest=source, tag=EXIT_)
                        else:
                            print('[controller] sending a task to worker {} with SNP {} {}'.format(
                                source,
                                a_task.snp_with_rsq_df.GENE[0],
                                a_task.snp_with_rsq_df.ID[0]
                            ))
                            comm.send(a_task, dest=source, tag=START_)
                    elif tag == DONE_:
                        # save the results
                        self.complete_snp_task_count += 1
                        self.task_complete(msg)
                        print('[controller] received a processed task from worker {}'.format(
                            source
                        ))
                    elif tag == EXIT_:
                        print('[controller] received exit message from worker {}'.format(source))
                        num_workers = num_workers - 1
                    elif tag == EXCEPTION_:
                        print('[controller] received exception message from worker {} with SNP {} {}'.format(
                            source,
                            msg.snp_with_rsq_df.GENE[0],
                            msg.snp_with_rsq_df.ID[0]
                        ))
                        self.task_failed(msg)
                        num_workers = num_workers - 1
                    else:
                        print('[controller] unrecognized message from source {} with tag {}:\n{}'.format(
                            source,
                            tag,
                            msg
                        ))
                print('[controller] all workers have exited')
        else:
            # this process is a worker
            self.initialize_worker()
            worker_t0 = time.time()
            task_count = 0
            name = MPI.Get_processor_name()
            print("[worker {}] running on {}".format(rank, name))
            while True:
                print('[worker {}] sending request for work'.format(rank))
                comm.send(None, dest=0, tag=READY_)
                print('[worker {}] waiting for work'.format(rank))
                t0 = time.time()
                my_task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                t1 = time.time()
                print('[worker {}] waited {:4.2f}s for work'.format(rank, t1-t0))
                tag = status.Get_tag()
                print('[worker {}] received message with tag {}'.format(rank, tag))
                if tag == START_:
                    try:
                        t0 = time.time()
                        my_task.do()
                        t1 = time.time()
                        print('[worker {}] task time: {:5.2f}s'.format(rank, t1-t0))
                        comm.send(my_task, dest=0, tag=DONE_)
                        task_count += 1
                    except Exception as e:
                        print('[worker {}] quitting with exception'.format(rank))
                        traceback.print_exc()
                        comm.send(my_task, dest=0, tag=EXCEPTION_)
                        break
                elif tag == EXIT_:
                    print('[worker {}] received exit message'.format(rank))
                    comm.send(None, dest=0, tag=EXIT_)
                    break
                else:
                    # what is this?
                    print('[worker {}] received an unrecognized message with tag {} - exiting'.format(rank, tag))
                    #raise Exception('worker {}: unexpected tag: {}'.format(rank, tag))
                    break

            worker_t1 = time.time()
            print('[worker {}] exiting after {} tasks in  {:6.2f}s'.format(rank, task_count, worker_t1-worker_t0))

    def get_task(self):
        """
        Return a SNP for testing.
        Reject SNPs that do not meet selection criteria.
        If there are no SNPs left return None.
        :return:
        """
        snp_task_count = 0
        vcf_reader = pd.read_csv(self.input_vcf_fp, sep='\t', chunksize=1)
        for snp_df in vcf_reader:
            # if self.snp_limit is -1 then all SNPs will be processed
            if self.snp_limit == snp_task_count:
                print('SNP limit {} has been reached'.format(snp_task_count))
                break

            aligned_snp_df, aligned_taxa_df = align_snp_and_taxa(
                snp_df,
                self.taxon_table_df
            )

            # check selection criteria
            new_snp_metadata_df, genotype_counter = self.get_new_snp_metadata_df(aligned_snp_df)
            if 'rsq_median' in snp_df.columns:
                print('snp_df already has metadata')
                snp_with_rsq_df = snp_df
            else:
                snp_with_rsq_df = pd.concat([new_snp_metadata_df, snp_df], axis=1)

            print('{} aligned samples for {} {} {}'.format(
                snp_with_rsq_df.aligned_sample_count[0],
                snp_df.CHROM[0],
                snp_df.GENE[0],
                snp_df.ID[0]
            ))

            snp_accepted = True

            if snp_with_rsq_df.aligned_sample_count[0] < 50:
                print('  fewer than 50 aligned samples for {} {}'.format(snp_df.GENE[0], snp_df.ID[0]))
                snp_accepted = False
            else:
                pass

            # greater than 5 of each genotype
            # note: aligned_snp_df dtypes should be 'int64'
            if any([not count > 5 for (genotype, count) in genotype_counter.most_common()]):
                print('  fewer than 6 samples for at least one genotype {} {}'.format(snp_df.GENE[0], snp_df.ID[0]))
                snp_accepted = False
            else:
                pass

            # maf >= 0.2
            if snp_with_rsq_df.maf[0] < 0.2:
                print('  maf {:3.2f} is too low for {} {}'.format(
                    snp_with_rsq_df.maf[0],
                    snp_df.GENE[0],
                    snp_df.ID[0]
                ))
                snp_accepted = False
            else:
                pass

            if snp_accepted:
                # this SNP has passed the selection criteria
                snp_task = SnpLassoTask(
                    aligned_snp_df=aligned_snp_df,
                    aligned_taxa_df=aligned_taxa_df,
                    snp_with_rsq_df=snp_with_rsq_df,
                    cv_count=self.cv_count,
                    permutation_method=self.permutation_method
                )
                snp_task_count += 1
                yield snp_task
            else:
                # write this SNP to the output file anyway
                # statistics will be NA
                # do not return anything, go to the next SNP in the input file
                self.write_snp_to_file(snp_with_rsq_df)

        # end of the input
        yield None

    def get_new_snp_metadata_df(self, aligned_snp_df):
        aligned_sample_count = aligned_snp_df.shape[1]
        genotype_counter = collections.Counter({0: 0, 1: 0, 2: 0})
        genotype_counter.update(aligned_snp_df.values[0])
        vaf = aligned_snp_df.sum(axis=0).sum() / (2.0 * aligned_snp_df.shape[1])
        maf = min(1.0 - vaf, vaf)
        new_snp_metadata_df = pd.DataFrame.from_items([
            ('aligned_sample_count', [aligned_sample_count]),
            ('aligned_count_0', [genotype_counter[0]]),
            ('aligned_count_1', [genotype_counter[1]]),
            ('aligned_count_2', [genotype_counter[2]]),
            ('vaf', [vaf]),
            ('maf', [maf]),
            # using np.nan gives these columns type float
            ('rsq_mean', [np.nan]),
            ('rsq_pibsp_mean_95ci_lo', [np.nan]),
            ('rsq_pibsp_mean_95ci_hi', [np.nan]),
            ('rsq_pibsp_mean_99ci_lo', [np.nan]),
            ('rsq_pibsp_mean_99ci_hi', [np.nan]),
            ('rsq_pibsp_median_95ci_lo', [np.nan]),
            ('rsq_pibsp_median_95ci_hi', [np.nan]),
            ('rsq_pibsp_median_99ci_lo', [np.nan]),
            ('rsq_pibsp_median_99ci_hi', [np.nan]),
            ('rsq_std', [np.nan]),
            ('rsq_sem', [np.nan]),
            ('rsq_median', [np.nan]),
            ('rsq_mad', [np.nan]),
            ('cv_skewness', [np.nan]),
            ('cv_kurtosis', [np.nan])
        ])
        return new_snp_metadata_df, genotype_counter

    def task_complete(self, snp_task):
        self.write_snp_to_file(snp_task.snp_with_rsq_df)
        self.write_cv_score_list_to_file(snp_task)

    def task_failed(self, snp_task):
        self.write_snp_to_file(snp_task.snp_with_rsq_df)

    def write_snp_to_file(self, snp_with_rsq_df):
        header = (self.output_line_count == 0)
        snp_with_rsq_df.to_csv(
            self.output_file,
            index=False,
            header=header,
            sep='\t',
            na_rep='NA',
            float_format='%6.4f'
        )
        self.output_line_count += 1

    def write_cv_score_list_to_file(self, snp_task):
        if self.complete_snp_task_count == 1:
            # write the column headers
            self.output_cv_scores_file.write('CHROM\tPOS\tID\tGENE\t')
            self.output_cv_scores_file.write(
                '\t'.join(['cv_s_{}'.format(n) for n in range(len(snp_task.cv_score_list))])
            )
            self.output_cv_scores_file.write('\n')
        self.output_cv_scores_file.write(snp_task.snp_with_rsq_df.CHROM[0])
        self.output_cv_scores_file.write('\t')
        self.output_cv_scores_file.write(str(snp_task.snp_with_rsq_df.POS[0]))
        self.output_cv_scores_file.write('\t')
        self.output_cv_scores_file.write(snp_task.snp_with_rsq_df.ID[0])
        self.output_cv_scores_file.write('\t')
        self.output_cv_scores_file.write(str(snp_task.snp_with_rsq_df.GENE[0]))
        self.output_cv_scores_file.write('\t')
        self.output_cv_scores_file.write('\t'.join(['{:9.6f}'.format(cv_score) for cv_score in snp_task.cv_score_list]))
        self.output_cv_scores_file.write('\n')


def read_taxon_file(taxon_file_path, transform=None):
    """
    Read a taxon table with taxa on the rows and samples on the columns.
    :param taxon_file_path:
    :param kwargs:
    :return: pandas.DataFrame
    """
    print('loading taxon table file {}'.format(taxon_file_path))
    if not os.path.exists(taxon_file_path):
        error_msg = 'file does not exist:\n  "{}"'.format(taxon_file_path)
        print(error_msg)
        raise Exception(error_msg)
    else:
        taxon_table = pd.read_csv(
            taxon_file_path,
            sep='\t',
            comment='#',
            index_col=0
        )
        print('  taxon table has {} rows'.format(len(taxon_table.index)))
        print('  taxon table has {} columns'.format(len(taxon_table.columns)))

        if transform is None:
            print('no transformation')
        elif transform == 'arcsinsqrt':
            print('applying arcsin sqrt transformation')
            taxon_table = taxon_table.apply(np.sqrt).apply(np.arcsin)
        elif transform == 'normalize':
            print('applying normalization transformation')
            taxon_table = taxon_table.div(taxon_table.sum())
        else:
            raise Exception('unrecognized transform {}'.format(transform))

        return taxon_table


def align_snp_and_taxa(snp_df, taxon_table_df):
    # use taxon_table_df column headers to exclude metadata from snp_df

    # snp_df looks like:
    #   <first 9 columns> '1234' '2345' '3456' '4567'
    #   . . . . . . . . .     0      1      2     NA
    # taxon_table_df looks like:
    #            '1234' '2345' '4567'
    #   taxon_0    0.1    0.2    0.4
    #   taxon_1    0.1    0.2    0.4
    #   taxon_2    0.1    0.2    0.4
    #   taxon_3    0.1    0.2    0.4

    # snp_aligned_to_taxa_df looks like:
    #   '1234' '2345' '4567'
    #       0      1     NA
    snp_aligned_to_taxa_df = snp_df[taxon_table_df.columns]

    # snp_aligned_to_taxa_dropna_df looks like:
    #   '1234' '2345'
    #       0      1
    snp_aligned_to_taxa_dropna_df = snp_aligned_to_taxa_df.dropna(axis=1)

    # taxa_aligned_to_snp_df looks like:
    #            '1234' '2345'
    #   taxon_0    0.1    0.2
    #   taxon_1    0.1    0.2
    #   taxon_2    0.1    0.2
    #   taxon_3    0.1    0.2
    taxa_aligned_to_snp_df = taxon_table_df[snp_aligned_to_taxa_dropna_df.columns]

    # taxa_aligned_to_snp_df is taxon-by-subject, or feature-by-sample
    # we need sample-by-feature for sklearn so return the transpose
    return snp_aligned_to_taxa_dropna_df, taxa_aligned_to_snp_df.T


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        'input_taxon_table_fp'
    )
    arg_parser.add_argument(
        'input_vcf_fp'
    )
    arg_parser.add_argument(
        'output_vcf_fp'
    )
    arg_parser.add_argument(
        'transform'
    )
    arg_parser.add_argument(
        'snp_limit',
        type=int
    )
    arg_parser.add_argument(
        'cv_count',
        type=int
    )
    arg_parser.add_argument(
        'permutation_method',
        type=str
    )

    args = arg_parser.parse_args()
    print(args)

    lasso_mpi = LassoMPI(
        **vars(args)
    )
    lasso_mpi.go()
