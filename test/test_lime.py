"""

"""
import logging
import os

from lime.lime import LassoMPI

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name=__name__)


def test_cmd_line_1(tmpdir, request):
    """
    Test using example data.
    :param request: pytest fixture with information about the path to this test file
    :return:
    """

    test_lime_fp = request.module.__file__
    log.info('test_lime_fp: {}'.format(test_lime_fp))
    test_lime_dir_path, _ = os.path.splitext(test_lime_fp)

    output_vcf_fp = os.path.join(str(tmpdir), 'lime_example_output.rvcf')

    lasso_mpi = LassoMPI(
        input_vcf_fp=os.path.join(test_lime_dir_path, 'lime_example_snp_input.rvcf'),
        input_taxon_table_fp=os.path.join(test_lime_dir_path, 'lime_example_taxon_table_input.txt'),
        output_vcf_fp=output_vcf_fp,
        permutation_method='no_permutation',
        transform='arcsinsqrt',
        snp_limit=-1,
        cv_count=100
    )
    lasso_mpi.go()

    assert os.path.exists(output_vcf_fp)
