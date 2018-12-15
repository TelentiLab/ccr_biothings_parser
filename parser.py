import os
import logging
from biothings.utils.dataload import dict_sweep

FILE_NOT_FOUND_ERROR = 'Cannot find input file: {}'  # error message constant
FILE_LINES = 8188410  # autosome
# FILE_LINES = 171987  # x chrom

# change following parameters accordingly
source_name = 'ccr'  # source name that appears in the api response
# file_name = 'ccrs.xchrom.v2.20180420.bed'  # xchrom file
file_name = 'ccrs.autosomes.v2.20180420.bed'  # autosome file
delimiter = '\t'  # the delimiter that separates each field

# configure logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', level=logging.INFO)
logger = logging.getLogger('ccr_logger')


def load_data(data_folder: str):
    """
    Load data from a specified file path. Parse each line into a dictionary according to the schema
    given by `data_schema`. Then process each dict by normalizing data format, remove null fields (optional).
    Append each dict into final result using its id.

    :param data_folder: the path(folder) where the data file is stored
    :return: a generator that yields data.
    """
    input_file = os.path.join(data_folder, file_name)
    # raise an error if file not found
    if not os.path.exists(input_file):
        logger.error(FILE_NOT_FOUND_ERROR.format(input_file))
        raise FileExistsError(FILE_NOT_FOUND_ERROR.format(input_file))

    with open(input_file, 'r') as file:
        logger.info(f'start reading file: {file_name}')
        count = 0
        for line in file:   # read and parse each line into a dict
            logger.info(f'reading line {count} ({(count / FILE_LINES * 100):.2f}%)')  # format to use 2 decimals
            count += 1
            # schema: (chrom, start, end, ccr_pct, gene, ranges, varflag, syn_density,
            #          cpg, cov_score, resid, resid_pctile, unique_key)
            try:
                (chrom, start, end, ccr_pct, gene, ranges, varflag, syn_density, cpg,
                    cov_score, resid, resid_pctile, unique_key) = line.strip().split(delimiter)
            except ValueError:
                logger.error(f'failed to unpack line {count}: {line}')
                continue    # skip error line
            _id = f'chr{chrom}:g.{start}_{end}'
            # enforce data type
            try:
                variant = {
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'ccr_pct': float(ccr_pct),
                    'gene': gene,
                    'ranges': ranges.split(','),
                    'varflag': varflag is 'VARTRUE',
                    'syn_density': float(syn_density),
                    'cpg': float(cpg),
                    'cov_score': float(cov_score),
                    'resid': float(resid),
                    'resid_pctile': float(resid_pctile),
                    'unique_key': int(unique_key),
                }
            except ValueError as e:
                logger.error(f'failed to cast type for line {count}: {e}')
                continue  # skip error line

            variant = dict_sweep(variant, vals=['', 'null', 'N/A', None, [], {}])   # remove empty fields

            yield {  # commit an entry by yielding
                "_id": _id,
                source_name: variant
            }
