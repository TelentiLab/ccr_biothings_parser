import os
import logging
import datetime
import time

FILE_NOT_FOUND_ERROR = 'Cannot find input file: {}'  # error message constant

# change following parameters accordingly
source_name = 'ccr'  # source name that appears in the api response
_files = [
    ('ccrs.autosomes.v2.20180420.bed', 8188410),  # autosome (file, lines)
    ('ccrs.xchrom.v2.20180420.bed', 171987),  # xchrom (file, lines)
]
delimiter = '\t'  # the delimiter that separates each field

# configure logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', level=logging.INFO)
logger = logging.getLogger('CCRs_parser')


def version(self):
    """
    This info will be shown in the API response under `version` field
    :return:
    """
    return 'v2.0.1'


def load_data(data_folder: str):
    """
    Load data from a specified file path. Parse each line into a dictionary according to the schema
    given by `data_schema`. Then process each dict by normalizing data format, remove null fields (optional).
    Append each dict into final result using its id.

    :param data_folder: the path(folder) where the data file is stored
    :return: a generator that yields data.
    """
    for file_name, file_lines in _files:
        input_file = os.path.join(data_folder, file_name)
        # raise an error if file not found
        if not os.path.exists(input_file):
            logger.error(FILE_NOT_FOUND_ERROR.format(input_file))
            raise FileExistsError(FILE_NOT_FOUND_ERROR.format(input_file))

        with open(input_file, 'r') as file:
            start_time = time.time()
            logger.info(f'start reading file: {file_name}')
            count = 0
            skipped = []
            for line in file:   # read and parse each line into a dict
                count += 1
                ratio = count / file_lines
                time_left = datetime.timedelta(seconds=(time.time() - start_time) * (1 - ratio) / ratio)
                # format to use 2 decimals for progress
                logger.info(
                    f'reading line {count} ({(count / file_lines * 100):.2f}%), estimated time left: {time_left}')

                # schema: (chrom, start, end, ccr_pct, gene, ranges, varflag, syn_density,
                #          cpg, cov_score, resid, resid_pctile, unique_key)
                try:
                    (chrom, start, end, ccr_pct, gene, ranges, varflag, syn_density, cpg,
                        cov_score, resid, resid_pctile, unique_key) = line.strip().split(delimiter)
                except ValueError:
                    logger.error(f'failed to unpack line {count}: {line}')
                    skipped.append(line)
                    continue    # skip error line
                _id = f'chr{chrom}:g.{start}_{end}'
                # enforce data type
                try:
                    variant = {
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'scores': [
                            {
                                'ccr_pct': float(ccr_pct),
                                'gene': gene,
                                'ranges': ranges.split(','),
                                'varflag': [each is 'VARTRUE' for each in varflag.split(',')],
                                'syn_density': float(syn_density),
                                'cpg': float(cpg),
                                'cov_score': float(cov_score),
                                'resid': float(resid),
                                'resid_pctile': float(resid_pctile),
                                'unique_key': int(unique_key),
                            }
                        ]
                    }
                except ValueError as e:
                    logger.error(f'failed to cast type for line {count}: {e}')
                    skipped.append(line)
                    continue  # skip error line

                yield {  # commit an entry by yielding
                    "_id": _id,
                    source_name: variant
                }
            logger.info(f'parse completed, {len(skipped)}/{count} lines skipped.')
            for x in skipped:
                logger.info(f'skipped line: {x}')
