import argparse
import logging.handlers
from logging import _nameToLevel
from typing import Dict, List, Union
import numpy as np

"""
    bins [(0, 0), (1, 9), (10, 49), (50, 99), (100, 999), (1000, 99999999), (100000000, 9223372036854775807)]
"""


if __name__ == '__main__':
    logger = logging.getLogger(__name__)
else:
    logger = None


class Converter:
    def __init__(self, depth_file_name: str, bins_file_name: str, reference: str, bins: List):
        self.reference = reference
        self.depth_file_name = depth_file_name
        self.depth = None
        self.bins_file_name = bins_file_name
        self.bins = bins

    def __enter__(self):
        self.depth = np.loadtxt(self.depth_file_name, delimiter='\t', usecols=[2], dtype=np.int)
        logger and logger.debug(f'depth shape {self.depth.shape}')
        return self

    def __exit__(self, *__args): ...

    def process(self):
        size = self.depth.shape[0]
        comparisons = None
        for e, b in enumerate(self.bins):
            low_boundary = b[0]
            c = (self.depth >= low_boundary) * (e + 1)
            c = np.reshape(c, (size, 1))
            if comparisons is None:
                comparisons = c
            else:
                comparisons = np.concatenate((comparisons, c), axis=1)

        logger and logger.debug(f'comparisons shape {comparisons.shape}')
        selection = (np.amax(comparisons, axis=1) - 1)
        logger and logger.debug(f'selection shape {selection.shape}')
        changes = np.where(np.diff(selection, prepend=np.nan))[0]
        logger and logger.debug(f'changes shape {changes.shape}')
        with open(self.bins_file_name, 'wt') as o:
            lag = None
            for c in changes.tolist():
                if lag:
                    o.write(f'{self.reference}\t{lag[0] + 1}\t{c}\t{self.bins[lag[1]][0]}\t{self.bins[lag[1]][1]}\n')
                lag = c, selection[c]
            o.write(f'{self.reference}\t{lag[0] + 1}\t{size}\t{self.bins[lag[1]][0]}\t{self.bins[lag[1]][1]}\n')


input_description = ('depth:path', 'ref:str')
output_description = ('depth_bins:path', )


def run(_input: Dict, _output: Dict, _config: Dict, _rule: Dict) -> int:
    """
    :param _input:
    :param _output:
    :param _config:
    :param _rule:
    :return:
    """
    bins = _config["bins"] if 'bins' in _config else [(0, 0), (1, 9), (10, 49), (50, 99), (100, 999), (1000, 99999999), (100000000, 9223372036854775807)]
    with Converter(_input["depth"], _output["depth_bins"], _input["ref"], bins) as c:
        c.process()

    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--i', help='path to input depth file', required=True)
    parser.add_argument('--o', help='path to output depth bins file', required=True)
    parser.add_argument('--r', help='reference', default="NC_045512.2")

    parser.add_argument('--log-level', help='log level', choices=list(_nameToLevel.keys()),
                        default='DEBUG')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s.%(msecs)03d %(module)s..%(funcName)s [%(levelname)s] %(message)s',
                        datefmt='%m-%d-%Y %H:%M:%S')

    run(_input={"depth": args.i, "ref": args.r}, _output={"depth_bins": args.o}, _config={}, _rule={})


if __name__ == '__main__':
    # for testing purposes
    main()
