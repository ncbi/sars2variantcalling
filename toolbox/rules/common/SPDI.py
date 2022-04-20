import os
import subprocess
import argparse
import logging.handlers
from logging import _nameToLevel
from typing import Dict, Tuple
from lxml.builder import E
from lxml import etree
from lxml.etree import Element
from xml.sax.saxutils import escape
import itertools
import re

"""
SPDI normalization of indels  
    - direct translation from Jimmy Jin's perl scripts 
      ConvertIndelsToSPDI.pl and SPDI.pm located in scf/toolbox/Scripts/ 
"""


if __name__ == '__main__':
    logger = logging.getLogger(__name__)
else:
    logger = None


class SPDI:
    alleles = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self):
        self.ref = []
        self.seq = None
        self.seq_length = 0

    def load(self, path):
        with open(path, 'rt') as _fa:
            for line in _fa:
                _l = line.rstrip()
                if not _l.startswith('>'):
                    self.ref.append(_l)
        self.seq = ''.join(self.ref)
        self.seq_length = len(self.seq)

        logger and logger.info(self.seq)

    def check_alleles_in_sequence(self, nuc_seq_ref: str):
        for ch in nuc_seq_ref:
            if ch not in self.alleles.keys():
                return 0
        return 1

    @staticmethod
    def trim_left(o_ref, o_alt):
        if not o_ref or not o_alt:
            return o_ref, o_alt, 0
        pos = 0
        while pos < len(o_ref) and pos < len(o_alt):
            if o_ref[pos] != o_alt[pos]:
                break
            pos += 1

        logger and logger.debug(f"left pos {pos}, ref {o_ref[pos:]}, alt {o_alt[pos:]}")
        return o_ref[pos:], o_alt[pos:], pos

    @staticmethod
    def trim_right(o_ref, o_alt):
        pos = 0
        while pos < len(o_ref) and pos < len(o_alt):
            if o_ref[-pos-1] != o_alt[-pos-1]:
                break
            pos += 1

        logger and logger.debug(f"right pos {pos}, ref {o_ref[:len(o_ref)-pos]}, alt {o_alt[:len(o_alt)-pos]}")
        return o_ref[:len(o_ref)-pos], o_alt[:len(o_alt)-pos], pos

    def blossom(self, o_pos, bud):
        bud_length = len(bud)

        logger and logger.debug(f"Blossoming bud {bud} at position {o_pos}")

        new_ref, new_alt = bud, ''
        pos = o_pos + bud_length

        while True:
            if pos > self.seq_length:
                break
            for i in range(0, bud_length):
                seq_nt = self.seq[pos - 1]
                bud_nt = bud[i]

                logger and logger.debug(f"pos {pos}, seq_nt {seq_nt}, bud_nt {bud_nt}")

                pos += 1


class VCF:
    def __init__(self):
        self.all_lines = []
        self.header = []
        self.format = None
        self.vcf = None

    def load(self, path):
        with open(path, 'rt') as _vcf:
            for line in _vcf:

                _l = line.rstrip()
                if _l.startswith('#'):
                    self.header.append(_l)
                else:
                    self.all_lines.append(_l.split('\t'))


input_description = ('ref:path', 'vcf:path')
output_description = ('vcf:path', 'summary:path')


def run(_input: Dict, _output: Dict, _config: Dict, _rule: Dict) -> int:
    """
    :param _input:
    :param _output:
    :param _config:
    :param _rule:
    :return:
    """

    with open(_input["vcf"], 'rt') as _i_vcf, open(_output["vcf"], 'wt') as _o_vcf:
        for _line in _i_vcf:
            print(_line)
    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--r', help='path to reference sequence fasta file', required=True)
    parser.add_argument('--i', help='path to input VCF file', required=True)
    parser.add_argument('--o', help='path to output VCF file', required=True)
    parser.add_argument('--s', help='path to output summary file', required=True)

    parser.add_argument('--log-level', help='log level', choices=list(_nameToLevel.keys()),
                        default='INFO')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s.%(msecs)03d %(module)s..%(funcName)s [%(levelname)s] %(message)s',
                        datefmt='%m-%d-%Y %H:%M:%S')


if __name__ == '__main__':
    main()
