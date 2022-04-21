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

from SPDI import SPDI

logger = logging.getLogger(__name__)

tests = [
    # testing ConvertIndelToSPDI
    (3301,  "AGA",   "AG",   3303,  "AG",       "G",       "",    "ConvertIndel test 1"),
    (3481,  "A",     "AG",   3482,  "G",        "GG",      "",    "ConvertIndel test 2"),
    (1001,  "",      "A",    1001,  "AAAA",     "AAAAA",   "",    "ConvertIndel test 3"),
    (1001,  "-",     "A",    1001,  "AAAA",     "AAAAA",   "",    "ConvertIndel test 4"),
    (11074, "CT",    "C",    11075, "TTTTTTTT", "TTTTTTT", "",    "ConvertIndel test 5"),
    (11082, "T",     "",     11075, "TTTTTTTT", "TTTTTTT", "",    "ConvertIndel test 6"),
    (14747, "A",     "",     14747, "AA",       "A",       "",    "ConvertIndel test 7"),
    (14746, "GAATT", "GATT", 14747, "AA",       "A",       "",    "ConvertIndel test 8"),
    (14746, "GAA",   "GA",   14747, "AA",       "A",       "",    "ConvertIndel test 9"),
    (14843, "AACT",  "A",    14844, "ACTACTA",  "ACTA",    "",    "ConvertIndel test 10"),
    (21881, "GGTAC", "GG",   21883, "TACTACT",  "TACT",    "",    "ConvertIndel test 11"),
    (9430,  "C",     "CGTA", 9431,  "GTA",      "GTAGTA",  "",    "ConvertIndel test 12"),
    (12297, "",      "A",    12297, "AA",       "AAA",     "",    "ConvertIndel test 13"),
    (11696, "TATGA", "TA",   11697, "ATGAT",    "AT",      "",    "ConvertIndel test 17"),
    (11696, "TATGAT","TA",   11698, "TGATT",    "T",       "",    "ConvertIndel test 18"),
    (18668, "GT",    "G",    18669, "TG",       "G",       "",    "ConvertIndel test 19"),
    (1,     "ATT",   "CATT", 1,     "",         "C",       "",    "ConvertIndel test 20"),
    (12297, "-",     "A",    12297, "AA",       "AAA",     "",    "ConvertIndel test 21"),
    (18669, "T",     "-",    18669, "TG",       "G",       "",    "ConvertIndel test 22"),
    (13200, "CGG",   "TGG",  13200, "C",        "T",       "SNP", "ConvertIndel test 14"),
    (13200, "CGG",   "CAG",  13201, "G",        "A",       "SNP", "ConvertIndel test 15"),
    (13201, "G",     "A",    13201, "G",        "A",       "SNP", "ConvertIndel test 21"),
    (13201, "G",     "G",    13201, "G",        "G",       "SNP", "ConvertIndel test 22"),

    # testing ConvertIndelToSPDI: unconvertable inputs
    (1001, "G",   "B",   0, "",     "",   "invalid alt 'B'", "ConvertIndel test 23"),
    (1001, "9",   "A",   0, "",     "",   "invalid ref '9'", "ConvertIndel test 24"),
    (1001, "G",   "AaT", 0, "",     "",   "invalid alt 'AaT'", "ConvertIndel test 25"),
    (1001, "T",   "A",   0, "",     "",   "", "ConvertIndel test 26"),
    (1001, "",    "",    0, "",     "",   "empty mutation", "ConvertIndel test 27"),
    (1001, "GAA", "CTT", 0, "GAA", "CTT", "invalid mutation", "ConvertIndel test 28"),
    (1001, "GAA", "CA",  0, "GA",  "C",   "invalid mutation", "ConvertIndel test 29"),
    (1001, "GA",  "CAT", 0, "GA",  "CAT", "invalid mutation", "ConvertIndel test 30"),
    (99999,"G",   "A",   0, "",    "",    "invalid pos 99999", "ConvertIndel test 31"),
]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--r',
                        help='path to reference sequence fasta file',
                        default='/panfs/traces01.be-md.ncbi.nlm.nih.gov/sra_review/scratch/kocherginai/static/reference/NC_045512.2.fa')

    parser.add_argument('--log-level', help='log level', choices=list(_nameToLevel.keys()),
                        default='DEBUG')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s.%(msecs)03d %(module)s..%(funcName)s [%(levelname)s] %(message)s',
                        datefmt='%m-%d-%Y %H:%M:%S')

    with SPDI(args.r) as spdi:
        print(f"{spdi.genome_length} == 29903")
        assert spdi.genome_length == 29903

        for t in tests:
            oldPos, oldRef, oldAlt, expPos, expRef, expAlt, expErr, testName = t
            newPos, newRef, newAlt, error = spdi.convert_indel_to_SPDI(oldPos, oldRef, oldAlt)

            print(f"{testName} {newPos} == {expPos}, {newRef} == {expRef}, {newAlt} == {expAlt}")

            assert newPos == expPos
            assert newRef == expRef
            assert newAlt == expAlt
            if expErr and error != expErr:
                print(f"error mismatch: {error} != {expErr}")


if __name__ == '__main__':
    main()
