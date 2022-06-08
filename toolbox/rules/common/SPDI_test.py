import argparse
import logging.handlers
from logging import _nameToLevel

from SPDI import SPDI

logger = logging.getLogger(__name__)

_tn = iter(range(1, 100))

tests = [
    # testing ConvertIndelToSPDI
   #(3301,  "AGA",   "AG",   3303,  "AG",       "G",       "",    "ConvertIndel test 1"),
    (3301,  "AGA",   "AG",   3302,  "GA",       "G",       "",    f"ConvertIndel test {next(_tn)}"),
    (3481,  "A",     "AG",   3482,  "G",        "GG",      "",    f"ConvertIndel test {next(_tn)}"),
    (1001,  "",      "A",    1001,  "AAAA",     "AAAAA",   "",    f"ConvertIndel test {next(_tn)}"),
    (1001,  "-",     "A",    1001,  "AAAA",     "AAAAA",   "",    f"ConvertIndel test {next(_tn)}"),
    (11074, "CT",    "C",    11075, "TTTTTTTT", "TTTTTTT", "",    f"ConvertIndel test {next(_tn)}"),
    (11082, "T",     "",     11075, "TTTTTTTT", "TTTTTTT", "",    f"ConvertIndel test {next(_tn)}"),
    (14747, "A",     "",     14747, "AA",       "A",       "",    f"ConvertIndel test {next(_tn)}"),
    (14746, "GAATT", "GATT", 14747, "AA",       "A",       "",    f"ConvertIndel test {next(_tn)}"),
    (14746, "GAA",   "GA",   14747, "AA",       "A",       "",    f"ConvertIndel test {next(_tn)}"),
    (14843, "AACT",  "A",    14844, "ACTACTA",  "ACTA",    "",    f"ConvertIndel test {next(_tn)}"),
    (21881, "GGTAC", "GG",   21883, "TACTACT",  "TACT",    "",    f"ConvertIndel test {next(_tn)}"),
    (9430,  "C",     "CGTA", 9431,  "GTA",      "GTAGTA",  "",    f"ConvertIndel test {next(_tn)}"),
    (12297, "",      "A",    12297, "AA",       "AAA",     "",    f"ConvertIndel test {next(_tn)}"),
    (11696, "TATGA", "TA",   11697, "ATGAT",    "AT",      "",    f"ConvertIndel test {next(_tn)}"),
    (11696, "TATGAT","TA",   11698, "TGATT",    "T",       "",    f"ConvertIndel test {next(_tn)}"),
   #(18668, "GT",    "G",    18669, "TG",       "G",       "",    "ConvertIndel test 19"),
    (18668, "GT",    "G",    18668, "GT",       "G",       "",    "ConvertIndel test {next(_tn)}"),
    (1,     "ATT",   "CATT", 1,     "",         "C",       "",    f"ConvertIndel test {next(_tn)}"),
    (12297, "-",     "A",    12297, "AA",       "AAA",     "",    f"ConvertIndel test {next(_tn)}"),
   #(18669, "T",     "-",    18669, "TG",       "G",       "",    "ConvertIndel test 22"),
    (18669, "T",     "-",    18668, "GT",       "G",       "",    f"ConvertIndel test {next(_tn)}"),
    (13200, "CGG",   "TGG",  13200, "C",        "T",       "SNP", f"ConvertIndel test {next(_tn)}"),
    (13200, "CGG",   "CAG",  13201, "G",        "A",       "SNP", f"ConvertIndel test {next(_tn)}"),
    (13201, "G",     "A",    13201, "G",        "A",       "SNP", f"ConvertIndel test {next(_tn)}"),
    (13201, "G",     "G",    13201, "G",        "G",       "SNP", f"ConvertIndel test {next(_tn)}"),

    # Two cases from https://jira.ncbi.nlm.nih.gov/browse/SRAB-396
    (11287, "GTCTGGTTTT",   "G", 11287, "GTCTGGTTTT", "G", "", f"ConvertIndel test {next(_tn)}"),
    (21764, "ATACATG",      "A", 21765, "TACATGT",    "T", "", f"ConvertIndel test {next(_tn)}"),

    # testing ConvertIndelToSPDI: unconvertable inputs
    (1001, "G",   "B",   0, "",     "",   "invalid alt 'B'",   f"ConvertIndel test {next(_tn)}"),
    (1001, "9",   "A",   0, "",     "",   "invalid ref '9'",   f"ConvertIndel test {next(_tn)}"),
    (1001, "G",   "AaT", 0, "",     "",   "invalid alt 'AaT'", f"ConvertIndel test {next(_tn)}"),
    (1001, "T",   "A",   0, "",     "",   "",                  f"ConvertIndel test {next(_tn)}"),
    (1001, "",    "",    0, "",     "",   "empty mutation",    f"ConvertIndel test {next(_tn)}"),
    (1001, "GAA", "CTT", 0, "GAA", "CTT", "invalid mutation",  f"ConvertIndel test {next(_tn)}"),
    (1001, "GAA", "CA",  0, "GA",  "C",   "invalid mutation",  f"ConvertIndel test {next(_tn)}"),
    (1001, "GA",  "CAT", 0, "GA",  "CAT", "invalid mutation",  f"ConvertIndel test {next(_tn)}"),
    (99999,"G",   "A",   0, "",    "",    "invalid pos 99999", f"ConvertIndel test {next(_tn)}"),
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

    for t in tests:
        print(t)

    with SPDI(args.r) as spdi:
        print(f"{spdi.genome_length} == 29903")
        assert spdi.genome_length == 29903

        for t in tests:
            oldPos, oldRef, oldAlt, expPos, expRef, expAlt, expErr, testName = t
            newPos, newRef, newAlt, error = spdi.convert_indel_to_SPDI(oldPos, oldRef, oldAlt)

            print(f"'{testName}' {newPos} == {expPos}, {newRef} == {expRef}, {newAlt} == {expAlt}")

            assert newPos == expPos
            assert newRef == expRef
            assert newAlt == expAlt
            if expErr and error != expErr:
                print(f"error mismatch: {error} != {expErr}")


if __name__ == '__main__':
    main()
