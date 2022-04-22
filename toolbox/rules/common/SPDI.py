import argparse
import logging.handlers
from logging import _nameToLevel
from typing import Dict
import json

"""
SPDI normalization of indels  
    - direct translation from Jimmy Jin's perl scripts 
        ConvertIndelsToSPDI.pl and SPDI.pm located in scf/toolbox/Scripts/
    - paper?: "SPDI: data model for variants and applications at NCBI" 
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7523648/
        https://www.researchgate.net/publication/330787665_SPDI_Data_Model_for_Variants_and_Applications_at_NCBI   
"""


if __name__ == '__main__':
    logger = logging.getLogger(__name__)
else:
    logger = None


class SPDI:
    alleles = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, path: str):
        self.ref = []
        self.ref_nucs = None
        self.genome_length = 0
        self.path = path

    def __enter__(self):
        with open(self.path, 'rt') as _fa:
            for line in _fa:
                _l = line.rstrip()
                if not _l.startswith('>'):
                    self.ref.append(_l)
        self.ref_nucs = ''.join(self.ref)
        self.genome_length = len(self.ref_nucs)

        logger and logger.debug(self.ref_nucs)

        return self

    def __exit__(self, *__args): ...

    def check_alleles_in_sequence(self, nuc_seq_ref: str):
        for ch in nuc_seq_ref:
            if ch not in self.alleles.keys():
                return False
        return True

    @staticmethod
    def trim_left(o_ref: str, o_alt: str):
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
    def trim_right(o_ref: str, o_alt: str):
        pos = 0
        while pos < len(o_ref) and pos < len(o_alt):
            if o_ref[-pos-1] != o_alt[-pos-1]:
                break
            pos += 1

        logger and logger.debug(f"right pos {pos}, ref {o_ref[:len(o_ref)-pos]}, alt {o_alt[:len(o_alt)-pos]}")
        return o_ref[:len(o_ref)-pos], o_alt[:len(o_alt)-pos], pos

    def blossom(self, o_pos: int, bud: str):
        bud_length = len(bud)

        logger and logger.debug(f"Blossoming bud {bud} at position {o_pos} to the right")

        new_ref, new_alt = bud, ''
        pos = o_pos + bud_length

        done = False
        while not done:
            for bud_nt in bud:
                if pos > self.genome_length:
                    done = True
                    break
                seq_nt = self.ref_nucs[pos - 1]

                logger and logger.debug(f"pos {pos}, seq_nt {seq_nt}, bud_nt {bud_nt}")
                if seq_nt == bud_nt:
                    new_ref += seq_nt
                    new_alt += seq_nt
                else:
                    done = True
                    break
                pos += 1

        logger and logger.debug(f"pos {pos}, new ref {new_ref}, new alt {new_alt}")

        logger and logger.debug(f"Blossoming bud {bud} at position {o_pos} to the left")

        pos = o_pos - 1

        done = False
        while not done:
            for i in range(0, bud_length):
                if pos < 1:
                    done = True
                    break
                seq_nt = self.ref_nucs[pos - 1]
                bud_nt = bud[-i-1]

                logger and logger.debug(f"pos {pos}, seq_nt {seq_nt}, bud_nt {bud_nt}")
                if seq_nt == bud_nt:
                    new_ref = seq_nt + new_ref
                    new_alt = seq_nt + new_alt
                else:
                    done = True
                    break

                pos -= 1

        logger and logger.debug(f"pos {pos}, new ref {new_ref}, new alt {new_alt}")

        return new_ref, new_alt, pos + 1

    def convert_indel_to_SPDI(self, o_pos: int, o_ref: str, o_alt: str):
        if o_pos < 1 or o_pos > self.genome_length:
            return 0, "", "", f"invalid pos {o_pos}"

        pos = o_pos
        ref = o_ref if o_ref != '-' else ''
        alt = o_alt if o_alt != '-' else ''

        logger and logger.debug(f"Checking mutation pos {pos}, ref {ref}, alt {alt}")

        # Make sure input ref and alt are string of A,T,G,C
        if not self.check_alleles_in_sequence(ref):
            return 0, "", "", f"invalid ref '{o_ref}'"
        if not self.check_alleles_in_sequence(alt):
            return 0, "", "", f"invalid alt '{o_alt}'"

        # Make sure ref is correct for the position
        seq = self.ref_nucs[pos - 1:pos - 1 + len(ref)]
        if seq != ref:
            return 0, "", "", f"reference '{o_ref}' does not match seq '{seq}'"

        if len(ref) == 1 and len(alt) == 1:
            return pos, o_ref, o_alt, "SNP"

        new_ref, new_alt, _ = self.trim_right(ref, alt)
        new_ref, new_alt, num_left_trims = self.trim_left(new_ref, new_alt)

        if not new_ref and new_alt:
            it_is_insert = True
            bud = new_alt
        elif new_ref and not new_alt:
            it_is_insert = False
            bud = new_ref
        elif new_ref and new_alt:
            if len(new_ref) == 1 and len(new_alt) == 1:
                logger and logger.debug(f"WARNING: mutation ref {o_ref}, alt {o_alt} is a SNP")
                return o_pos + num_left_trims, new_ref, new_alt, "SNP"
            else:
                logger and logger.debug(f"WARNING: mutation ref {o_ref}, alt {o_alt} is not valid")
                return 0, new_ref, new_alt, "invalid mutation"
        else:  # if not new_ref and not new_alt:
            logger and logger.debug(f"WARNING: mutation ref {o_ref}, alt {o_alt} is not valid: ref == alt")
            return 0, '', '', "empty mutation"

        mut_ref, mut_alt, mut_pos = self.blossom(o_pos + num_left_trims, bud)

        out_ref, out_alt = (mut_alt, mut_ref) if it_is_insert else (mut_ref, mut_alt)

        if not out_ref and mut_pos > 1:
            mut_pos -= 1
            prev_nuc = self.ref_nucs[mut_pos - 1]
            out_ref = prev_nuc
            out_alt = prev_nuc + out_alt

        if not out_alt:
            ref_len = len(out_ref)
            if mut_pos and mut_pos + ref_len < self.genome_length:
                next_nuc = self.ref_nucs[mut_pos + ref_len - 1]
                out_ref += next_nuc
                out_alt = next_nuc

        logger and logger.debug("Final indel call:")
        logger and logger.debug(f"ref {ref} alt {alt} pos {o_pos}")
        logger and logger.debug("to")
        logger and logger.debug(f"ref {out_ref} alt {out_alt} pos {mut_pos}")

        return mut_pos, out_ref, out_alt, ''


class VCF:
    I_POS = 1
    I_REF = 3
    I_ALT = 4

    def __init__(self, path):
        self.variations = []
        self.header = []
        self.path = path

    def __exit__(self, *args): ...

    def __enter__(self):
        with open(self.path, 'rt') as _vcf:
            for line in _vcf:

                _l = line.rstrip()
                if _l.startswith('#'):
                    self.header.append(_l)
                else:
                    self.variations.append(_l.split('\t'))
        return self


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

    with VCF(_input["vcf"]) as vcf, SPDI(_input["ref"]) as spdi:
        with open(_output["vcf"], 'wt') as o_vcf, open(_output["summary"], 'wt') as o_sum:
            for h in vcf.header:
                o_vcf.write(h)
                o_vcf.write('\n')
            unique_positions = set()
            converted = 0
            failed = 0
            left_shift = 0
            for v in vcf.variations:
                pos = int(v[VCF.I_POS])
                ref = v[VCF.I_REF].replace('-', '')
                alt = v[VCF.I_ALT].replace('-', '')

                unique_positions.add(pos)

                if len(ref) != 1 or len(alt) != 1:
                    n_pos, n_ref, n_alt, err = spdi.convert_indel_to_SPDI(pos, ref, alt)
                    if n_pos:
                        v[VCF.I_POS] = str(n_pos)
                        v[VCF.I_REF] = n_ref
                        v[VCF.I_ALT] = n_alt
                        message = f", message '{err}'" if err else ''
                        shift = ''
                        if pos > n_pos:
                            left_shift += 1
                            shift = '!'
                        o_sum.write(f"converted: {pos} -> {n_pos}{shift}, {ref} -> {n_ref}, {alt} -> {n_alt}{message}\n")
                        converted += 1
                    else:
                        o_sum.write(f"failed to convert: {pos}, {ref}, {alt}, message {err}\n")
                        failed += 1

                o_vcf.write('\t'.join(v))
                o_vcf.write('\n')

            print(json.dumps({"total-variants": len(vcf.variations),
                              "unique-positions": len(unique_positions),
                              "converted": converted,
                              "failed": failed,
                              "left-shifts": left_shift
                              }), end='')

    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--r',
                        help='path to reference sequence fasta file',
                        required=True)
    parser.add_argument('--i',
                        help='path to input VCF file',
                        required=True)
    parser.add_argument('--o',
                        help='path to output VCF file',
                        required=True)
    parser.add_argument('--s',
                        help='path to output summary file',
                        required=True)

    parser.add_argument('--log-level', help='log level', choices=list(_nameToLevel.keys()),
                        default='INFO')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s.%(msecs)03d %(module)s..%(funcName)s [%(levelname)s] %(message)s',
                        datefmt='%m-%d-%Y %H:%M:%S')

    run(_input={"ref": args.r, "vcf": args.i},
        _output={"vcf": args.o, "summary": args.s},
        _config={}, _rule={})


if __name__ == '__main__':
    # for testing purposes
    main()
