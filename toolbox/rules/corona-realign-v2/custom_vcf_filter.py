import argparse
import logging.handlers
import sys
import os
from logging import _nameToLevel
from typing import Dict
import json
import numpy as np
import collections

"""
custom VCF filter, which could not be or hard to implement using GATK VariantFiltration or Selection. 
see SRAB-308
    && $ad_alt / $DP_unfilt >= 0.15 && $dep >= 10
        ignoring those variants with AF < 0.15;  
        "ad_alt" for ALT allele counts from FORMAT field, 
        "DP_unfilt" for DP from INFO field;  
        "dep" for the coverage of this variant POS from bam_genomecoverage
    && $dep_min_next10bps >= 10
        ignoring this variant if the minimum coverage of NEXT 10bps after this POS is 10, 
        for getting rid of some bad calls from read ends
    && $dp_ori / $DP_unfilt >= 0.5 && $dep / $dp_ori >= 0.5
        trying to filter out bad calls from regions with too much differences between DP from INFO field (DP_unfilt), 
        "dp_ori" (DP from FORMAT), and "dep" (DP from bam_genomecoverage).
"""

if __name__ == '__main__':
    logger = logging.getLogger(__name__)
else:
    logger = None


class VCF:
    row = collections.namedtuple('row', 'pos ref alt qual info formats')
    val = collections.namedtuple('val', 'ad_alt dp_unfilt dp_ori del_sz pos_d dep dep_min_next10bps')

    def __init__(self, vcf_file, depth_file):
        self.meta = collections.defaultdict(list)
        self.locations = []
        self.header = []
        self.format = None
        self.columns = None
        self.vcf = None
        self.depth = None
        self.vcf_file = vcf_file
        self.depth_file = depth_file

    def __exit__(self, *args): ...

    def __enter__(self):
        self.load_depth(self.depth_file)
        self.load_vcf(self.vcf_file)
        return self

    def parse(self):
        def _map(pair):
            if len(pair) == 1:
                pair.append(True)
            return pair
        for loc in self.locations:
            format_def = loc[8].split(':')
            yield self.row._make([int(loc.pos),
                                  loc.ref,
                                  loc.alt,
                                  loc.qual,
                                  dict([_map(a.split('=')) for a in loc.info.split(';')]),
                                  [dict(zip(format_def, f.split(':'))) for f in loc[9:]]
                                 ]), loc

    def get_min_depth_4_window(self, pos, window=10):
        return self.depth[pos - 1:pos + window - 1].min(initial=10)

    def get_depth_window(self, pos, window=10):
        return self.depth[pos - 1:pos + window - 1].tolist()

    def get_depth(self, pos):
        return self.depth[pos - 1]

    def load_depth(self, depth):
        self.depth = np.loadtxt(depth, delimiter='\t', usecols=[2])

    def load_vcf(self, vcf):
        with open(vcf) as fp:
            for e, l in enumerate(fp):

                l = l.rstrip()
                if self.vcf is None:
                    if l.startswith('##fileformat='):
                        self.vcf = l.split('=')[1]
                        self.header.append(l)
                        continue
                    else:
                        raise Exception('It does not look like VCF file')
                elif l.startswith('##'):
                    self.header.append(l)
                    pair = l[2:].split('=', 1)
                    if pair[0] in ['INFO', 'FORMAT', 'FILTER']:
                        self.meta[pair[0]].append(pair[1][4:-1].split(',')[0])
                    else:
                        self.meta[pair[0]]=pair[1]
                    continue
                elif l.startswith('#'):
                    self.header.append(l)
                    cols = ' '.join([c.lower() if c.isalnum() else os.path.basename(c).split('.')[0].lower() for c in l[1:].split('\t')])
                    self.columns = collections.namedtuple('vcf', cols)
                    continue
                else:
                    self.locations.append(self.columns._make(l.rstrip().split('\t')))


input_description = ('vcf:path', 'coverage:path')
output_description = ('vcf:path', )


def run(_input: Dict, _output: Dict, _config: Dict, _rule: Dict) -> int:
    """
    :param _input:
    :param _output:
    :param _config:
    :param _rule:
    :return:
    """

    snps = 0
    filtered = 0
    dropouts = collections.defaultdict(int)
    with VCF(_input["vcf"], _input["coverage"]) as vcf:
        with open(_output["vcf"], 'wt') as o_vcf:
            for h in vcf.header:
                o_vcf.write(h)
                o_vcf.write('\n')
            for v, r in vcf.parse():
                del_sz = len(v.ref) - len(v.alt)
                pos_d = v.pos + del_sz + 1 if del_sz > 0 else v.pos + len(v.ref) + 1

                calc = vcf.val(ad_alt=int(v.formats[0]['AD'].split(',')[1]),  # for ALT allele counts from FORMAT field
                               dp_unfilt=int(v.info['DP']),                   # for DP from INFO field
                               dp_ori=int(v.formats[0]['DP']),                # DP from FORMAT
                               del_sz=del_sz, pos_d=pos_d,
                               dep=vcf.get_depth(v.pos),
                               dep_min_next10bps=vcf.get_min_depth_4_window(pos_d))

                if v.alt != '*' and r.filter == 'PASS' \
                        and calc.dep_min_next10bps >= 10 \
                        and calc.dp_ori / calc.dp_unfilt >= 0.5 \
                        and calc.dep / calc.dp_ori >= 0.5 \
                        and calc.ad_alt / calc.dp_unfilt >= 0.15 \
                        and calc.dep >= 10:
                    o_vcf.write('\t'.join(r))
                    o_vcf.write('\n')
                    snps += 1
                else:
                    filters = dict()
                    filtered += 1
                    if r.filter != 'PASS':
                        filters[r.filter] = ''
                    if v.alt == '*':
                        filters['altStar'] = '*'
                    if calc.dep_min_next10bps < 10:
                        filters['lowCovTail10'] = calc.dep_min_next10bps
                    if calc.dp_ori / calc.dp_unfilt < 0.5:
                        filters['lowRatioInfoDP2fmtDP.5'] = calc.dp_ori / calc.dp_unfilt
                    if calc.dep / calc.dp_ori < 0.5:
                        filters['lowRatioCov2infoDP.5'] = calc.dep / calc.dp_ori
                    if calc.ad_alt / calc.dp_unfilt < 0.15:
                        filters['lowRatioCov2infoDP.15'] = calc.ad_alt / calc.dp_unfilt
                    if calc.dep < 10:
                        filters['lowCov10'] = calc.dep
                    dropouts[';'.join(sorted(filters.keys()))] += 1

                    logger and logger.debug(f"pos = {v.pos}, filter = {filters}, variant = {r}, format = {dict(zip(r.format.split(':'),r[-1].split(':')))}, calc = {calc}")

    print(json.dumps({'passed': snps,
                      'filtered': filtered,
                      'filters': dropouts}), end='')

    if not snps:
        sys.stderr.write('no-SNPs')
        return 1

    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--c',
                        help='path to depth file, built using bedtools genomecov -ibam $bam -g $ref -d > $depth',
                        required=True)
    parser.add_argument('--i',
                        help='path to input VCF file',
                        required=True)
    parser.add_argument('--o',
                        help='path to output VCF file',
                        required=True)

    parser.add_argument('--log-level', help='log level', choices=list(_nameToLevel.keys()),
                        default='INFO')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s.%(msecs)03d %(module)s..%(funcName)s [%(levelname)s] %(message)s',
                        datefmt='%m-%d-%Y %H:%M:%S')

    run(_input={"vcf": args.i, "coverage": args.c},
        _output={"vcf": args.o},
        _config={}, _rule={})


if __name__ == '__main__':
    # for testing purposes
    main()
