#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--outWell_imputed", help="")
parser.add_argument("--inWell_imputed", help="")

args = parser.parse_args()


def well_imputed_by_maf(inWell_imputed, outWell_imputed, group=''):
    """

    :return:
    """
    datas = {}
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_out_1 = open(outWell_imputed + "_summary.tsv", 'w')
    outWell_imputed_out_1.writelines('\\t'.join([group, 'MAF<1%', 'MAF 1-5%', 'MAF>5%', 'TOTAL']) + '\\n')
    outWell_imputed_out.writelines('\\t'.join([group, 'MAF<1%', 'MAF 1-5%', 'MAF>5%']) + '\\n')
    info_datas = open(inWell_imputed).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['rare'] = []
            datas[dataset]['moderate'] = []
            datas[dataset]['common'] = []
            datas[dataset]['total'] = 0
        if "SNP" in line and "Rsq" in line:
            idx_exp_freq_a1 = data.index('MAF')
        else:
            maf = float(data[idx_exp_freq_a1])
            if maf >= 0.5:
                maf = 1 - maf
            if maf < 0.01:
                datas[dataset]['rare'].append(maf)
                datas[dataset]['total'] += 1
            elif maf > 0.05:
                datas[dataset]['common'].append(maf)
                datas[dataset]['total'] += 1
            elif maf <= 0.05 and maf >= 0.01:
                datas[dataset]['moderate'].append(maf)
                datas[dataset]['total'] += 1
    try:
        datasets = [str(it) for it in sorted([int(it) for it in datas])]
    except:
        datasets = sorted(datas)
    for dataset in datasets:
        tot = datas[dataset]['total']
        rare = len(datas[dataset]['rare'])
        moderate = len(datas[dataset]['moderate'])
        common = len(datas[dataset]['common'])
        if rare >= 1000000:
            rare_ = str(format(rare / 1000000., '0,.1f')) + 'M ('
        else:
            rare_ = str(format(rare, '0,')) + ' ('
        if moderate >= 1000000:
            moderate_ = str(format(moderate / 1000000., '0,.1f')) + 'M ('
        else:
            moderate_ = str(format(moderate, '0,')) + ' ('
        if common >= 1000000:
            common_ = str(format(common / 1000000., '0,.1f')) + 'M ('
        else:
            common_ = str(format(common, '0,')) + ' ('
        if tot == 0:
            outWell_imputed_out.write("dataset {} is empty (tot=0)".format(dataset))
            outWell_imputed_out_1.write("dataset {} is empty (tot=0)".format(dataset))
        else:
            outWell_imputed_out.writelines('\\t'.join([dataset, \
                                                       str(len(datas[dataset]['rare'])), \
                                                       str(len(datas[dataset]['moderate'])), \
                                                       str(len(datas[dataset]['common']))]) + '\\n')
            outWell_imputed_out_1.writelines('\\t'.join([dataset, \
                                                         rare_ + str(rare * 100 / tot) + '%)', \
                                                         moderate_ + str(moderate * 100 / tot) + '%)', \
                                                         common_ + str(common * 100 / tot) + '%)', \
                                                         str(format(tot / 1000000., '0,.1f')) + 'M']) + '\\n')
    outWell_imputed_out_1.close()
    outWell_imputed_out.close()


args.inWell_imputed = "${inWell_imputed}"
args.outWell_imputed = "${outWell_imputed}"
args.group = "${group}"
if args.inWell_imputed and args.outWell_imputed:
    well_imputed_by_maf(args.inWell_imputed, args.outWell_imputed, args.group)
