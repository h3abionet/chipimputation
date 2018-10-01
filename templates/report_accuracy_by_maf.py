#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inSNP_acc", help="")
parser.add_argument("--report_acc", help="")
args = parser.parse_args()


def acc_by_maf(inSNP_acc, outSNP_acc):
    """
    :return:
    """
    datas = {}

    outSNP_acc_out = open(outSNP_acc, 'w')
    outSNP_acc_out.writelines('\\t'.join(['CHRM', 'MAF<1%', 'MAF 1-5%', 'MAF>5%']) + '\\n')
    info_datas = open(inSNP_acc).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['rare'] = []  # < 1%
            datas[dataset]['moderate'] = []  # 1-5%
            datas[dataset]['common'] = []  # >5%
            datas[dataset]['total'] = 0
        if "SNP" in line and "Rsq" in line:
            idx_exp_freq_a1 = data.index('MAF')
            # idx_conc = data.index("concord_type0")
            idx_conc = data.index("LooRsq")
        else:
            maf = float(data[idx_exp_freq_a1])
            acc = float(data[idx_conc])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1 - maf
            if maf < 0.01:
                datas[dataset]['rare'].append(acc)
            elif maf > 0.05:
                datas[dataset]['common'].append(acc)
            elif maf <= 0.05 and maf >= 0.01:
                datas[dataset]['moderate'].append(acc)
    try:
        datasets = [str(it) for it in sorted([int(it) for it in datas])]
    except:
        datasets = sorted(datas)
    for dataset in datasets:
        tot = datas[dataset]['total']

        if len(datas[dataset]['common']) == 0:
            maf_5 = 0.0
        else:
            maf_5 = sum(datas[dataset]['common']) / float(len(datas[dataset]['common']))
        if len(datas[dataset]['rare']) == 0:
            maf_1 = 0.0
        else:
            maf_1 = sum(datas[dataset]['rare']) / float(len(datas[dataset]['rare']))
        if len(datas[dataset]['moderate']) == 0:
            maf_1_5 = 0.0
        else:
            maf_1_5 = sum(datas[dataset]['moderate']) / float(len(datas[dataset]['moderate']))

        outSNP_acc_out.write("{}\\t{:3.3f}\\t{:3.3f}\\t{:3.3f}\\n".format(dataset, maf_1, maf_1_5, maf_5))

    outSNP_acc_out.close()


args.inSNP_acc = "${inSNP_acc}"
args.outSNP_acc = "${outSNP_acc}"
if args.inSNP_acc and args.outSNP_acc:
    acc_by_maf(args.inSNP_acc, args.outSNP_acc)
