#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--infoFiles", help="")
parser.add_argument("--outWell_imputed", help="")
parser.add_argument("--outSNP_acc", help="")
parser.add_argument("--infoCutoff", help="")

args = parser.parse_args()


def filter_info(infoFiles, infoCutoff, outWell_imputed, outSNP_acc):
    """
    Return:
        well_imputed: certainy >= 1
        SNP_concordance: concord_type0 != -1
    """
    well_imputed = {}
    SNP_concordance = {}
    count = 0
    infoFiles = infoFiles.split(',')
    header = []
    outWell_imputed_out = open(outWell_imputed + ".tsv", 'w')
    outWell_imputed_snp_out = open(outWell_imputed + "_snp.tsv", 'w')
    outSNP_accuracy_out = open(outSNP_acc + ".tsv", 'w')
    for infoFile in infoFiles:
        infoFile = infoFile.strip().split('==')
        dataset = infoFile[0]
        info = infoFile[1]
        well_imputed[dataset] = []
        SNP_concordance[dataset] = []
        print info
        for line in open(info):
            data = line.strip().split()
            if "SNP" in line and "Rsq" in line:
                if len(header) == 0:
                    header = data
                    info_idx = header.index("Rsq")
                    conc_idx = header.index("EmpRsq")
                    outWell_imputed_out.writelines(' '.join([dataset] + data) + '\\n')
                    outWell_imputed_snp_out.writelines(data[1] + '\\n')
                    outSNP_accuracy_out.writelines(' '.join([dataset] + data) + '\\n')
            else:
                print info_idx, data
                if float(data[info_idx]) >= float(infoCutoff):
                    outWell_imputed_out.writelines(' '.join([dataset] + data) + '\\n')
                    outWell_imputed_snp_out.writelines(data[1] + '\\n')
                if data[conc_idx] != '-':
                    outSNP_accuracy_out.writelines(' '.join([dataset] + data) + '\\n')
                count += 1
    outWell_imputed_out.close()
    outWell_imputed_snp_out.close()
    outSNP_accuracy_out.close()


args.infoFiles = "${infos}"
args.infoCutoff = "${impute_info_cutoff}"
args.outWell_imputed = "${well_out}"
args.outSNP_acc = "${acc_out}"
if args.infoFiles and args.infoCutoff:
    filter_info(args.infoFiles, args.infoCutoff, args.outWell_imputed, args.outSNP_acc)

