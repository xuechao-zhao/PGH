#! /usr/bin/python
# _*_ coding: utf-8 _*_

import os
import sys
import pprint
import numpy as np

def logRdev(arr):
    std = np.std(arr)
    return std

def callrate(allele1, allele2):
    gt_status = 0
    if allele1 != "-" and allel2 != '-':
        gt_status = 1
    return gt_status


def out_qc(filename, outfile):

    #读入位点信息
    file_handle = open(filename, 'r')
    out = open(outfile, 'w')
    line_count = 0
    heads = []
    info_dict = {}
    logRatio_dict = {}
    callrate_dict = {}
    callrate_dict_all = {}
    for line in file_handle:
        line = line.strip()
        line_count += 1
        cols = line.split('\t')
        if line_count < 10:
            continue
        elif line_count == 10:
            heads = cols
        else:
            try:
                chipid = cols[heads.index('Sample ID')]
                chrom = cols[heads.index('Chr')]
                allele1 = cols[heads.index('Allele1 - AB')]
                allele2 = cols[heads.index('Allele2 - AB')]
                logRatio = cols[heads.index('Log R Ratio')]
            except Exception as e:
                print('-----PGH_QC.py-----')
                traceback.print_exc()
                print(type(e), e)
                exit(1)
            if chipid not in logRatio_dict.keys():
                logRatio_dict[chipid] = []
                callrate_dict[chipid] = 0
                callrate_dict_all[chipid] = 0
            callrate_dict_all[chipid] +=1
            if allele1 != "-" and allele2 != "-":
                callrate_dict[chipid] += 1            
            if chrom != "0" and chrom != "MT" and chrom != "X" and chrom != "Y" and chrom != "XY":
                if logRatio != "NaN" and logRatio != "nan" and logRatio != "Nan" and logRatio != 'inf':
                    logRatio_dict[chipid].append(float(logRatio))

    #pprint.pprint(callrate_dict)
    #pprint.pprint(logRatio_dict)

    for chipid in callrate_dict.keys():
        out.write('%s\t%.4f\t%.2f\n' % (chipid, float(callrate_dict[chipid])/float(callrate_dict_all[chipid]), logRdev(logRatio_dict[chipid])))
                
if __name__ == '__main__':
    out_qc(sys.argv[1], sys.argv[2])
