#! /usr/bin/python3

import os
import sys
import re
import logging
from argparse import ArgumentParser
import traceback
from vcf.parser import Writer, Reader
import pprint
import time

selffile = os.path.abspath(__file__)
selfdir = os.path.dirname(selffile)
sys.path.append(selfdir +  '/../PGH_API_lims14')
from config import config


def get_logger(log_file):
    logger = logging.getLogger('Convert VCF file to Report format')
    logger.setLevel(logging.DEBUG)

    if not logger.handlers:

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        handler1 = logging.StreamHandler()
        handler1.setLevel(logging.INFO)
        handler1.setFormatter(formatter)
        logger.addHandler(handler1)

        if not (log_file is None):
            handler2 = logging.FileHandler(log_file, mode = 'a')
            handler2.setLevel(logging.DEBUG)
            handler2.setFormatter(formatter)
            logger.addHandler(handler2)

    return logger

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def bcfquery(vcf):
    query1 = ["%CHROM","%POS","%ID","%REF","%ALT","%ALLELE_A","%ALLELE_B", "%GC_SCORE"]
    query2 = ["%GT","%GQ","%IGC","%BAF","%LRR","%X","%Y"]
    query1_out = '\t'.join(query1)
    query2_out = '\t'.join(query2)
    query1_len = len(query1)
    query2_len = len(query2)
    qout = vcf + '.format.txt'
    try:
        cmd = os.system("%s query -H -f '%s[\t%s]\n' %s -o %s" % (bcftools, query1_out, query2_out, vcf, qout))
    except Exception as e:
        logger.error(e.args)
        logger.debug(traceback.format_exc(e))

    if cmd != 0:
        logger.error("cmd :" + str(cmd))
        sys.exit(1)

def getindex(indexfile):
    indexh = open(indexfile, 'r')
    indexdict = {}
    for line in indexh:
        line = line.strip()
        cols = line.split('\t')
        indexdict[cols[0]] = cols[3]
    return indexdict


def to_report(vcf, outh, samples, indexdict):
    qout = vcf + '.format.txt'
    vcfh = open(qout, 'r')
    linecount = 0
    hdict1 = {}
    hdict2 = {}
    for line in vcfh:
        line = line.strip()
        cols = line.split('\t')
        if linecount == 0:
            for cv in cols:
                searchobj1 = re.search(r'\[(\d+)\]([A-Z_]+)$', cv)
                searchobj2 = re.search(r'\[(\d+)\]([A-Z0-9_]+):([A-Z_]+)', cv)
                if searchobj1:
                    hdict1[searchobj1.group(2)] = int(searchobj1.group(1)) - 1
                elif searchobj2:
                    if searchobj2.group(2) not in hdict2.keys():
                        hdict2[searchobj2.group(2)] = {}
                    hdict2[searchobj2.group(2)][searchobj2.group(3)] = int(searchobj2.group(1)) - 1
        else:
            ref = cols[hdict1["REF"]]
            alt = cols[hdict1["ALT"]].split(',')
            alt1 = ""
            alt2 = ""
            hash_Allele = {}
            #AB
            hash_Allele[cols[hdict1["ALLELE_A"]]] = 'A'
            hash_Allele[cols[hdict1["ALLELE_B"]]] = 'B'
            
            if len(alt) == 1:
                alt1 = alt[0]
                alt2 = ""
            elif len(alt) == 2:
                alt1 = alt[0]
                alt2 = alt[1]
            # ref alt
            hash_base = {}
            if cols[hdict1["ALLELE_A"]] == '0':
                hash_base['0'] = ref
            elif cols[hdict1["ALLELE_A"]] == '1':
                hash_base['1'] = alt1
            elif cols[hdict1["ALLELE_A"]] == '2':
                hash_base['2'] = alt2
            if cols[hdict1["ALLELE_B"]] == '0':
                hash_base['0'] = ref
            elif cols[hdict1["ALLELE_B"]] == '1':
                hash_base['1'] = alt1
            elif cols[hdict1["ALLELE_B"]] == '2':
                hash_base['2'] = alt2
            #for i in range(query1_len, len(cols), query2_len):
            for key in samples:
                if cols[hdict2[key]["GT"]] == "./.":
                    rp_ref = "-"
                    rp_alt = "-"
                    rp_allA = "-"
                    rp_allB = "-"
                elif cols[hdict2[key]["GT"]] == "0/1":
                    rp_ref = hash_base[cols[hdict1["ALLELE_A"]]]
                    rp_alt = hash_base[cols[hdict1["ALLELE_B"]]]
                    rp_allA = 'A'
                    rp_allB = 'B'
                elif cols[hdict2[key]["GT"]] == "0/0":
                    rp_ref = hash_base['0']
                    rp_alt = hash_base['0']
                    rp_allA = hash_Allele['0']
                    rp_allB = hash_Allele['0']
                elif cols[hdict2[key]["GT"]] == "1/1":
                    rp_ref = hash_base['1']
                    rp_alt = hash_base['1']
                    rp_allA = hash_Allele['1']
                    rp_allB = hash_Allele['1']
                elif cols[hdict2[key]["GT"]] == "1/2":
                    rp_ref = hash_base[cols[hdict1["ALLELE_A"]]]
                    rp_alt = hash_base[cols[hdict1["ALLELE_B"]]]
                    rp_allA = 'A'
                    rp_allB = 'B'
                elif cols[hdict2[key]["GT"]] == "2/2":
                    rp_ref = hash_base['2']
                    rp_alt = hash_base['2']
                    rp_allA = hash_Allele['2']
                    rp_allB = hash_Allele['2']
                else:
                    print("Error : unexpected genotype")
                    sys.exit(-1)
                rp_gc = cols[hdict2[key]["IGC"]]
                rp_baf = cols[hdict2[key]["BAF"]]
                rp_lrr = cols[hdict2[key]["LRR"]]
                rp_x = cols[hdict2[key]["X"]]
                rp_y = cols[hdict2[key]["Y"]]
                if rp_ref != "":
                    outh.write('\t'.join([cols[hdict1["ID"]], key, rp_ref, rp_alt, rp_gc, cols[hdict1["CHROM"]], cols[hdict1["POS"]], rp_x, rp_y, "-", "-", rp_baf, rp_lrr, "-", "-", rp_allA, rp_allB, indexdict[cols[hdict1["ID"]]], cols[hdict1["GC_SCORE"]], "-", "-", "-"]))
                    outh.write('\n')
        if linecount % 200000 == 0:
            logger.info("Step3: %s %d line completed" % (vcf, linecount))
        linecount += 1
    vcfh.close()

def main():
    parser = ArgumentParser(description="Convert VCF file to Report format")
    parser.add_argument("--vcf", dest="vcf", nargs="+", required=True, help="One or more VCF files to process")
    parser.add_argument("--outfile", dest="outfile", required=True, help="Path for generation of Report output")
    parser.add_argument("--log", dest="log_file", default=None, required=False, help="File to write logging information (optional)")
    args = parser.parse_args()
    
    global logger, indexfile, bcftools
    indexfile = config('database', 'indexfile')
    bcftools = config('software', 'bcftools')
    logger = get_logger(args.log_file)
    
    #####VCF Header => Report Header #####
    vcfsamples = AutoVivification()
    samples = []
    outh = open(args.outfile, 'w')
    bpm = ""
    version = ""
    command = ""
    date = ""
    for vcf in args.vcf:
        vcfsamples[vcf] = []
        vcfh = open(vcf, 'r')
        for cc in range(1, 150):
            vcfline = vcfh.readline()
            vcfline = vcfline.strip()
            linesearch1 = re.match("##BPM=(.*)", vcfline)
            linesearch2 = re.match("##bcftools.*gtc2vcfVersion=(.*)", vcfline)
            linesearch3 = re.match("##bcftools.*gtc2vcfCommand=(.*) Date=(.*)", vcfline)
            linesearch4 = re.match("#CHROM", vcfline)
            if linesearch1 and bpm == "":
                bpm = linesearch1.group(1)
            elif linesearch2 and version == "":
                version = linesearch2.group(1)
            elif linesearch3 and command == "":
                command = linesearch3.group(1)
                date = linesearch3.group(2)
            elif linesearch4:
                cols = vcfline.split('\t')
                for v in cols[9:]:
                    vcfsamples[vcf].append(v)
                    samples.append(v)
        vcfh.close()
    pprint.pprint(vcfsamples)

    ####out head#####
    outh.write("[Header]\n")
    outh.write("bcftools_+gtc2vcfVersion %s\n" % (version))
    outh.write("Processing Date          %s\n" % (date))
    outh.write("Content                  %s\n" % bpm)
    outh.write("Num SNPs                 -\n")
    outh.write("Total SNPs               -\n")
    outh.write("Num Samples              %s\n" % len(samples))
    outh.write("Total Samples            %s\n" % len(samples))
    outh.write("[Data]\n")
    outh.write("SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tChr\tPosition\tX\tY\tX Raw\tY Raw\tB Allele Freq\tLog R Ratio\tCNV Value\tCNV Confidence\tAllele1 - AB\tAllele2 - AB\tSNP Index\tGT Score\tCluster Sep\tTheta\tR\n")
    logger.info("Step1 : vcf head  ==> report head completed")

    #####index get#####
    indexdict = getindex(indexfile)

    #####process#####
    for vcf in args.vcf:
        bcfquery(vcf)    #####bcftools query
        logger.info("Step2 : %s bcfquery completed" % vcf)
        to_report(vcf, outh, vcfsamples[vcf], indexdict)  #####vcf to report
        logger.info("End : %s" % vcf)
    outh.close()

if __name__ == "__main__":
    main()
