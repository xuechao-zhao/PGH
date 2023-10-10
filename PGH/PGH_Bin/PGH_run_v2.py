#! /usr/bin/python3

import os
import re
import sys
import pprint
import logging
from argparse import ArgumentParser
import traceback
import time
from config import config
import shutil
import math
from multiprocessing import Pool
import getpass
from PGH_QC_V1 import PGH_QC_report
VERSION = "1.0.0"


def get_logger(log_file):
    logger = logging.getLogger('PGH RUN')
    logger.setLevel(logging.DEBUG)

    if not logger.handlers:
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        handler1 = logging.StreamHandler()
        handler1.setLevel(logging.DEBUG)
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

def check_data_exists(sampleidlist):
    try:
        pathdict = AutoVivification()
        for home, dirs, files in os.walk(datadir):
            for filename in files:
                for sampleid in sampleidlist:
                    if sampleid not in pathdict[family].keys():
                        pathdict[family][sampleid]['type1'] = []
                        pathdict[family][sampleid]['type2'] = []
                    if filename.startswith(sampleid) and filename.endswith(postfix1):
                        fullname = os.path.join(home, filename)
                        pathdict[family][sampleid]['type1'].append(fullname)
                    elif filename.startswith(sampleid) and filename.endswith(postfix2):
                        fullname = os.path.join(home, filename)
                        pathdict[family][sampleid]['type2'].append(fullname)
                    else:
                        pass
        for sampleid in pathdict[family].keys():
            if len(pathdict[family][sampleid]['type1']) == int(mincount) and len(pathdict[family][sampleid]['type2']) == int(mincount):
                #post('log', family, '%s 原始数据查找成功' % sampleid, logfile, 0)
                logger.info('%s 原始数据查找成功' % sampleid)
            else:
                raise Exception("测序数据1与测序数据2文件数有误")
        return pathdict
    except Exception as e:
        #post('log', family, '原始数据查找失败: %s' % e.args, logfile, 0)
        logger.info('%s 原始数据查找失败 : %s' % (sampleid, e.args))

def read_sheet(sheet):
    sheetdt = AutoVivification()
    idmapdt = AutoVivification()
    refdt = AutoVivification()
    DNA = []
    MDA = []

    sheetfd = open(sheet, 'r')
    for line in sheetfd:
        line = line.strip()
        chipid, sampleid, samplerole, samplename, labid, sampletype = line.split('\t')

        if sampletype == "DNA" and chipid not in DNA:
            DNA.append(chipid)
        elif sampletype == "MDA" and chipid not in MDA:
            MDA.append(chipid)

        if sampleid != 'reference':
            idmapdt[sampleid] = chipid
            sheetdt[chipid]['sampleid'] = sampleid
            sheetdt[chipid]['samplerole'] = samplerole
            sheetdt[chipid]['samplename'] = samplename
            sheetdt[chipid]['labid'] = labid
            sheetdt[chipid]['sampletype'] = sampletype
        else:
            refdt[chipid]['samplerole'] = samplerole
            refdt[chipid]['samplename'] = samplename
            refdt[chipid]['labid'] = labid
            refdt[chipid]['sampletype'] = sampletype

    for refid in refdt.keys():
        OUTH = open('%s/samplesheet.%s.txt' % (outdir, refid), 'w')
        for chipid in sheetdt.keys():
            OUTH.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chipid, sheetdt[chipid]['sampleid'], sheetdt[chipid]['samplerole'], sheetdt[chipid]['samplename'], sheetdt[chipid]['labid'], sheetdt[chipid]['sampletype']))
        OUTH.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (refid, 'reference', refdt[refid]['samplerole'], refdt[refid]['samplename'], refdt[refid]['labid'], refdt[refid]['sampletype']))
        OUTH.close()

    return sheetdt, refdt, idmapdt, DNA, MDA


def get_gtc(pathdict):
    try:
        threads = len(pathdict[family].keys())
        if threads > 10:
            threads = 10 #控制进程数目
        ##创建目录
        if os.path.exists(chipdir):
            shutil.rmtree(chipdir)
        os.mkdir(chipdir)
        if os.path.exists(reportdir):
            shutil.rmtree(reportdir)
        os.mkdir(reportdir)

        ##链接数据
        for chipid in pathdict[family]:
            for datatype in pathdict[family][chipid]:
                for datafile in pathdict[family][chipid][datatype]:
                    os.system('ln -s %s %s' % (datafile, chipdir))

        ##idat转gtc
        write_log('ALL', 'IDAT2GTC', 'Start')
        cmd = os.system('%s gencall %s %s %s -f %s -t %s -g' % (iaapcli, bpm, egt, reportdir, chipdir, threads))
        if cmd == 0:
            logger.info('IDAT --> GTC 成功')
            write_log('ALL', 'IDAT2GTC', 'Finish')
        else:
            raise Exception('IDAT --> GTC 失败')
        
    except Exception as e:
        logger.error(e.args)
        logger.error('=====================================================')
        logger.error(traceback.format_exc())
        logger.error('=====================================================')
        exit(1)
        
def get_vcf_one(chipid):
    try:
        write_log(chipid, 'GTC2VCF', 'Start')
        gtc = '%s/%s.gtc' % (reportdir, chipid)
        gtclist = '%s/%s.gtc.list' % (reportdir, chipid)
        vcf = '%s/%s.vcf' % (reportdir, chipid)
        gtcfd = open(gtclist, 'w')
        gtcfd.write(gtc + '\n')
        gtcfd.close()
        cmd = os.system('%s +gtc2vcf -b %s -e %s -f %s --threads 2 -o %s/%s.vcf -O v -g %s/%s.gtc.list' % (bcftools, bpm, egt, fasta, reportdir,chipid, reportdir, chipid))
        if cmd == 0:
            logger.info('GTC --> VCF成功[ %s ]' % chipid)
            write_log(chipid, 'GTC2VCF', 'Finish')
        else:
            raise Exception('GTC --> VCF失败[ %s ]' % chipid)
        return '%s : AAA' % chipid
    except Exception as e:
        logger.error(e.args)
        logger.error('=====================================================')
        logger.error(traceback.format_exc())
        logger.error('=====================================================')
        exit(1)

def get_report_one(chipid):
    try:
        write_log(chipid, 'VCF2Report', 'Start')
        vcf = '%s/%s.vcf' % (reportdir, chipid)
        report = '%s/%s.report.txt' % (reportdir, chipid)
        cmd = os.system('python3 %s/vcf2report.py --vcf %s/%s.vcf --outfile %s/%s.report.txt' % (script, reportdir, chipid, reportdir, chipid))
        if cmd == 0:
            logger.info('VCF --> Report成功[ %s ]' % chipid)
            write_log(chipid, 'VCF2Report', 'Finish')
        else:
            raise Exception('VCF --> Report失败[ %s ]' % chipid)
        return '%s : BBB' % chipid
    except Exception as e:
        logger.error(e.args)
        logger.error('=====================================================')
        logger.error(traceback.format_exc())
        logger.error('=====================================================')
        exit(1)



def gender(Xhet, Xhom, Ynumber):
    gender_status = 0
    gender_t1 = 0
    gender_t2 = 0
    #判断X染色体
    #print(str(Xhet)+'\t'+str(Xhom)+'\t'+str(Ynumber)+'\n')
    if Xhet + Xhom > 20000:
        if float(Xhet)/float(Xhom) > 0.1:
            gender_t1 = 2
        else:
            gender_t1 = 1
    else:
        gender_t1 = 0
    #判断Y染色体
    if Ynumber > 3500:
        gender_t2 = 1
    else:
        gender_t2 = 2
    #结合X染色体与Y染色体
    if gender_t1 == 2 and gender_t2 == 2:
        gender_status = 'female'
    elif gender_t1 == 1 and gender_t2 == 1:
        gender_status = 'male'
    else:
        gender_status = 'unknown'

    return gender_status

def get_gender(report_file):
    hash_gender = {}
    file = open(report_file,'r')
    line_count = 0
    heads = []
    for line in file.readlines():
        line = line.strip()
        data = line.split('\t')
        line_count +=1
        if line_count <= 10:
            heads = data
        else:
            chr = data[heads.index('Chr')]
            sample_ID = data[heads.index('Sample ID')]
            Allele1 = data[heads.index('Allele1 - Top')]
            Allele2 = data[heads.index('Allele2 - Top')]
            if sample_ID not in hash_gender:
                hash_gender[sample_ID]={}
                hash_gender[sample_ID]['hetX'] = 0
                hash_gender[sample_ID]['homX'] = 0
                hash_gender[sample_ID]['Ynum'] = 0
            if chr == 'X' and Allele1 != '-' and Allele2 != '-':
                if Allele1 != Allele2:
                    hash_gender[sample_ID]['hetX']+=1
                else:
                    hash_gender[sample_ID]['homX']+=1
            if chr == 'Y' and Allele1 != '-' and Allele2 != '-':
                hash_gender[sample_ID]['Ynum']+=1    
    file.close()
    hash_result = {}
    for kk in hash_gender:
        hetx = hash_gender[kk]['hetX']
        hmox = hash_gender[kk]['homX']
        Ynum = hash_gender[kk]['Ynum']
        key_gender = gender(hetx,hmox,Ynum)
        hash_result[kk] = key_gender 
    return hash_result



def readbreak():
    breakdict = AutoVivification()
    BREAKH = open(breakfile, 'r')
    for line in BREAKH:
        line = line.strip()
        chrom, band, size = line.split('\t')
        size_num = 5000000
        matchObj = re.match(r'(\d+)M', size)
        if matchObj:
            size_num = int(matchObj.group(1))*1000000
        breakdict[chrom][band] = size_num
    BREAKH.close()
    return breakdict

def write_log(chipid, module, state):
    SLH = open(steplog, 'a')
    dateinfo = os.popen('echo `date`')
    date = dateinfo.read().split('\n')[0]
    SLH.write('[Process:] %s %s %s %s\n' % (date, chipid, module, state))
    SLH.close()

def main():
    parser = ArgumentParser(description = "PGH RUN")
    parser.add_argument("--family", dest="family", required=True, help="input family name") ###nargs="+"
    parser.add_argument("--sheet", dest="sheet", required=True, help="input samplesheet.txt") ###nargs="+"
    parser.add_argument("--outdir", dest="outdir", required=True, help="outdir")
    parser.add_argument("--step", dest="step", required=True, help="step", choices=['all', 'report'])
    parser.add_argument("--breakfile", dest="breakfile", required=False, help="breakpoint file")
    parser.add_argument("--log", dest="log_file", default=None, required=False, help="File to write logging information (optional)")
    args = parser.parse_args()

    ##读取配置文件
    global bpm,egt,fasta,iaapcli,bcftools,datadir,script,postfix1,postfix2,mincount,gc,cnvscript_embryo,cnvscript_dna,hg19_dict,cytoband,site1,site2
    bpm = config('database', 'bpm')
    egt = config('database', 'egt')
    fasta = config('database', 'fasta')
    hg19_dict = config('database', 'hg19_dict')
    cytoband = config('database', 'cytoband')
    gc = config('database', 'gc20kb')
    iaapcli = config('software', 'iaap-cli')
    bcftools = config('software', 'bcftools')
    datadir = config('path', 'data')
    script = config('path', 'script')
    postfix1 = config('basic', 'postfix1')
    postfix2 = config('basic', 'postfix2')
    mincount = config('basic', 'mincount')

    global family, outdir, chipdir, reportdir, logger, breakfile, steplog
    family = args.family
    outdir = os.path.abspath(args.outdir)
    step = args.step
    sheet = os.path.abspath(args.sheet)
    if args.log_file:
        logfile = os.path.abspath(args.log_file)
    else:
        logfile = args.log_file
    if args.breakfile:
        breakfile = os.path.abspath(args.breakfile)
    else:
        breakfile = args.breakfile
    chipdir = outdir + '/Chip'
    reportdir = outdir + '/Report'
    steplog = outdir + '/step.%s.log' % step
    logger = get_logger(logfile)
    
    sheetdt, refdt, idmapdt, DNA, MDA = read_sheet(sheet)
    pprint.pprint(sheetdt)
    pprint.pprint(refdt)
    pprint.pprint(idmapdt)
    pprint.pprint(DNA)
    pprint.pprint(MDA)

    chipids = DNA + MDA
    chipids = list(set(chipids))
    pathdict = check_data_exists(chipids)
    pprint.pprint(pathdict)
    threads = len(pathdict[family].keys())
    if threads > 10:
        threads = 10 #控制进程数目
    if step == 'all' or step == 'report':
        ##IDAT转GTC
        get_gtc(pathdict)
        
        ##GTC转VCF
        p = Pool(threads)
        for chipid in pathdict[family].keys():
            p.apply_async(get_vcf_one, args=(chipid,))
        p.close()
        p.join()
        
        ##VCF转Report
        p = Pool(threads)
        for chipid in pathdict[family].keys():
            p.apply_async(get_report_one, args=(chipid,))
        p.close()
        p.join()
    
        vcfs = ''
        for chipid in pathdict[family].keys():
            vcfs = vcfs + ' %s/%s.vcf' % (reportdir, chipid)
        cmd = os.system('python3 %s/vcf2report.py --vcf %s --outfile %s/%s.report.txt' % (script, vcfs, reportdir, family))
        if cmd == 0:
            logger.info('VCF --> Report成功[ %s ]' % family)

        cmd = os.system('python3 %s/PGH_QC.py %s/%s.report.txt %s/%s.report.qc.txt' % (script, reportdir, family, reportdir, family))
        if cmd == 0:
            logger.info('Report --> QC成功[ %s ]' % family)


if __name__ == "__main__":
    main()









