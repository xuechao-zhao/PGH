#!usr/bin/env python
#coding:utf-8
#version:1.0
#从LIMS系统获取信息生成质控报告
import re,sys,os
import pandas as pd
import numpy as np
import pprint
import logging
from argparse import ArgumentParser
import traceback
import datetime
from config import config
import requests
import base64
import json
import time
import glob

def get_logger(log_file):
    logger = logging.getLogger('LIMS Result Post')
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




def cal_QC(report,sample,cnvfile,outfile):
    df_report = pd.read_table(report,sep='\t',header=9,low_memory=False)
    hash_log = {}
    hash_call = {}
    sample_set = set()
    #read cnvfile
    hash_chr = {}
    switch1 = 0
    switch2 = 0
    mother_id = ''
    father_id = ''
    embryo = set()
    if cnvfile != 'NA':
        file = open(cnvfile,'r')
        for line in file.readlines():
            line = line.strip()
            data = line.split('\t')
            if data[0] not in hash_chr.keys():
                hash_chr[data[0]] = set()
            hash_chr[data[0]].add(data[1])
        file.close()
    #read sample file
    file = open(sample,'r')
    result = open(outfile,'w+')
    result.write('Chip_ID\tCall_Rate\tlogRdev\tGender\tBAF_SD\tBAD_MAPD\tADO\tADI\tHetero_loci\n')
    for line in file.readlines():
        line = line.strip()
        data = line.split('\t')
        sample_set.add(data[0])
        if data[1] == 'mother':
            mother_id = data[0]
            switch1 = 1
        elif data[1] == 'father':
            father_id = data[0]
            switch2 = 1
        if data[2] == 'embryo':
            embryo.add(data[0])
    file.close()
    
    temp_data1 = []
    for x in range(1,23):
        m = str(x)
        temp_data1.append(m)
        
    #mother and father allele
    df_mother = pd.DataFrame(df_report[(df_report['Sample ID'] == mother_id) & (df_report['Allele1 - Top'] != '-') & (df_report['Allele2 - Top'] !=
 '-') & (df_report['Allele1 - Top'] == df_report['Allele2 - Top'])],columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
    df_mother.columns= ['SNP Name','Chr','Position','mother_Allele1','mother_Allele2']
    df_father = pd.DataFrame(df_report[(df_report['Sample ID'] == father_id)& (df_report['Allele1 - Top'] != '-') & (df_report['Allele2 - Top'] !=
 '-') & (df_report['Allele1 - Top'] == df_report['Allele2 - Top'])],columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
    df_father.columns= ['SNP Name','Chr','Position','father_Allele1','father_Allele2']
    df_parent = pd.merge(df_father,df_mother,on =['SNP Name','Chr','Position'],how = 'inner')
    for kk in sample_set:
        df1 = df_report[df_report['Sample ID'] == kk]
        ALL_num = df1['SNP Name'].count()
        if ALL_num == 0:
            sys.exit()
        df3 = df1[(df1['Log R Ratio'] != 'NaN') & (df1['Log R Ratio'] != 'nan') & (df1['Log R Ratio'] != 'NA') & (df1['Chr'] !='0') & (df1['Chr'] !='X') & (df1['Chr'] !='Y') & (df1['Chr'] !='XY') & (df1['Chr'] !='MT')]
        logR = df3['Log R Ratio'].std()
        logR = '%.2f' %(logR)
        df2 = df1[(df1['Allele1 - AB']!='-') & (df1['Allele2 - AB']!='-')]
        call = df2['SNP Name'].count()
        rate1 = float(call)/float(ALL_num)
        rate1 = '%.4f'%(rate1)
        Ynum = df2[df2['Chr'] == 'Y']['SNP Name'].count()
        hetx = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] != df2['Allele2 - AB'])]['SNP Name'].count()
        hmox = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] == df2['Allele2 - AB'])]['SNP Name'].count()
        key_gender = gender(hetx,hmox,Ynum)
        df_autosome = df2[(df2['Chr'] !='0') & (df2['Chr'] !='X') & (df2['Chr'] !='Y') & (df2['Chr'] !='XY') & (df2['Chr'] !='MT')]
        #计算杂合位点数目
        temp=[]
        for key in temp_data1:
            df_a1 = df_autosome[(df_autosome['Chr'] == key) & (df_autosome['Allele2 - AB'] != df_autosome['Allele1 - AB'])]
            count = df_a1['SNP Name'].count()
            temp2 = str(key)+':'+str(count)
            temp.append(temp2)
        value1 = '{'+kk+':{'+','.join(temp)+'}}'
        #计算BAF离散程度以及ADO/ADI
        df_filter = df_autosome
        if kk in hash_chr.keys():
            region = hash_chr[kk]
            for key in region:
                if not re.findall('-',key):
                    df_filter = df_filter[df_filter['Chr'] != key]
                else:
                    temp1 = key.split(':')
                    chr1 = temp1[0]
                    temp2 = temp1[1].split('-')
                    start1 = int(temp2[0])
                    end1 = int(temp2[1])
                    df_filter.loc[(df_filter.Chr == chr1) & (df_filter.Position >= start1) & (df_filter.Position <= end1),'Allele1 - AB'] = '-'

        #计算BAF
        df_het = df_filter[(df_filter['Allele1 - AB'] != '-') & (df_filter['Allele2 - AB'] != '-') & (df_filter['Allele1 - AB'] != df_filter['Allele2 - AB'])]
        df_BAF = df_het[(df_het['B Allele Freq'] != 'Nan') & (df_het['B Allele Freq'] != 'nan')]
        #计算MAPD
        temp_MAF = 0
        MAPD_data = []
        df_BAF['BAF'] = df_BAF['B Allele Freq']
        for key in df_BAF.itertuples():
            BAF_Value = getattr(key,'BAF')
            BAF1 = float(BAF_Value)-float(temp_MAF)
            if BAF1 < 0:
                BAF1 = -BAF1
            MAPD_data.append(BAF1)
            temp_MAF = BAF_Value
        BAF_SD = df_BAF['B Allele Freq'].std()
        BAF_SD = '%.4f'%(BAF_SD)
        MAPD = np.median(MAPD_data)
        MAPD = '%.4f'%(MAPD)
        if kk not in embryo:
            result.write(kk+'\t'+str(rate1)+'\t'+str(logR)+'\t'+key_gender+'\t'+str(BAF_SD)+'\t'+str(MAPD)+'\tNA\tNA\t'+value1+'\n')
            continue
            
        #计算ADO/ADI
        if switch1 != 1 or switch2 != 1:
            continue
        df_embryo = df_filter[(df_filter['Allele1 - AB'] != '-') & (df_filter['Allele2 - AB'] != '-')]
        df_embryo = pd.DataFrame(df_embryo,columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
        df_embryo.columns= ['SNP Name','Chr','Position','embryo_Allele1','embryo_Allele2']
        df_ref = pd.merge(df_parent,df_embryo,on =['SNP Name','Chr','Position'],how = 'inner')
        num_ADO_ALL = 0
        num_ADO = 0
        num_ADI_ALL = 0
        num_ADI = 0
        
        for key in df_ref.itertuples():
            mother_Allele1 = getattr(key,'mother_Allele1')
            mother_Allele2 = getattr(key,'mother_Allele2')
            father_Allele1 = getattr(key,'father_Allele1')
            father_Allele2 = getattr(key,'father_Allele2')
            embryo_Allele1 = getattr(key,'embryo_Allele1')
            embryo_Allele2 = getattr(key,'embryo_Allele2')
            if mother_Allele1 != father_Allele1:
                num_ADO_ALL +=1
                if embryo_Allele1 == embryo_Allele2:
                    num_ADO +=1
            if mother_Allele1 == father_Allele1:
                num_ADI_ALL +=1
                if embryo_Allele1 != mother_Allele1 or embryo_Allele2 != mother_Allele1:
                    num_ADI +=1
        ADO = 'Nan'
        ADI = 'Nan'
        if num_ADO_ALL != 0:
            ADO = float(num_ADO)/float(num_ADO_ALL)*100
            ADO = '%.2f'%(ADO)
        if num_ADI_ALL != 0:
            ADI = float(num_ADI)/float(num_ADI_ALL)*100
            ADI = '%.2f'%(ADI)
        result.write(kk+'\t'+str(rate1)+'\t'+str(logR)+'\t'+key_gender+'\t'+str(BAF_SD)+'\t'+str(MAPD)+'\t'+str(ADO)+'\t'+str(ADI)+'\t'+value1+'\n')
    result.close()
    
def get_info(file1):
    hash2={}
    file = open(file1,'r')
    hash_loci = {}
    hash_base={}
    hash_id = {}
    temp_data = []
    for x in range(1,23):
        m = 'chr'+str(x)
        temp_data.append(m)
    for line in file.readlines():
        line = line.strip()
        data = line.split('\t')
        if 'mother' in data[2] or 'father' in data[2] or 'reference' in data[2]:
            continue
        if data[0] not in temp_data:
            continue
        sample_id = re.findall('(.*?)\_[F|M]',data[2])[0]
        key_Value = data[0]+'\t'+data[1]+'\t'+data[2]
        key_id = data[0]+'\t'+data[1]+'\t'+sample_id
        hash_id[key_id]=1
        hash_base[key_Value] = data[3]
        hash_loci[key_Value] = data[4]
    file.close()
    hash_num={}
    
    for key in hash_id.keys():
        key1 = key+'_F'
        key2 = key+'_M'
        temp = key.split('\t')
        embryo_id = temp[2]
        base1 = hash_base[key1]
        base2 = hash_base[key2]
        if base1 == base2:
            continue
        value1 = hash_loci[key1]
        value2 = hash_loci[key2]
        switch1 = 0
        switch2 = 0
        if value1 == 'unphased' or value1 == 'grey':
            switch1 = 1
        if value2 == 'unphased' or value2 == 'grey':
            switch1 += 1
        if value1 !='unphased' and value1 != 'grey':
            switch2 =1
        if value2 !='unphased' and value2 != 'grey':
            switch2 += 1
        if temp[2] not in hash_num.keys():
            hash_num[temp[2]]={}

        if switch1 == 1 and switch2 ==1:
            temp[0] = temp[0].replace('chr','')
            temp[0] = int(temp[0])
            if temp[0] not in hash_num[temp[2]].keys():
                hash_num[temp[2]][temp[0]]=1
            else:
                hash_num[temp[2]][temp[0]] = hash_num[temp[2]][temp[0]]+1
    
    
    return hash_num


def PGH_QC_report(json_file,report_file,sheet,user_name,merge1_file,chrX_merge1_file,cnv_merge_report_file,outdir,log_file,script,index):
    try:
        logger = get_logger(log_file)
        date = datetime.datetime.now().strftime('%Y.%m.%d')
        index+=1
        hash_result = {}
        #获得Call_Rate,logR和染色体结果
        with open(json_file) as fd:
            js = json.load(fd)
        family_path = js['data']['family_path']
        family_name = js['data']['name']
        session_id = js['data']['session']
        Host = config('path', 'url')
        sample = js['data']['sample']
        analyze_id = js['data']['id']
        companyid_name = js['data']['companyid_name']
        serial_number = js['data']['name']
        company_id = js['data']['companyid']
        url = Host + "/pgtsr/result/search"
        sample_chip_position = js['data']['sample_chip_position']
        hash_sample_type = {}
        hash_sample_role = {}
        hash_CNV = {}
        hash_sample = {}
        #获取样本的CNV结果
        file_CNV = open(cnv_merge_report_file,'r')
        for line_CNV in file_CNV.readlines():
            line_CNV = line_CNV.strip()
            data = line_CNV.split('\t')
            if data[0] not in hash_CNV.keys():
                hash_CNV[data[0]] = []
            hash_CNV[data[0]].append(data[1])
        file_CNV.close()
        #获取样本的sample名
        #找出重复的sibling和embryo
        file_sample = open(sheet,'r')
        hash_tags = {}
        for line_sample in file_sample.readlines():
            line_sample = line_sample.strip()
            data = line_sample.split('\t')
            if data[0] not in hash_tags.keys():
                hash_tags[data[0]] = 1
            else:
                hash_tags[data[0]] += 1
        file_sample.close()
        #计算有效位点数目
        hash_info_num = get_info(merge1_file)
        file_sample = open(sheet,'r')
        hash_info = {}
        for line_sample in file_sample.readlines():
            line_sample = line_sample.strip()
            data = line_sample.split('\t')
            if hash_tags[data[0]] == 2 and data[2] == 'sibling':
                continue
            hash_sample[data[0]] = data[1]+'_'+data[2]
            hash_sample_role[data[0]] = data[2]
            value1 = 'NA'
            if data[1] in hash_info_num.keys():
                value1 = hash_info_num[data[1]]
            value1 = str(value1)
            value1 = value1.replace('chr','')
            hash_info[data[0]] = '{'+data[0]+':{'+value1+'}}'
        file_sample.close()
        for kk in js['data']['sample']:
            chip= kk['chip_position']
            back_true = kk['is_pass_back']
            ana_role = kk['ana_role']
            if ana_role == 'mother':
                family_name = kk['sj_name']#获得家系名
            sample_type = kk['smp_ana_type']
            hash_sample_type[chip] = sample_type
            if back_true == False:
                del sample_chip_position[chip]#删除FALSE的芯片ID
        for chipid in sample_chip_position:
            sample_type = hash_sample_type[chipid]
            sample_ID = sample_chip_position[chipid]
            payload = {'company_id': company_id,'sample_id': sample_ID}
            response = requests.request("POST", url, headers={'Cookie': session_id}, data=payload)
            outcode = json.loads(response.text)['http_code']
            if outcode == 200: #表明该样本的结果已经存在，获取基本质控值
                result_text = json.loads(response.text)
                chip_id = result_text['data'][0]['chip_position'][1]
                call_rate = result_text['data'][0]['char4']
                logR = result_text['data'][0]['char5']
                cnv_result = result_text['data'][0]['char10']
                trans = result_text['data'][0]['char13']
                remark = result_text['data'][0]['char15']
                if chip_id not in hash_result.keys():
                    hash_result[chip_id] = {}
                cnv_result = cnv_result.replace('Euploid','正常')
                hash_result[chip_id]['call_rate'] = call_rate
                hash_result[chip_id]['logR_dev'] = logR
                hash_result[chip_id]['Transplant_Result'] = remark #移植结果
                hash_result[chip_id]['Phenotype'] = cnv_result #CNV分析结果
                hash_result[chip_id]['sampletype']= hash_sample_type[chipid] #样本类型
                hash_result[chip_id]['Result']= trans #分型结果
        #计算ADO/ADI，female
        hash_ado = get_ADO_ADI_female(report_file,sheet,cnv_merge_report_file,outdir,log_file)
        #亲缘关系鉴定
        command = "python3 %s/%s %s %s"%(script,'array_select_sample_map.py', report_file, sheet)
        os.system(command)
        file_sample = os.path.basename(sheet+'_FinalReport.txt')
        command = "python %s/%s %s %s %s"%(script,'array_kinship.py', outdir,file_sample, os.path.basename(sheet))
        os.system(command)
        file_kin = os.path.join(outdir,sheet+'_FinalReport.txt.kin.json.txt')
        file_kin_R = open(file_kin,'r')
        json_txt = ''
        for line in file_kin_R.readlines():
            line = line.strip()
            json_txt = line
        kin_file = os.path.join(outdir,sheet+'_FinalReport.txt.kin.txt')
        
        #生成QC.txt
        fileout = open(report_file+'.chipid','w+')
        for chip_id in sample_chip_position:
            outfile = os.path.join(outdir,chip_id+'.QC.txt')
            result = open(outfile,'w+')
            chip = re.findall('(.*)\_',chip_id)[0]
            ado = hash_ado[chip_id]['ado']
            adi = hash_ado[chip_id]['adi']
            BAF_SD = hash_ado[chip_id]['BAF_SD']
            BAF_MAPD = hash_ado[chip_id]['BAF_MAPD']
            gender = hash_ado[chip_id]['gender']
            sampletype = hash_result[chip_id]['sampletype']
            if sampletype == 'DNA':
                Phenotype = 'NA'
                pheno_result = 'NA'
                Transplant_Result = 'NA'
            else:
                Phenotype = hash_result[chip_id]['Phenotype']
                pheno_result = hash_result[chip_id]['Result'] #分型结果
                Transplant_Result = hash_result[chip_id]['Transplant_Result']
            Remark = '家系样本'
            report = 'TRUE'
            heteo_num = hash_ado[chip_id]['heteo_num'] 
            relationship = json_txt
            Exp_num = serial_number
            Location = companyid_name
            familyid = family_name
            call_rate = hash_result[chip_id]['call_rate']
            logR_dev = hash_result[chip_id]['logR_dev']
            cnv_result = 'NA'
            if chip_id in hash_CNV.keys():
                cnv_result = ','.join(hash_CNV[chip_id])
                cnv_result = '{'+chip_id+':{'+cnv_result+'}}'
            sample = hash_sample[chip_id]
            samplename = hash_sample_role[chip_id]
            info_num = hash_info[chip_id]
            result.write('sampleid'+ '\t' +chip_id+ '\n')
            result.write('familyid'+ '\t'+ familyid+'\n')
            result.write('samplename'+ '\t'+ samplename+'\n')
            result.write('sampletype'+ '\t'+ sampletype+'\n')
            result.write('chipid'+ '\t'+ chip+'\n')
            result.write('logR_dev'+ '\t'+ logR_dev+'\n')
            result.write('call_rate'+ '\t'+ call_rate+'\n')
            result.write('gender'+ '\t'+ gender+'\n')
            result.write('Result'+ '\t'+ pheno_result+'\n')
            result.write('relationship'+ '\t'+ relationship+'\n')
            result.write('kinfile'+ '\t'+ kin_file+'\n')
            result.write('ado'+ '\t'+ ado+'\n')
            result.write('adi'+ '\t'+ adi+'\n')
            result.write('CNV'+ '\t'+ cnv_result+'\n')
            result.write('sample'+ '\t'+ sample+'\n')
            result.write('Location'+ '\t'+ Location+'\n')
            result.write('Phenotype'+ '\t'+ Phenotype+'\n')
            result.write('Remark'+ '\t'+ Remark+'\n')
            result.write('BAF_MAPD'+ '\t'+ BAF_MAPD+'\n')
            result.write('BAF_SD'+ '\t'+ BAF_SD+'\n')
            result.write('Exp_num'+ '\t'+ Exp_num+'\n')
            result.write('report'+ '\t'+ report+'\n')
            result.write('heteo_num'+ '\t'+ heteo_num+'\n')
            result.write('info_num'+ '\t'+ info_num+'\n')
            result.write('Transplant_Result'+ '\t'+ Transplant_Result+'\n')
            result.write('PGS'+ '\t'+ 'NA'+'\n')
            result.close()
            fileout.write(chip_id+'\n')
        fileout.close()
        #生成QC PDF报告
        command = "python3 %s/%s %s %s %s %s %s %s"%(script,'UpdateMySQLdb.py',outdir,report_file+'.chipid',user_name,sheet,merge1_file,chrX_merge1_file)
        os.system(command)
    except Exception as e:
        if index <= 5:
            time.sleep(10)
            PGH_QC_report(json_file,report_file,sheet,user_name,merge1_file,chrX_merge1_file,cnv_merge_report_file,outdir,log_file,script,index)
        else:
            logger.error("%s : Failed PGH_QC Report" % (index))
            logger.error(e.args)
            logger.error('=====================================================')
            logger.error(traceback.format_exc())
            logger.error('=====================================================')
            exit(1)

def get_ADO_ADI_female(report_file,sheet,cnv_merge_report_file,outdir,log_file):
    df_report = pd.read_table(report_file,sep='\t',header=9,low_memory=False)
    sample_set = set()
    switch1 = 0
    switch2 = 0
    mother_id = ''
    father_id = ''
    embryo = set()
    hash_chr = {}
    hash_result = {}
    #读取CNV文件
    file = open(cnv_merge_report_file, 'r')
    for line in file.readlines():
        line = line.strip()
        data = line.split('\t')
        if data[0] not in hash_chr.keys():
            hash_chr[data[0]] = set()
        hash_chr[data[0]].add(data[1])
    file.close()
    #读取samplesheet文件
    file = open(sheet,'r')
    for line in file.readlines():
        line = line.strip()
        data = line.split('\t')
        sample_set.add(data[0])
        if data[0] not in hash_result.keys():
            hash_result[data[0]] = {}
            hash_result[data[0]]['gender'] = 'unknown'
            hash_result[data[0]]['BAF_SD'] = 'nan'
            hash_result[data[0]]['BAF_MAPD'] = 'nan'
            hash_result[data[0]]['ado'] = 'NA'
            hash_result[data[0]]['adi'] = 'NA'
        if data[1] == 'mother':
            mother_id = data[0]
            switch1 = 1
        elif data[1] == 'father':
            father_id = data[0]
            switch2 = 1
        if data[2] == 'embryo':
            embryo.add(data[0])
    file.close()
    temp_data1 = []
    for x in range(1,23):
        m = str(x)
        temp_data1.append(m)
    if mother_id != '' and father_id != '':
        #mother and father allele
        df_mother = pd.DataFrame(df_report[(df_report['Sample ID'] == mother_id) & (df_report['Allele1 - Top'] != '-') & (df_report['Allele2 - Top'] !=
     '-') & (df_report['Allele1 - Top'] == df_report['Allele2 - Top'])],columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
        df_mother.columns= ['SNP Name','Chr','Position','mother_Allele1','mother_Allele2']
        df_father = pd.DataFrame(df_report[(df_report['Sample ID'] == father_id)& (df_report['Allele1 - Top'] != '-') & (df_report['Allele2 - Top'] !=
     '-') & (df_report['Allele1 - Top'] == df_report['Allele2 - Top'])],columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
        df_father.columns= ['SNP Name','Chr','Position','father_Allele1','father_Allele2']
        df_parent = pd.merge(df_father,df_mother,on =['SNP Name','Chr','Position'],how = 'inner')
        for kk in sample_set:
            df1 = df_report[df_report['Sample ID'] == kk]
            df2 = df1[(df1['Allele1 - AB']!='-') & (df1['Allele2 - AB']!='-')]
            ALL_num = df1['SNP Name'].count()
            if ALL_num == 0:
                sys.exit()
            Ynum = df2[df2['Chr'] == 'Y']['SNP Name'].count()
            hetx = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] != df2['Allele2 - AB'])]['SNP Name'].count()
            hmox = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] == df2['Allele2 - AB'])]['SNP Name'].count()
            key_gender = gender(hetx,hmox,Ynum)
            hash_result[kk]['gender'] = key_gender
            df_autosome = df2[(df2['Chr'] !='0') & (df2['Chr'] !='X') & (df2['Chr'] !='Y') & (df2['Chr'] !='XY') & (df2['Chr'] !='MT')]
         #计算杂合位点数目
            temp=[]
            for key in temp_data1:
                df_a1 = df_autosome[(df_autosome['Chr'] == key) & (df_autosome['Allele2 - AB'] != df_autosome['Allele1 - AB'])]
                count = df_a1['SNP Name'].count()
                temp2 = str(key)+':'+str(count)
                temp.append(temp2)
            value1 = '{'+kk+':{'+','.join(temp)+'}}'
            hash_result[kk]['heteo_num'] = value1
            #计算BAF离散程度以及ADO/ADI
            #将CNV区域的allele 改为-
            df_filter = df_autosome
            if kk in hash_chr.keys():
                region = hash_chr[kk]
                for key in region:
                    if not re.findall('-',key):
                        df_filter = df_filter[df_filter['Chr'] != key]
                    else:
                        temp1 = key.split(':')
                        chr1 = temp1[0]
                        temp2 = temp1[1].split('-')
                        start1 = int(temp2[0])
                        end1 = int(temp2[1])
                        df_filter.loc[((df_filter.Chr == chr1) & (df_filter.Position >= start1) & (df_filter.Position <= end1)),'Allele1 - AB'] = '-'
            #计算BAF
            df_het = pd.DataFrame(df_filter[(df_filter['Allele1 - AB'] != '-') & (df_filter['Allele2 - AB'] != '-') & (df_filter['Allele1 - AB'] != df_filter['Allele2 - AB'])],columns = ['SNP Name','Chr','Position','B Allele Freq'])
            df_BAF = df_het.dropna()
            #计算MAPD
            temp_MAF = 0
            MAPD_data = []
            df_BAF.rename(columns={'B Allele Freq':'BAF'},inplace=True) 
            for key in df_BAF.itertuples():
                BAF_Value = getattr(key,'BAF')
                BAF1 = float(BAF_Value)-float(temp_MAF)
                if BAF1 < 0:
                    BAF1 = -BAF1
                MAPD_data.append(BAF1)
                temp_MAF = BAF_Value
            BAF_SD = df_BAF['BAF'].std()
            BAF_SD = '%.4f'%(BAF_SD)
            MAPD = np.median(MAPD_data)
            BAF_MAPD = '%.4f'%(MAPD)
            hash_result[kk]['BAF_SD'] = BAF_SD
            hash_result[kk]['BAF_MAPD'] = BAF_MAPD
            #计算ado adi
            if kk not in embryo:
                hash_result[kk]['ado'] = 'NA'
                hash_result[kk]['adi'] = 'NA'
                continue
            df_embryo = df_filter[(df_filter['Allele1 - AB'] != '-') & (df_filter['Allele2 - AB'] != '-')]
            df_embryo = pd.DataFrame(df_embryo,columns = ['SNP Name','Chr','Position','Allele1 - Top','Allele2 - Top'])
            df_embryo.columns= ['SNP Name','Chr','Position','embryo_Allele1','embryo_Allele2']
            df_ref = pd.merge(df_parent,df_embryo,on =['SNP Name','Chr','Position'],how = 'inner')
            num_ADO_ALL = 0
            num_ADO = 0
            num_ADI_ALL = 0
            num_ADI = 0
            for key in df_ref.itertuples():
                mother_Allele1 = getattr(key,'mother_Allele1')
                mother_Allele2 = getattr(key,'mother_Allele2')
                father_Allele1 = getattr(key,'father_Allele1')
                father_Allele2 = getattr(key,'father_Allele2')
                embryo_Allele1 = getattr(key,'embryo_Allele1')
                embryo_Allele2 = getattr(key,'embryo_Allele2')
                if mother_Allele1 != father_Allele1:
                    num_ADO_ALL +=1
                    if embryo_Allele1 == embryo_Allele2:
                        num_ADO +=1
                if mother_Allele1 == father_Allele1:
                    num_ADI_ALL +=1
                    if embryo_Allele1 != mother_Allele1 or embryo_Allele2 != mother_Allele1:
                        num_ADI +=1
            ADO = 'Nan'
            ADI = 'Nan'
            if num_ADO_ALL != 0:
                ADO = float(num_ADO)/float(num_ADO_ALL)*100
                ADO = '%.2f'%(ADO)
            if num_ADI_ALL != 0:
                ADI = float(num_ADI)/float(num_ADI_ALL)*100
                ADI = '%.2f'%(ADI)
            hash_result[kk]['ado'] = ADO
            hash_result[kk]['adi'] = ADI
    else:
        for kk in sample_set:
            df1 = df_report[df_report['Sample ID'] == kk]
            df2 = df1[(df1['Allele1 - AB']!='-') & (df1['Allele2 - AB']!='-')]
            ALL_num = df1['SNP Name'].count()
            if ALL_num == 0:
                sys.exit()
            Ynum = df2[df2['Chr'] == 'Y']['SNP Name'].count()
            hetx = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] != df2['Allele2 - AB'])]['SNP Name'].count()
            hmox = df2[(df2['Chr'] == 'X') & (df2['Allele1 - AB'] == df2['Allele2 - AB'])]['SNP Name'].count()
            key_gender = gender(hetx,hmox,Ynum)
            hash_result[kk]['gender'] = key_gender
            df_autosome = df2[(df2['Chr'] !='0') & (df2['Chr'] !='X') & (df2['Chr'] !='Y') & (df2['Chr'] !='XY') & (df2['Chr'] !='MT')]
            #计算BAF离散程度以及ADO/ADI
            #将CNV区域的allele 改为-
            df_filter = df_autosome
            if kk in hash_chr.keys():
                region = hash_chr[kk]
                for key in region:
                    if not re.findall('-',key):
                        df_filter = df_filter[df_filter['Chr'] != key]
                    else:
                        temp1 = key.split(':')
                        chr1 = temp1[0]
                        temp2 = temp1[1].split('-')
                        start1 = int(temp2[0])
                        end1 = int(temp2[1])
                        df_filter.loc[((df_filter.Chr == chr1) & (df_filter.Position >= start1) & (df_filter.Position <= end1)),'Allele1 - AB'] = '-'
            #计算BAF
            df_het = pd.DataFrame(df_filter[(df_filter['Allele1 - AB'] != '-') & (df_filter['Allele2 - AB'] != '-') & (df_filter['Allele1 - AB'] != df_filter['Allele2 - AB'])],columns = ['SNP Name','Chr','Position','B Allele Freq'])
            df_BAF = df_het.dropna()
            #计算MAPD
            temp_MAF = 0
            MAPD_data = []
            df_BAF.rename(columns={'B Allele Freq':'BAF'},inplace=True) 
            for key in df_BAF.itertuples():
                BAF_Value = getattr(key,'BAF')
                BAF1 = float(BAF_Value)-float(temp_MAF)
                if BAF1 < 0:
                    BAF1 = -BAF1
                MAPD_data.append(BAF1)
                temp_MAF = BAF_Value
            BAF_SD = df_BAF['BAF'].std()
            BAF_SD = '%.4f'%(BAF_SD)
            MAPD = np.median(MAPD_data)
            BAF_MAPD = '%.4f'%(MAPD)
            hash_result[kk]['BAF_SD'] = BAF_SD
            hash_result[kk]['BAF_MAPD'] = BAF_MAPD
            hash_result[kk]['ado'] = 'Nan'
            hash_result[kk]['adi'] = 'Nan'
    #print(hash_result)
    return hash_result
if __name__ == '__main__':
    parser = ArgumentParser(description = "LIMS Result Post")
    parser.add_argument("--in", dest="input", required=True, help="input file") ###nargs="+"
    parser.add_argument("--report", dest="report_file", required=True, help="report file")
    parser.add_argument("--sheet", dest="sheetfile", required=True, help="samplesheet file")
    parser.add_argument("--user", dest="user_name", required=True, help="user name")
    parser.add_argument("--merge1", dest="merge1_file", required=True, help="merge1 file")
    parser.add_argument("--chrX_merge1", dest="chrX_merge1_file", required=True, help="chrX merge1 file")
    parser.add_argument("--cnv", dest="cnv_file", required=True, help="cnv file")
    parser.add_argument("--out", dest="outdir", required=True, help="out dir")
    parser.add_argument("--script", dest="script_dir", required=True, help="script dir")
    parser.add_argument("--log", dest="log_file", default=None, required=False, help="File to write logging information (optional)")
    args = parser.parse_args()
    json_file = args.input
    report_file = args.report_file
    sheet = args.sheetfile
    user_name = args.user_name
    merge1_file = args.merge1_file
    chrX_merge1_file = args.chrX_merge1_file
    cnv_merge_report_file = args.cnv_file
    outdir = args.outdir
    log_file = args.log_file
    PGH_QC_report(json_file,report_file,sheet,user_name,merge1_file,chrX_merge1_file,cnv_merge_report_file,outdir,log_file,script,0)