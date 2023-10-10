#! /usr/bin/python3

import os
import sys
import pprint
from argparse import ArgumentParser
import traceback
import configparser

VERSION = "1.0.0"

selffile = os.path.abspath(__file__)
selfdir = os.path.dirname(selffile)
config_path = selfdir + '/config.ini'
api_dir = selfdir + '/../PGH_Bin'
database_dir = selfdir + '/../PGH_Database'
scripts_dir = selfdir + '/../PGH_Scripts'

class myconf(configparser.ConfigParser):
    def __init__(self,defaults=None):
        configparser.ConfigParser.__init__(self,defaults=None)
    def optionxform(self, optionstr):
        return optionstr

def config(section_name, option_name):
    config_reader = myconf()
    config_reader.read(config_path)

    result = ''
    
    if section_name in config_reader.sections():
        if option_name in config_reader.options(section_name):
            result = config_reader.get(section_name, option_name)
            if section_name == 'database':
                result = os.path.join(database_dir, result)
        
    if section_name == 'path' and option_name == 'api':
        result = api_dir
    if section_name == 'path' and option_name == 'script':
        result = scripts_dir
    return result

def configitems(section_name):
    config_reader = myconf()
    config_reader.read(config_path)
    return config_reader.items(section_name)

def main():
    parser = ArgumentParser(description = "Config")
    #parser.add_argument("--in", dest="input", required=True, help="input file") ###nargs="+"
    #parser.add_argument("--out", dest="output", required=True, help="output file")
    #parser.add_argument("--log", dest="log_file", default=None, required=False, help="File to write logging information (optional)")
    args = parser.parse_args()

    #a = AutoVivification()
    #a[1][2][3] = 4

    print(configitems('database'))






if __name__ == "__main__":
    main()
