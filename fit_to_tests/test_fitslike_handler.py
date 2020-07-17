#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:55:45 2020

@author: debbio

Fitslike Handler objects testing
"""

import pdb
import logging
import os
import shutil
import sys
import argparse

from fit_to_class import fitslike_commons
from fit_to_class import fitslike_handler


def run(p_folder, p_scan_type, p_outPath, p_feed, p_parallelism):
    """Scan input file test"""
    l_fh= fitslike_handler.Fitslike_handler('fitszilla', p_scan_type, p_feed, p_parallelism)        
    tail, head = os.path.split(p_folder)                 
    l_fh.setOutputPath(p_outPath)
    l_fh.scan_data(p_folder)
    l_fh.group_on_off_cal()        
    l_fh.normalize()
    l_fh.ClassFitsAdaptations()        
    #l_fh.classfitsWrite('raw')
    #l_fh.classfitsWrite('on_off')
    #l_fh.classfitsWrite('cal')
    #pdb.set_trace()
                
if __name__ == "__main__":
    # Logger
    l_commons = fitslike_commons.Fitslike_commons()    
    l_logger = logging.getLogger(l_commons.logger_name())
    l_formatter =logging.Formatter('[%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s')
    l_logCh = logging.StreamHandler()            
    l_logCh.setFormatter(l_formatter)
    l_logger.setLevel(logging.INFO)
    l_logger.info("Testing fitslike handler unit")
    if not len(l_logger.handlers):
        l_logger.addHandler(l_logCh)
    # Parser
    parser= argparse.ArgumentParser()
    parser.add_argument("-f", "--feed", help= "Work only this feed or all", type= int, default= None)
    parser.add_argument("-i", "--folder", help= "Scan input folder to test basic fit_to_class class", type= str)
    parser.add_argument("-t", "--type", help= "Scan type", choices= ['on_off', 'nodding', 'map'])    
    parser.add_argument("-o", "--output", help= "Output folder to test basic fit_to_class class", type= str, default='')
    parser.add_argument("-p", "--parallelism", help= "Number of input file to be processed at once", type= int, default=1)
    
    args= parser.parse_args()
    # Input File
    if not args.folder:
        l_logger.error("Missing input folder argument")
        sys.exit(0)
    if not os.path.exists(args.folder):
        l_logger.error("Input folder not exists")
        sys.exit(0)
    # Output dir
    # It has to be different from input folder!!
    # eitherway you'll erase input data!
    l_outpath=''
    if args.output == args.folder:
        l_outpath= os.path.join(l_outpath, 'fit_to_class')
    elif args.output == '':
        l_outpath= os.path.join(args.folder, 'fit_to_class')
    # Cleanign destination
    if os.path.exists(l_outpath):
        shutil.rmtree(l_outpath)
    os.mkdir(l_outpath)
    # Run test
    l_fhu= run(args.folder,args.type, l_outpath, args.feed, args.parallelism)    
        


