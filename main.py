#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
converter entry point
"""

import argparse
import os
import logging
from fit_to_class import fitslike_handler
from fit_to_class import fitslike_commons

output_fname_suffix= '_classfits'
logger= None

def scan(p_path, p_scan_type):
    """Scan input file test
    p_path: string
        subscan folder path
    p_scan_type: string
        subscan type
    """
    global logger
    # Output path building
    l_fh= fitslike_handler.Fitslike_handler('fitszilla', p_scan_type)        
    tail, head = os.path.split(p_path)         
    l_outPath= tail+ head+ output_fname_suffix +'/'        
    l_fh.setOutputPath(l_outPath)
    # Data scan
    l_fh.scan_data(p_path)
    # Data conversion
    l_fh.group_on_off_cal()        
    l_fh.normalize()
    l_fh.ClassFitsAdaptations()        
    # Data to disk
    l_fh.classfitsWrite('raw')
    l_fh.classfitsWrite('on_off')
    l_fh.classfitsWrite('cal')

def run():
    global logger
    # Arg parser
    parser= argparse.ArgumentParser()
    parser.add_argument("folder", help= "Subscan folder")
    parser.add_argument("-t","--type", help= "Subscan type [nod, on_off, otf]",\
                        choices= ['nod', 'on_off', 'otf'], default= 'nod')    
    args= parser.parse_args()    
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
    # Scan and conversion
    l_conversion= scan(args.folder, args.type)
    


