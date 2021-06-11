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

def scan(p_path, p_raw, p_scan_type, p_feed, p_files_per_scan):
    """Scan input file test
    p_path: string
        subscan folder path
    p_raw: bool
        True implies not grouping bdata by on off cal
    p_scan_type: string
        subscan type
    p_feed: int
        Feed to be parsed, None
    p_files_per_scan: int
        input file parsed at a time
    """
    global logger
    # Output path building
    l_fh= fitslike_handler.Fitslike_handler( p_raw, p_scan_type, p_feed, p_files_per_scan)        
    tail, head = os.path.split(p_path)         
    l_outPath= tail+ head+ output_fname_suffix +'/'        
    l_fh.setOutputPath(l_outPath)
    # Data scan
    l_fh.scan_data(p_path)
    if p_raw:
    # Data conversion
        l_fh.group_on_off_cal()
        l_fh.ClassFitsAdaptations()
    else:
        l_fh.group_on_off_cal()
        l_fh.normalize()
        l_fh.ClassFitsAdaptations()


def run():
    global logger
    # Arg parser
    parser= argparse.ArgumentParser()
    parser.add_argument("folder", help= "Subscan folder")
    parser.add_argument("-t","--type", help= "Subscan type [nod, on_off, otf]",\
                        choices= ['nod', 'on_off', 'otf'], default= 'nod')
    parser.add_argument("-p", "--parallels", help= "How many input file parsed at a time", type= int)
    parser.add_argument("-fd", "--feed", help="Which feed to be parsed( it'll parse just one feed)", type=str)
    parser.add_argument("-r", "--raw", help= "Avoid grouping data by on off cal, convert all data in raw mode", action= "store_true")
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
    scan(args.folder.rstrip('/'), args.raw, args.type, args.feed, args.parallels)
    


