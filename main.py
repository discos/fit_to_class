#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
converter entry point
"""

import argparse
import os
import logging
from fit_to_class import fitslike_commons
import scanoptions
import scanpipeline

def run():
    """
    Running program
    """
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

    # Scan options
    scan_options= scanoptions.ScanOptions(l_logger)

    # Arg parser
    parser= argparse.ArgumentParser()
    # Scan file location
    parser.add_argument("folder", help= "Subscan folder")
    # Type and geometry
    parser.add_argument("-t","--type", help= "Subscan types: " + str(scan_options.available_scan_types) ,\
                        choices= scan_options.available_scan_types, default= 'nod')
    parser.add_argument("-g", "--geometry", help="Define scan geometry rules: " + scan_options.scan_geo_rules , type= str, default= "")    
    # Specific feed 
    parser.add_argument("-fd", "--feed", help="Which feed to be parsed( it'll parse just one feed)", type=int)    
    # Grouping options
    parser.add_argument("-r", "--raw", help= "Avoid grouping data by on off cal, convert all data in raw mode", action= "store_true")
    # Computing oprions
    parser.add_argument("-p", "--parallels", help= "How many input file parsed at a time", type= int)
    # Parsing
    args= parser.parse_args()

    scan_options.folder= args.folder
    scan_options.raw= args.raw
    scan_options.type= args.type
    scan_options.geometry= args.geometry
    scan_options.parallel= args.parallels

    # Scan pipeline
    pipeline= scanpipeline.ScanPipeline(l_logger)
    pipeline.set_scan_options(scan_options)
    pipeline.scan()
    


