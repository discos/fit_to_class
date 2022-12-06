#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
converter entry point
"""

import argparse
from fit_to_class.log_formatter import LogFormatter
import logging
import warnings
from astropy.utils.exceptions import AstropyWarning
from fit_to_class import fitslike_commons
from fit_to_class import scanoptions
from fit_to_class import scanpipeline

def run():
    """
    Running program
    """    
    # Soppressione astropy warning riguardo ai frame usati nella trasformazione di coordinate
    warnings.simplefilter('ignore', AstropyWarning)
    # Logger
    l_commons = fitslike_commons.Fitslike_commons()    
    l_logger = logging.getLogger(l_commons.logger_name())    
    l_logCh = logging.StreamHandler()            
    l_logCh.setFormatter(LogFormatter())
    l_logger.setLevel(logging.WARNING)
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
                        choices= scan_options.available_scan_types)
    parser.add_argument("-g", "--geometry", help="Define scan geometry rules: " + scan_options.scan_geo_rules , type= str, default= "")    
    # Specific feed 
    parser.add_argument("-fd", "--feed", help="Which feed to be parsed( it'll parse just one feed)", type=int, default= None)    
    # Grouping options
    parser.add_argument("-r", "--raw", help= "Avoid grouping data by on off cal, convert all data in raw mode", action= "store_true")
    # Computing oprions
    parser.add_argument("-p", "--parallels", help= "How many input file parsed at a time", type= int, default= 4)
    # Pipeline conf file
    parser.add_argument("-s", "--scanconf", help="Pipeline json conf file path" , type= str, default= "scan_conf/conf.json")    
    # Parsing
    args= parser.parse_args()
    
    # Build scan options
    try:
        scan_options.folder= args.folder
        scan_options.raw= args.raw
        scan_options.type= args.type
        scan_options.geometry= args.geometry
        scan_options.parallel= args.parallels
        scan_options.scan_conf_path= args.scanconf
    except ValueError as e:
        l_logger.error("Got wrong input parameter :\n" + str(e))

    l_logger.info("\n" + str(scan_options))
    
    # Scan pipeline
    pipeline= scanpipeline.ScanPipeline(l_logger)
    pipeline.set_scan_options(scan_options)
    pipeline.pipeline_start()
    


