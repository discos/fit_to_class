# -*- coding: utf-8 -*-
"""Fitslike objects testing
"""

from astropy.io import fits
import pdb
import logging
import argparse
import os
import sys

from fit_to_class  import fitslike_commons
from fit_to_class  import fitslike
from fit_to_class  import awarness_fitszilla


def main(p_file, p_type):
    """Parse a fitszilla through Awarness_fitszilla
    p_file: input file path
    p_type: scan type
    """
    l_fits = fits.open(p_file)        
    l_aware = awarness_fitszilla.Awarness_fitszilla(l_fits, p_file)        
    l_intermediateDict = l_aware.parse()
    l_processedDict = l_aware.process()         
    l_fitsLike = fitslike.Fitslike(l_processedDict)        
    pdb.set_trace()

if __name__ == "__main__":
    # Set Logger
    l_commons = fitslike_commons.Fitslike_commons()    
    l_logger = logging.getLogger(l_commons.logger_name())
    l_formatter =logging.Formatter('[%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s')
    l_logCh = logging.StreamHandler()            
    l_logCh.setFormatter(l_formatter)
    l_logger.setLevel(logging.INFO)
    l_logger.info("Testing fitslike unit")
    if not len(l_logger.handlers):
        l_logger.addHandler(l_logCh)
    # Parser
    parser= argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help= "Single file input to test basic fit_to_class class", type= str)
    parser.add_argument("-t", "--type", help= "Scan type", choices= ['on_off', 'nodding', 'map'])
    args= parser.parse_args()
    # Input File
    if not args.file:
        l_logger.error("Missing input file argument")
        sys.exit(0)
    if not os.path.exists(args.file):
        l_logger.error("Input file not exists")
        sys.exit(0)
    # Execution
    main(args.file, args.type)
    
    
    
    
    
    
    
