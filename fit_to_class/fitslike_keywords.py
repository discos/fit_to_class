#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:26:18 2020

@author: debbio
"""
import json
import os
import logging
import pdb

from fit_to_class import fitslike_commons

class Keyword_json():
    """Json keyword reader

    It has to instantiated for every input files
    eg fitszilla, mbfits..
    
    TODO Only keyword now
    """

    def __init__(self, p_fitsType):
        """
        Open and decodes json file for the selected type
        and the json fitszilla representation

        Parameters
        ----------
        p_fitsType : fits type
            'fitszilla',
            ...

        Returns
        -------
        None

        """
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_jsonDefs = None
        self.m_components = None
        # Input fits type, reading json parsing definition
        l_jsonFileName = p_fitsType + '.json'       
        l_jsonPath= os.path.join(os.path.dirname(__file__), l_jsonFileName)
        if os.path.exists(l_jsonPath):
            with open(l_jsonPath, 'r') as l_jsonFile:
                try:
                    # Json definition loaded
                    self.m_jsonDefs = json.load(l_jsonFile)
                except IOError:
                    self.m_logger.error('Json file definition not found for '
                                  + p_fitsType)
                    return  
        else:
            self.m_logger.error("{} path not found!".format(l_jsonPath))
            return        
        # Decoding json definitions
        self.m_components = self.m_jsonDefs['components'].keys()

    def fitslike_components(self):
        """
        Getter fitslike components keys

        Parameters
        ----------
        l_component : string
            Component name

        Returns
        -------
        List of fitslike ccomponent list

        """
        return self.m_components

    def inputs_components(self):
        """Returns input file semantic components"""
        return self.m_components

    def input_tables_dict(self):
        """
        Returns input file table name list

        Returns
        -------
        Input file table name list

        """
        if 'tables' in self.m_jsonDefs.keys():
            return self.m_jsonDefs['tables']
        return {}

    def parser(self, p_component):
        """
        Returns appropriate parser for p_component
        Parser binds input file keywords to fitslike keywords

        Parameters
        ----------
         p_component : p_component
            Choosen component
            eg. 'observation', 'map', ..

        Returns
        -------
        Parser dictionary

        """
        if p_component in self.m_jsonDefs['components'].keys():
            return self.m_jsonDefs['components'][p_component]['parser']
        return []


class Keyword_helper():
    """Helper keyword handling

    It works only with given json format
    """

    @staticmethod
    def decode_entries(p_entry):
        """Return table for given parse dictionary and given key"""        
        if len(p_entry) < 2:
            return None, None, None
        return p_entry[0], p_entry[1], p_entry[2]
    

    @staticmethod
    def get_table_index(p_tableDict, p_tableName):
        """Decode table name into table index"""
        if p_tableName not in p_tableDict:
            return None
        return p_tableDict[p_tableName]

    @staticmethod
    def is_header_or_data(p_tableName):
        """Decode table name, header or data ?
           Returns True for 'HEADER', false for 'DATA'
        """
        if p_tableName == 'HEADER':
            return True
        else:
            return False
