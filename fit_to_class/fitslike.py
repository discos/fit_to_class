# -*- coding: utf-8 -*-
"""
Root representation of fits like data set

It's the contatiner for a generic scan data set representation
"""

import logging
import numpy as np
import pdb

from fit_to_class import fitslike_commons

class Fitslike():
    """Base container to represent generic scan data

    It's a semantic representation of scan data fitted in fits like files
    It means important data are grouped by semantic field:
        - Generic observation data
        - Rf inputs
        - Data and their context
        - Map
    """

    def __init__(self, p_representation):
        """
        Istantiate Fitslike object
        Takes an already processed representation from *fits input file.
        This high level rep might be processed again here.
        
        
            Parameters:
                p_representation: dict
                    Generic fits like representation                
                    
        """
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_inputRepr = p_representation.copy()
        if 'summary' in  self.m_inputRepr.keys():
            self.m_isSummary= True
        else:
            self.m_isSummary= False
        self.m_processedRepr={}
        
    def is_summary(self):
        return self.m_isSummary
        
    def get_inputRepr(self):
        """
        Getter input representaion

        Returns
        -------
        None.

        """
        return self.m_inputRepr
    
    def data_channel_integration(self):
        """
        Data and relative coordinates mean channel per channel
        The average is done along one channel over different spectrums
        Bins per bins ( vertically along data matrix )
        
        Warning: It doesn't take into account if data needs to be averaged.
        i.e. it doesn't check if this is a map or not
        It means it modifies:
            - coordinates : [commons to every data table entry]    
                "data_time"       	
    			"data_ra"
    			"data_dec"
    			"data_az"
    			"data_el"
            - spectrum
        
        New key on processed data, same as obove
        "integrated_data"={
            "data_time"       	
    		"data_ra"
    		"data_dec"
    		"data_az"
    		"data_el"
            "spectrum"
        }
        
        New keyword 'integrated_data' is added to inputRepr keyword
        
        Returns
        -------
        A new keyword 'channel_integration' added to m_processedRepr with:
            
        Numpy array with the same bin number, where every bin is averaged
        with the corresponding binf from other spectrum data entries.
        
        Averaged az, el, ra, dec

        """    
        l_processedDictKeys = ["data_time", "data_mjd", "weather", "data_ra", "data_dec",
                               "data_az", "data_el", "spectrum", "data_integration"]         
        for l_feed in self.m_inputRepr.keys(): 
            l_stokes= {}
            for l_ch in self.m_inputRepr[l_feed].keys():                              
                l_chx= self.m_inputRepr[l_feed][l_ch]               
                " group by flag_cal first, integrate on this groups "
                print(str(l_chx))                
                try:
                    l_coord_time= l_chx['coordinates']['data_time'] 
                    l_coord_time_mjd= l_chx['coordinates']['time_mjd']
                    l_weather= l_chx['extras']['weather']
                    l_coord_ra= l_chx['coordinates']['data_ra']
                    l_coord_dec= l_chx['coordinates']['data_dec']
                    l_coord_az= l_chx['coordinates']['data_az']
                    l_coord_el= l_chx['coordinates']['data_el']
                    l_spectrum= l_chx['spectrum']['data']
                    l_integration= len(l_chx['coordinates']['data_time'])
                except KeyError as e:
                    self.m_logger.error("key not found :" + e.args[0])
                    return {}
                l_coord_time= np.mean(l_coord_time, axis= 0)
                l_coord_time_mjd= l_coord_time_mjd[0]
                l_weather= np.mean(l_weather)
                l_coord_ra= np.mean(l_coord_ra, axis= 0)
                l_coord_dec= np.mean(l_coord_dec, axis= 0)
                l_coord_az= np.mean(l_coord_az, axis= 0)
                l_coord_el= np.mean(l_coord_el, axis= 0)
                if l_chx['backend']['data_type']== 'stokes':
                    ""
                elif l_chx['backend']['data_type']== 'spectrum':
                    l_spectrum= np.mean(l_spectrum, axis= 0)                
                l_integration *= l_chx['backend']['integration_time']
                l_intDataList= [l_coord_time, l_coord_time_mjd, l_weather, l_coord_ra, l_coord_dec,
                                l_coord_az, l_coord_el, l_spectrum, l_integration]
                l_intDataDict= dict(zip(l_processedDictKeys,l_intDataList))
                " Adding stokes  data dict"
                if l_chx['backend']['data_type'] == 'stokes':
                    l_intDataDict['stokes']= l_stokes
                " Creating integrated data repr. "
                self.m_processedRepr[l_ch]= {}
                self.m_processedRepr[l_ch]['integrated_data']= \
                    l_intDataDict.copy()
                self.m_inputRepr[l_feed][l_ch]['integrated_data']= \
                    l_intDataDict.copy()
                
        return self.m_processedRepr
        
    def dump(self):
        """Object contents dumping"""        
        print(self.m_inputRepr)
        
    def dump_keys(self):
        """Nested keys dumping"""
        
        def _recursive(p_dict):
            """Iteration """
            for l_key, l_value in p_dict.items():
                if type(l_value) is dict:
                    yield from _recursive(l_value)
                else:
                    yield (l_key, l_value)
        
        # iterations
        for l_key, l_value in _recursive(self.m_inputRepr):
            print(l_key)
            
            