# -*- coding: utf-8 -*-

import os
import re
import logging
from astropy.io import fits
import numpy as np
import astropy.units as unit
from astropy.table import QTable, vstack, Column
from astropy.time import Time
import astropy.constants as const
from multiprocessing import Pool
import sys
import os
import shutil
import pdb

from fit_to_class import fitslike
from fit_to_class  import awarness_fitszilla
from fit_to_class import fitslike_commons
from fit_to_class import fitslike_scangeometry


def _async_subscan(p_logger, p_ftype, p_feed, p_inpath, p_outpath) -> dict:
    """
    Embeds Fitslike subscan parsing

    Parameters
    ----------

    p_logger: obj
        Logger di sistema
    p_ftype : string
        fits file type.
    p_path : string
        path fits file.

    Returns
    -------
    ....

    """
    "todo comporre il parsing attraverso il fitslike ?"
    p_logger.info("Scanning " + p_inpath)
    l_path, l_filename= os.path.split(p_inpath)
    p_logger.info("fitszilla parsing: " + p_inpath)
    try:
        l_fits = fits.open(p_inpath)
    except Exception as e:
        p_logger.error("Skipping input file {} due to an excpetion: \n{}".format(p_inpath, e))
        return {}
    l_aware= awarness_fitszilla.Awarness_fitszilla(l_fits, p_inpath, l_filename , p_feed)
    l_aware.setOutputPath(p_outpath)
    l_aware.parse()
    l_repr= l_aware.process()
    l_errors= l_aware.getErrorList()
    l_fits.close()
    l_fitslike= fitslike.Fitslike(l_repr)
    l_aware = None
    l_fits =None
    l_repr= l_fitslike.get_inputRepr()
    l_repr['file_name']= l_filename
    l_repr['errors']= l_errors
    for er in l_errors:
        p_logger.error(er)
    return l_repr

class Fitslike_handler():
    """
    Subscan master operation handler.
    Based on input data type starts appropiate parsing/processing

    Reading
        Takes a path as input directory with source'subscans
        Spawns multiprocessing parsing of input subscan files
        In particular add a process to Pool for every fits input file through
        global(at module level) functions:
            _envelope_subscan
            _subscan_callback
        This approch comes in handy to overcome  Pickle serialization
        limitations/complexity when dealing with class methods/isntances
        _subscan_callback imposes parsing data are returned as dictionary
        i.e. as python built-in type


    Processing
        It takes fitslike subscan data carrying:
            - data
            - coordinates
            - observation data ( son off cal, site etc..)
        It applies normalization and calibration
        Providing data ready to be studied or to be written to disc (through
        appropriate output conversion module)
    """

    def __init__(self, p_raw, p_scan_geometry, p_scanType, p_feed, p_files_per_scan= 3):
        """
        Internal struct definition

        Parameters
        ----------

        p_raw : bool
            Avoid grouping data by on off cal
        p_scan_geometry:
            List of tuples (N,ON/OFF/CAL)
        p_scanType: string
            Input scan data type
            Nod on off map..
        p_feed: string
            Feed to be processed, None all feeds
        p_files_per_scan: int
            Input file to be processed per pool

        Returns
        -------
        None.

        """
        self._critical_error= False
        self.m_raw= p_raw
        self.m_geometry= p_scan_geometry
        self.m_geometry_scan_list= []
        self._geo_group_on_off_cal= []
        self.m_files_per_pool= p_files_per_scan
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_inputType = 'fitszilla'
        self.m_scanType= p_scanType
        self.m_subscans=[]
        self.m_dataDir =""
        self.m_feed= p_feed
        self.m_group_on_off_cal={}
        self.m_summary= {}
        self.m_outputPath= ''
        self.m_no_cal= False

    def setOutputPath(self, p_path):
        """ output path setter"""
        self.m_outputPath= p_path
        # Output dir creation
        try:
            if os.path.exists(self.m_outputPath):
                shutil.rmtree(self.m_outputPath)
            os.mkdir(self.m_outputPath)
        except Exception as e:
            self.m_logger.error("Error creating tmp dir {}".format(e))

    def get_summary(self):
        """Getter summary"""
        return self.m_summary

    def has_critical_error(self) -> bool:
        return self._critical_error

    def scan_data(self, p_dataDir):
        """
        Takes data input directory and lunchs subscan data conversion
        to fitslike (intermediate rapresentation)

        Parameters
        ----------
        p_dataDir : string
            Input data directory path

        Returns
        -------
        None.

        """
        global g_subscans

        def chunker(p_seq, p_size):
            return (p_seq[pos:pos + p_size] for pos in range(0,len(p_seq), p_size))

        # Parsing dir
        self.m_dataDir = p_dataDir
        if not os.path.isdir(p_dataDir):
            self.m_logger.error("Input data dir is not valid")
            return
        # Processing
        l_filt= self.m_commons.filter_input_files(self.m_inputType)
        l_inputFiles= []
        for item in os.listdir(self.m_dataDir):
            if os.path.isfile(os.path.join(self.m_dataDir,item)):
                l_inputFiles.append(item)
        l_inputFiles= [f for f in l_inputFiles if re.match(l_filt,f)]
        # Split parsing multiprocessing
        self.m_results=[]
        l_poolSize= self.m_files_per_pool
        for l_group in chunker(l_inputFiles, l_poolSize ):
            l_results=[]
            self.m_pool=Pool(l_poolSize)
            for l_fPath in l_group:
                l_results.append(self.m_pool.apply_async(
                        _async_subscan,
                        [self.m_logger,
                        'fitszilla',
                        self.m_feed,
                        p_dataDir +'/'+ l_fPath,
                        self.m_outputPath
                        ])
                    )
            self.m_pool.close()
            self.m_pool.join()
            self.m_subscans= self.m_subscans + [x.get() for x in l_results]
        self.m_logger.info("subscan numbers " + str(len(self.m_subscans)))

    def geometry_group(self):
        """
        Group scans by supposed geometry
        First group is made by file sequence
        Geometry has fixed number of file for every 
        """
        l_loop_len= fitslike_scangeometry.get_geo_loop_len(self.m_geometry)
        # Subscan filename sorting without the summary file
        self.m_subscans_geo= sorted( ( f for f in self.m_subscans if 'summary' not in f['file_name']),\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))
        # Ceck if number of scans is multiple of geo file len
        if len(self.m_subscans_geo) % l_loop_len != 0:
            self.m_logger.error("we have {} subscans and a geometry loop of {}, MCM criteria not satisfied!".format(len(self.m_subscans_geo), l_loop_len))
            self.critical_error= True
            return
        # Group sorted file by geo loop len        
        num_geo_groups= len(self.m_subscans_geo) / l_loop_len
        self.m_logger.info("Grouping subscans in {} repetitions".format(num_geo_groups))
        self.m_geometry_scan_list= [s.tolist() for s in np.array_split(self.m_subscans_geo, l_loop_len) ]
        for gr in self.m_geometry_scan_list:
            self.m_logger.debug("Geometry subscan group :\n{}".format(gr))

    def geometry_check(self):
        """
        Step to verify grouped subscans if they fit to given geometry
        """
        for gr in self.m_geometry_scan_list:
            group_geo_on_off_cal= self.group_on_off_cal(gr)
            check= self._geometry_check_group(group_geo_on_off_cal)
            if check:
                self._geo_group_on_off_cal.append(group_geo_on_off_cal)
        # Debug        
        self.m_logger.debug("Geometry groups by loop lenght:\n {}".format(self._geo_group_on_off_cal))

    def _geometry_check_group(self, p_group_dict) ->bool:
        """
        Given a grouped and worked subscan dictionary by on off cal,
        a check to match geometry is made

        p_group_dict: dict grouper per feed in a single geo loop
        
        Returns:
            True check ok
        """
        

    def group_on_off_cal(self, p_group_list= None):
        """
        Creates a new dict for :
            - On
            - Off
            - Cal_on
            - Cal_off
            For every ch_x

        it means:
            feed{
             ch_x{
                on:[]
                off:[]
                cal_on:[]
                cal_off:[]
                }
             }

        g_on_off_cal_dict will store data groups per type
        Subscan collection parsing per feed
        Output groups are separeted basing on scan type,
        but generally per feed
        """
        subscan_to_work= None
        if p_group_list:
            subscan_to_work= p_group_list
        else:
            # Subscan filename sorting
            self.m_subscans= sorted(self.m_subscans,\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))
            subscan_to_work= self.m_subscans
                
                
        # Recognition and grouping for ON,OFF,CAL
        for subscan in subscan_to_work:
            l_file_name= ''
            if 'file_name' in subscan:
                l_file_name= subscan['file_name']
                self.m_logger.info(l_file_name)
            # Summary file doesn't group
            if 'summary' in subscan.keys():
                self.m_summary= subscan
                self.m_logger.info("summary.fits excluded from on off cal")
                continue
            # Analyzing every subscan
            for l_feed in subscan:
                # Work only selected feed
                if self.m_feed:
                    if not self.m_feed in str(l_feed):
                        continue
                # Prepare on off group keys hirearchy
                if l_feed not in self.m_group_on_off_cal.keys():
                    self.m_group_on_off_cal[l_feed]={}
                # Section (l_chx) navigation
                for l_chx in subscan[l_feed]:
                    if 'ch_' not in l_chx:
                        continue
                    l_chObj= subscan[l_feed][l_chx]
                    # on off group keys hierarchy
                    if l_chx not in self.m_group_on_off_cal[l_feed].keys():
                        # replicate section data to group on off cal to retrieve them easely
                        self.m_group_on_off_cal[l_feed][l_chx]= l_chObj
                    l_feed= l_chObj['frontend']['feed']
                    if 'pol_tables_dict' not in l_chObj.keys():
                        continue
                    l_pol_tables_dict= l_chObj['pol_tables_dict']
                    if 'pols' not in self.m_group_on_off_cal[l_feed][l_chx].keys():
                        self.m_group_on_off_cal[l_feed][l_chx]['pols']={}
                    l_pols_dict= self.m_group_on_off_cal[l_feed][l_chx]['pols']
                    for l_pol in l_pol_tables_dict.keys():
                        if l_pol not in l_pols_dict.keys():
                            l_pols_dict[l_pol]= {'signal': [],'reference': [],'cal_on':[],'cal_off': []}
                        try:
                            l_pol_table= QTable.read(l_pol_tables_dict[l_pol], memmap= True)
                            # ID data type by columns, creating a list made by
                            if 'cal_on' in l_pol_table.colnames:
                                l_pols_dict[l_pol]['cal_on'].append(l_pol_tables_dict[l_pol])
                            elif 'cal_off' in l_pol_table.colnames:
                                l_pols_dict[l_pol]['cal_off'].append(l_pol_tables_dict[l_pol])
                            elif 'signal' in l_pol_table.colnames:
                                l_pols_dict[l_pol]['signal'].append(l_pol_tables_dict[l_pol])
                            elif 'reference' in l_pol_table.colnames:
                                l_pols_dict[l_pol]['reference'].append(l_pol_tables_dict[l_pol])
                        except Exception as e:
                            self.m_logger.error("{}-{}-{}-reading disk table exception : {}".format(l_feed, l_chx, l_pol,str(e)))

        # Feed base work dir cleaning
        self.m_group_path= os.path.join(self.m_outputPath, 'fits_groups')
        if os.path.exists(self.m_group_path):
            shutil.rmtree(self.m_group_path)
        os.mkdir(self.m_group_path)
        # Join tables from on off group
        for l_feed in self.m_group_on_off_cal.keys():
            for l_section in self.m_group_on_off_cal[l_feed].keys():
                if 'ch_' not in l_section:
                    continue
                 # Work dir creation per feed
                l_this_group_path= os.path.join(self.m_group_path, 'group_feed_{}'.format(l_feed))
                if not os.path.exists(l_this_group_path):
                    os.mkdir(l_this_group_path)
                for l_pol in self.m_group_on_off_cal[l_feed][l_section]['pols'].keys():
                    # Join by pol and cal signal
                    for l_type in self.m_group_on_off_cal[l_feed][l_section]['pols'][l_pol].keys():
                        try:
                            l_table_files= self.m_group_on_off_cal[l_feed][l_section]['pols'][l_pol][l_type]
                            if not l_table_files:
                                continue
                            if len(l_table_files)==0:
                                continue
                            l_opened_tables= [QTable.read(t, memmap= True) for t in l_table_files]
                            joined= vstack(l_opened_tables)
                            l_fname= "{}_{}_{}_{}.fits".format(l_feed, l_section, l_pol, l_type)
                            l_gr_table_path= os.path.join(l_this_group_path, l_fname)
                            self.m_group_on_off_cal[l_feed][l_section]['pols'][l_pol][l_type]= l_gr_table_path
                            joined.write(l_gr_table_path)
                            # Deleting input tables
                            # for input_file in l_table_files:
                            #     try:
                            #         os.remove(input_file)
                            #     except IOError as e:
                            #         self.m_logger.error("Can't delete intermediate file {} : {}".format(input_file, e))
                            self.m_logger.info("{}_{}_{}_{}".format(l_feed, l_section, l_pol, l_type))
                            self.m_logger.info("Grouping table into {}".format(l_gr_table_path))
                        except Exception as e:
                            self.m_logger.error("{}-{}-{}-{} error writing group table on disk : {}".format(l_feed, l_section, l_pol, l_type, str(e)))
        # Deleting intermediate folder
                

    def normalize(self):
        """
        On - Off - Cal normalize

        cal on off, off , on are grouped:
            nod, on off:
                subscan - feed - pol -  on off calon caloff  [astropy table group]
            map : ?

        In general terms calculations resemble the following:

            off, average
            cal_on, average
            cal_off, average

            (on - offAvg) / offAvg

            (calOn - on[0])/ on[0] [it should be onRef]
            (calOff - off[0])/ off[0] [it should be onRef]

            Check calOn calOff Normalized , nan, inf

            if [calOnNorm, calOffNorm] is a long array. average / median the array

            calibrationFactor= 1 / [calOnNorm, calOffNorm] * cal_temp

            normalized= on * calibrationFactor

        If cal_on/cal_off are missing or off are missing simply it skips calculations,
        and provides only signal data table to:
            - output_path/fits_normalized
        Data output will be grouped by feed_sectiìon_pol

        """

        def _tableToDict(p_table) -> dict:
            """
            Disgregate :) and build a single row dictionary
            Gets unit from column names _u_unit
            """
            "  "
            l_out= {}
            if len(p_table) == 0:
                self.m_logger.error("Empty table..why?")
                return l_out
            for field in p_table.colnames:
                splits= field.split('_u_')
                if len(splits) > 1:
                    l_out[splits[0]]= p_table[field].data.data
                    l_out[splits[0]] *= unit.Unit(splits[1])
                else:
                    try:
                        l_out[field]= p_table[field].data.data
                    except TypeError:
                        # Cannot average strings
                        l_out[field]= p_table[field][0]
            return l_out


        # Norm output dir creation
        self.m_norm_path= os.path.join(self.m_outputPath, 'fits_normalized')
        if os.path.exists(self.m_norm_path):
            shutil.rmtree(self.m_norm_path)
        os.mkdir(self.m_norm_path)
        # Feed traverse
        for l_feed in self.m_group_on_off_cal.keys():
            # Work only selected feed, or All
            if self.m_feed:
                if self.m_feed not in str(l_feed):
                    continue
            
            if not str(l_feed).isdigit():
                continue
            # Feed folder
            l_feed_path= os.path.join(self.m_norm_path, 'normalized_feed_{}'.format(l_feed))
            if not os.path.exists(l_feed_path):
                os.mkdir(l_feed_path)
            # Section traverse
            for ch in self.m_group_on_off_cal[l_feed].keys():
                l_section= self.m_group_on_off_cal[l_feed][ch]
                l_section_pols= self.m_group_on_off_cal[l_feed][ch]['pols']
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    l_calMarkTemp= None
                    if pol not in l_section_pols.keys(): continue
                    # Get calibration mark temperature if available
                    try:
                        l_calMarkTemp= l_section['frontend']['cal_mark_temp']
                    except:
                        self.m_logger.warning("[{}][{}][{}] no cal mark temp available".format(l_feed, ch, pol))

                    # Here i have feed - section - pol - on/off/cals
                    # Read tables on off cal_on cal_off
                    l_opened_tables={
                        'signal': None,
                        'reference' : None,
                        'cal_on': None,
                        'cal_off': None
                        }
                    l_data={
                        'signal': QTable(),
                        'reference' : QTable(),
                        'cal_on': QTable(),
                        'cal_off': QTable()
                        }                    
                    for on_off_cal in l_opened_tables.keys():
                        if on_off_cal not in l_section_pols[pol].keys():
                            continue
                        try:                            
                            l_table_path= l_section_pols[pol][on_off_cal]
                            if not len(l_table_path):
                                continue                            
                            l_opened_tables[on_off_cal]= QTable.read(l_table_path, memmap= True)
                        except Exception as e:
                            self.m_logger.warning("[{}][{}][{}][{}]".format(l_feed, ch, pol, on_off_cal))
                            self.m_logger.error("Error reading goruped table {} - {}".format(l_table_path, e))

                    # Now i've loaded all tables needed to calibrate one section
                    # Processing reference mean
                    try:                        
                        # Signal whole data set
                        l_data['signal']= l_opened_tables['signal']['data']
                        if not l_opened_tables['signal']:                            
                            continue
                        if not len(l_opened_tables['signal']):
                            continue                        
                        if l_opened_tables['reference'] != None:
                            if len(l_opened_tables['reference']) != 0:
                                l_opened_tables['reference']= l_opened_tables['reference'].group_by(['pol']).groups.aggregate(np.mean)
                                l_data['reference']= l_opened_tables['reference']['data']
                        # Processing cal_on mean norm
                        if l_opened_tables['cal_on'] != None and l_opened_tables['reference'] != None:
                            if len(l_opened_tables['cal_on']) != 0 and len(l_opened_tables['reference']) != 0 :
                                l_opened_tables['cal_on']= l_opened_tables['cal_on'].group_by(['pol']).groups.aggregate(np.mean)
                                l_data['cal_on']= l_opened_tables['cal_on']['data']
                                l_data['cal_on']= (l_data['cal_on'] - l_data['signal'][0]) / l_data['signal'][0]
                        # Processing cal_off mean norm
                        if l_opened_tables['cal_off'] != None and l_opened_tables['reference'] != None:
                            if len(l_opened_tables['cal_off']) != 0 and len(l_opened_tables['reference']) != 0 :
                                l_opened_tables['cal_off']= l_opened_tables['cal_off'].group_by(['pol']).groups.aggregate(np.mean)
                                l_data['cal_off']= l_opened_tables['cal_off']['data']                                            
                                l_data['cal_off']= (l_data['cal_off'] - l_data['reference']) / l_data['reference']
                    except KeyError as e:
                        self.m_logger.warning("[{}][{}][{}]".format(l_feed,ch,pol))
                        self.m_logger.error("Missing mandatory keyword {}".format(e))
                    except Exception as e:
                        self.m_logger.warning("[{}][{}][{}]".format(l_feed,ch,pol))
                        self.m_logger.error("Exception averaging grouped table {}".format(e))

                    # Have we data?
                    l_data_raw= False
                    l_on_off_possible = False
                    l_calibration_possible= False
                    if len(l_data['signal']) != 0:
                        l_data_raw= True
                        self.m_logger.info("[{}][{}][{}] Signal data are present".format(l_feed,ch,pol))
                        # Sub calculations, on - off
                        if len(l_data['reference']) and len(l_data['signal']):                            
                            # Check if calculations are possible
                            # TODO implementare solo conteggi senzi cal
                            l_on_off_possible= True
                            self.m_logger.info("[{}][{}][{}] Reference data are present".format(l_feed,ch,pol))
                            l_calibration_possible= False
                            if l_calMarkTemp:
                                if len(l_data['cal_on']) and len(l_data['cal_off']):
                                    l_calibration_possible= True
                                    self.m_logger.info("[{}][{}][{}] Calibration data are present".format(l_feed,ch,pol))
                    #
                    if not l_data_raw:
                        self.m_logger.warning("[{}][{}][{}] No input data, skipping".format(l_feed,ch,pol))
                        continue
                    self.m_logger.warning("[{}][{}][{}] Applying calibration".format(l_feed,ch,pol))
                    # TODO subs signal table with calibrated data
                    # Calc with units..review
                    # Keep data as small as possible (numpy.float32)
                    try:                        
                        if l_calibration_possible:
                            # Complete calibration
                            # On Off
                            l_signal= l_data['signal']
                            l_reference= l_data['reference']
                            l_data['on_off']= (l_signal -l_reference) / l_reference                            
                            # Calibration factor
                            l_cal = np.concatenate((l_data['cal_on'], l_data['cal_off']))                            
                            if len(l_cal) == 0:
                                continue
                            good =  ~np.isnan(l_cal) & ~np.isinf(l_cal)
                            l_cal = l_cal[good]                            
                            if len(l_cal) > 0:
                                meancal = np.median(l_cal) if len(l_cal) > 30 else np.mean(l_cal)                                
                                calibration_factor = 1 / meancal * l_calMarkTemp             
                                calibration_factor= np.float32(calibration_factor.value)* calibration_factor.unit
                            else:
                                continue                              
                            l_data['calibrated']= l_data['on_off'].astype(np.float32) * calibration_factor                    
                            # Replace data on table                             
                            l_opened_tables['signal']['data']= l_data['calibrated']
                            # Remove useless columns
                            if 'signal' in l_opened_tables['signal'].colnames:
                                del l_opened_tables['signal']['signal']                            
                            # Rewrite table on disk, overwriting signal table
                            l_file_name= "{}_{}_{}_calibrated.fits".format(l_feed, ch, pol)
                            l_file_path= os.path.join(l_feed_path , l_file_name)                                   
                            l_section_pols[pol]['calibrated']= l_file_path
                            l_opened_tables['signal'].write(l_file_path, overwrite= True)                                                        
                        elif l_on_off_possible:
                            # Only counts
                            l_signal= l_data['signal']
                            l_reference= l_data['reference']
                            l_data['on_off']= (l_signal -l_reference) / l_reference                            
                            l_opened_tables['signal']['data']= l_data['on_off']
                            # Remove useless columns
                            if 'signal' in l_opened_tables['signal'].colnames:
                                del l_opened_tables['signal']['signal']                
                            # Writing to disk 
                            l_file_name= "{}_{}_{}_counts.fits".format(l_feed, ch, pol)
                            l_file_path= os.path.join(l_feed_path , l_file_name)                            
                            l_section_pols[pol]['counts']= l_file_path
                            l_opened_tables['signal'].write(l_file_path, overwrite= True)    
                        elif len(l_opened_tables['signal']): 
                            # Signal with no calculations                                                                                
                            # Remove useless columns
                            if 'signal' in l_opened_tables['signal'].colnames:
                                del l_opened_tables['signal']['signal']                
                            # Writing to disk 
                            l_file_name= "{}_{}_{}_signal.fits".format(l_feed, ch, pol)
                            l_file_path= os.path.join(l_feed_path , l_file_name)                            
                            l_section_pols[pol]['signal']= l_file_path
                            l_opened_tables['signal'].write(l_file_path, overwrite= True)    
                        else:
                            # Nothing to do..
                            self.m_no_cal= True
                            self.m_logger.error("[{}][{}][{}] No appropriate data to normalize this dataset".format(l_feed, ch, pol))   
                        # Write data
                    except IOError as e:
                        self.m_no_cal= True                        
                        self.m_logger.error("[{}][{}][{}] Exception applying calibration to this data set : {}".format(l_feed,ch,pol,e))                           
                    except Exception as e:
                        self.m_no_cal= True            
                        self.m_logger.error("[{}][{}][{}] Exception applying calibration to this data set : {}".format(l_feed,ch,pol,e))                   
                    
        # Delete group data
        # if os.path.exists(self.m_group_path):
        #     shutil.rmtree(self.m_group_path, ignore_errors= True)
                        
                    
    def ClassFitsAdaptations(self):
        """
        Raw data groups stored as :
            self.m_group_on_off_cal[l_feed][l_section]['pols'][LL LR RL RR][on off cal_on cal off]

        On disk we have a folder relative to oputput path:
            output_path/fits_groups, folder where table data (disk part) are grouped by feed section pol on off cal
            output_path/fits_norms, folder where grouped data produce a calibrated measures ( if cal is present, counts otherwise )
            output_path/classfits, single, calibrated classfits file boxing every feed and pol

        Conversion tips:
             desired cooord in az, el o ra, dec
            observed coord in crdelt2,3


        """

        def getQTableColWithUnit(p_table, p_col, p_unit):
            """ Get QTable[col_u_unit] data adding unit to returned data """
            l_col_names= p_table.colnames
            for col_name in l_col_names:
                if p_col in col_name:
                    l_field_unit= col_name.split('_u_')
                    if len(l_field_unit) == 2:
                        data= p_table[col_name] * unit.Unit(l_field_unit[1])
                        data= data.to(p_unit).value
                        return data

        # Work folder
        self.m_class_path= os.path.join(self.m_outputPath, 'classfits')
        if os.path.exists(self.m_class_path):
            shutil.rmtree(self.m_class_path)
        os.mkdir(self.m_class_path)
        # Feed traverse
        for l_feed in self.m_group_on_off_cal:
            # Feed filter added by user ?
            if self.m_feed:
                # Check if selected feed is the feed we want to work
                if self.m_feed not in str(l_feed):
                    continue
            # Only feed, avoid extra data on this level
            if not str(l_feed).isdigit():
                continue            
            # Feed folder
            l_feed_path= os.path.join(self.m_class_path, 'class_feed_{}'.format(l_feed))
            if not os.path.exists(l_feed_path):
                os.mkdir(l_feed_path)
            # Sections traverse
            for l_ch in self.m_group_on_off_cal[l_feed]:
                l_section= self.m_group_on_off_cal[l_feed][l_ch]
                l_section_pols= self.m_group_on_off_cal[l_feed][l_ch]['pols']
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    if pol not in l_section_pols.keys(): continue                    
                    # Check work on data on only ( Raw Mode )                    
                    calibration_type= ''
                    l_polarization_table= None                                          
                    # ONLY 1 table per pol!
                    try:                                                             
                        # Data retreiving one level up from 'on' data
                        if 'calibrated' in l_section_pols[pol].keys():
                            calibrated_data_path= l_section_pols[pol]['calibrated']                                
                            self.m_logger.info("Loading normalized {}-{}-{} path {}".format(l_feed, l_ch, pol,calibrated_data_path))
                            if calibrated_data_path: 
                                l_polarization_table= QTable.read(calibrated_data_path, memmap= True)
                                calibration_type= 'calibrated'
                            else:
                                self.m_logger.info("calibrated for {}-{}-{} is empty".format(l_feed, l_ch, pol))                                    
                                
                        if  not l_polarization_table and 'on_off' in l_section_pols[pol].keys():
                            on_off_data_table_path= l_section_pols[pol]['on_off']
                            self.m_logger.info("Loading on_off {}-{}-{} path {}".format(l_feed, l_ch, pol,on_off_data_table_path))
                            if on_off_data_table_path:
                                l_polarization_table= QTable.read(on_off_data_table_path, memmap= True)
                                calibration_type= 'on_off'                                
                            else:
                                self.m_logger.info("on_off data for {}-{}-{} is empty".format(l_feed, l_ch, pol))    
                                
                        if not l_polarization_table and 'signal' in l_section_pols[pol].keys():
                            signal_data_path= l_section_pols[pol]['signal']
                            self.m_logger.info("Loading signal {}-{}-{} path {}".format(l_feed, l_ch, pol,signal_data_path))
                            if signal_data_path:
                                l_polarization_table= QTable.read(signal_data_path, memmap= True)
                                calibration_type= 'signal'
                            else:
                                self.m_logger.info("signal data for {}-{}-{} is empty".format(l_feed, l_ch, pol))                                    
                                
                        if  not l_polarization_table and 'reference' in l_section_pols[pol].keys():
                            reference_data_path= l_section_pols[pol]['reference']
                            self.m_logger.info("Loading reference {}-{}-{} path {}".format(l_feed, l_ch, pol,reference_data_path))
                            if reference_data_path:
                                l_polarization_table= QTable.read(reference_data_path, memmap= True)
                                calibration_type= 'reference'
                            else:
                                self.m_logger.info("reference data for {}-{}-{} is empty".format(l_feed, l_ch, pol))                                                                    
                                
                        # No available input data options..
                        if not l_polarization_table:
                            self.m_logger.warning("{}_{}_{} No data suitable for class conversion".format(l_feed, l_ch, pol))
                            continue
                        
                    except IOError as e:
                        self.m_logger.error("Exception retrieving normalized data from disk tables {}".format(e))
                        continue
                    except Exception as e:
                        self.m_logger.error("Exception retrieving normalized data from disk tables {}".format(e))
                        continue
                    #pdb.set_trace()
                    # Generic observation data copy to dedicated dict, more copies
                    # below during calculations
                    self.m_obs_general_data={}
                    self.m_obs_general_data['ra']= l_section['scheduled']['ra'].to(unit.deg).value
                    self.m_obs_general_data['dec']= l_section['scheduled']['dec'].to(unit.deg).value
                    self.m_obs_general_data['source']= l_section['scheduled']['source']
                    self.m_obs_general_data['date-red']= Time.now().to_datetime().strftime('%d/%m/%y')
                    # classfits new dict filling process
                    l_class_data= {}
                    # ch by ch
                    try:
                        # Data first to get shapes                            
                        l_class_data['SPECTRUM']= l_polarization_table['data']
                        # data shape they must be equals ( on_off cal on )
                        data_shape= l_class_data['SPECTRUM'].shape
                        l_class_data['CRPIX1']=  l_section['backend']['bins'] // 2 + 1
                        # Lavoro con i dati integrati
                        # ut
                        l_tMjd= l_polarization_table['data_time_mjd']
                        l_timeMjd= Time(l_tMjd, format='mjd', scale='utc')                        
                        l_class_data['UT']= ( l_tMjd - np.floor(l_tMjd)) * 86400                        
                        # date
                        l_class_data['DATE-OBS']= l_timeMjd.strftime('%d/%m/%y')
                        # lsts
                        l_lsts= l_timeMjd.sidereal_time('apparent', \
                                  fitslike_commons.Fitslike_commons.\
                                      get_site_location(l_section['scheduled']['antenna']).lon)
                        l_lsts= l_lsts.value * unit.hr
                        # infos
                        l_class_data['OBJECT']= l_section['scheduled']['source']
                        l_class_data['LINE']= "F{}-{:3.3f}-MHz"\
                            .format(l_section['frontend']['feed'], l_section['backend']['bandwidth'])
                        self.m_obs_general_data['line']= l_class_data['LINE']
                        try:
                            #pdb.set_trace()
                            l_class_data['TELESCOP']=\
                                self.m_commons.class_telescope_name(l_section,self.m_summary['summary']['backend_name'], pol, l_ch)
                        except ValueError as e:
                            self.m_logger.error("Missing summary !? TELESCOP set to empty value \n Excpetion: "+ str(e))
                            l_class_data['TELESCOP']= ""
                        l_mH2O= l_polarization_table['weather']
                        l_class_data['MH2O']= l_mH2O
                        # temp
                        l_class_data['TSYS']= 1.0
                        l_class_data['CALTEMP']= l_section['frontend']['cal_mark_temp'].value
                        # time
                        l_class_data['LST'] = l_lsts.to('s').value
                        # CDELT
                        l_class_data['CDELT1']= (l_section['frontend']['bandwidth'] / l_section['backend']['bins']).to('Hz')
                        # freq and velocity
                        #pdb.set_trace()
                        l_class_data['RESTFREQ']= self.m_summary['summary']['restfreq'].to(unit.Hz).value                        
                        # Caso semplificato (<1) se la rest freq è a 0 uso la frequenza del canale bins/2+1
                        if l_class_data['RESTFREQ'] != None:
                            if l_class_data['RESTFREQ'] < 1:
                                l_class_data['RESTFREQ']= l_section['frontend']['local_oscillator'].to('Hz') +\
                                                    l_class_data['CDELT1'] * (l_section['backend']['bins']/2 +1)
                                l_class_data['RESTFREQ']= l_class_data['RESTFREQ'].value
                        self.m_obs_general_data['restfreq']= l_class_data['RESTFREQ']
                        l_class_data['VELOCITY']= l_section['scheduled']['vlsr'].to("m/s").value
                        l_df= (l_section['backend']['bandwidth'] / l_section['backend']['bins']).to('Hz')
                        l_class_data['CDELT1']= l_df.value
                        self.m_obs_general_data['cdelt1']= l_class_data['CDELT1']
                        l_deltav= - l_df/ l_class_data['RESTFREQ'] * const.c
                        l_class_data['DELTAV']= l_deltav.value
                        # LOG test
                        #self.m_logger.warn("RESTFREQ {}".format(l_class_data['RESTFREQ']))
                        # Objects Coordinates   
                        # CRVAL 2,3: target object coordiantes                        
                        l_class_data['CRVAL2']= self.m_summary['summary']['target_ra'].to(unit.deg).value
                        l_class_data['CRVAL3']= self.m_summary['summary']['target_dec'].to(unit.deg).value
                        # CDELT 2,3: observer - target
                        l_observed_ra= getQTableColWithUnit(l_polarization_table,'data_ra', unit.deg)                        
                        l_offset_ra= (l_observed_ra - l_class_data['CRVAL2'])* np.cos(l_class_data['CRVAL3'])
                        l_class_data['CDELT2']= l_offset_ra
                        l_observed_dec= getQTableColWithUnit(l_polarization_table,'data_dec', unit.deg)                        
                        l_offset_dec= (l_observed_dec - l_class_data['CRVAL3'])
                        l_class_data['CDELT3']= l_offset_dec
                        # az el deg
                        l_class_data['AZIMUTH']= getQTableColWithUnit(l_polarization_table,'data_az', unit.deg)
                        l_class_data['ELEVATIO']= getQTableColWithUnit(l_polarization_table,'data_el', unit.deg)
                        # data
                        l_class_data['OBSTIME'] = l_section['backend']['integration_time'].to('s').value
                        l_class_data['MAXIS1'] = l_section['backend']['bins']
                        self.m_obs_general_data['maxis1']= l_class_data['MAXIS1']
                        # we have to shape classfits data properly according to data shape
                        # es data shape is 12, 16384 we have to replicate data this shape
                        rows= data_shape[0]
                        for k in l_class_data.keys():
                            l_value= l_class_data[k]
                            if "SPECTRUM" in k:
                                continue
                            try:
                                # this prevent already shaped (rows,) to be warped
                                l_class_data[k]= np.full((rows,), l_value)
                            except Exception  as e:
                                self.m_logger.error("Error reshaping classfits data: " +str(e))
                        # Transform dictionary based converted table in QTable on disk
                        l_class_table= QTable()
                        for col in l_class_data.keys():
                            try:
                                # Skip empty columns
                                l_class_table.add_column( Column(data= l_class_data[col], name= col) )
                            except ValueError as e:
                                self.m_logger.error("Exception adding column [{}] to classfit converted table {}".format(col, e))
                        # Table name generation                        
                        l_class_name= '{}_{}_{}_{}_class.fits'.format(l_feed, l_ch, pol, calibration_type)       
                        l_class_path= os.path.join(l_feed_path, l_class_name)
                        self.m_logger.info("class table path: {}".format(l_class_path))
                        try:                            
                            self.classfitsWrite(l_class_path, l_class_table)
                        except ValueError as e:                            
                            self.m_logger.error("Error writing small classfit table to disk {}".format(e))   
                    except KeyError as e:
                        self.m_logger.error("[KeyError] Error preparing class data: {}".format(e))                                      
                    except TypeError as e:                        
                        self.m_logger.error("[TypeError] Error preparing class data: {}".format(e))
                    except Exception as e:                        
                        self.m_logger.error("[Exception] Error preparing class data: {}".format(e))
        # Delete norm data
        # try:
        #     if os.path.exists(self.m_norm_path) and self.m_norm_path:
        #         shutil.rmtree(self.m_norm_path, ignore_errors= True)
        # except Exception as e:
        #         self.m_logger.error("Error cleaning normalized data folder: {}\n {}\n".format(self.m_norm_path, e))


    def classfitsWrite(self, p_file_name, p_class_table):
        """
        Convert Qtable classfits column to fits file format adding appropriate
        header data

        Parameters
        ---------
        p_file_name: string
            Output file name
        p_class_table: string
            QTable with clasfits data            
        """              
        l_newCols=[]
        # for every column expressed in classfits definition ..
        for classCol in self.m_commons.getClassfitsColumnsZip():            
            # [ name, form, unit ] column by column data building 
            # fill one column looking into every feed[on], and builds column data 
            l_columnFound = False
            # conditionals, some fields needs dedicated approach            
            # converted fits data matches with classfits columns? 
            # some columns need special care
            l_inferredCol= classCol[0]
            if l_inferredCol in p_class_table.colnames:                
                # found match, add data to column data 
                l_columnFound= True                
            try:
                # adding column to classfits if fitszilla representation matches it
                if l_columnFound:
                    # some fields needs dedicated approach 
                    if classCol[0] == "SPECTRUM":
                        l_rows= p_class_table[l_inferredCol][0].shape[0]
                        l_newCols.append(fits.Column(array= p_class_table[l_inferredCol],name= classCol[0],format= "{}D".format(l_rows),unit= classCol[2]))
                    else:
                        l_newCols.append(fits.Column(array= p_class_table[l_inferredCol],name= classCol[0],format= classCol[1],unit= classCol[2]) )

            except Exception as e:
                self.m_logger.error("classfits column creation exception: {} - {}".format(classCol, e))                

        l_hdData= self.m_obs_general_data
        # header 
        l_hdu= fits.PrimaryHDU()        
        # data 
        try:
            # TEST 
            #for col in l_newCols:
            #    print(col.name + " " + str(col.array.shape))
            l_cdefs= fits.ColDefs(l_newCols)            
            l_hdu= fits.BinTableHDU().from_columns(l_cdefs)
        except Exception as e:            
            self.m_logger.error("Exception creating classfits model file {} - {}".format(p_file_name, e))            
            return
        #header 
        try:
            # TODO Rivedere i campi DELT, PIX ...
            l_hdu.header['EXTNAME'] = "MATRIX"
            l_hdu.header['EXTVER'] = 1
            l_hdu.header['MAXIS'] = 4
            if 'maxis1' in l_hdData.keys():
                l_hdu.header['MAXIS1'] = l_hdData['maxis1']
            else:
                l_hdu.header['MAXIS1'] = 0.0
            l_hdu.header['MAXIS2'] = 4
            l_hdu.header['MAXIS3'] = 4
            l_hdu.header['MAXIS4'] = 4            
            l_hdu.header['CTYPE1']= "FREQ"
            l_hdu.header['CRVAL1']= 0.0
            if 'cdelt1' in l_hdData.keys():
                l_hdu.header['CDELT1'] = l_hdData['cdelt1']
            else:
                l_hdu.header['CDELT1'] = 0.0
            l_hdu.header['CTYPE2']= "RA"
            if 'ra' in l_hdData.keys():
                l_hdu.header['CRVAL2']= l_hdData['ra']
            else:
                l_hdu.header['CRVAL2']= 0.0
            l_hdu.header['CDELT2']= 0.0
            l_hdu.header['CRPIX2']= 0.0
            l_hdu.header['CTYPE3']= "DEC"
            if 'dec' in l_hdData.keys():
                l_hdu.header['CRVAL3']= l_hdData['dec']
            else:
                l_hdu.header['CRVAL3']= 0.0            
            l_hdu.header['CDELT3']= 0.0
            l_hdu.header['CRPIX3']= 0.0
            l_hdu.header['CTYPE4']= "STOKES"
            l_hdu.header['CRVAL4']= 0.0
            l_hdu.header['CDELT4']= 0.0
            l_hdu.header['CRPIX4']= 0.0
            l_hdu.header['SUBSCAN']= 1
            if 'line' in l_hdData.keys():
                l_hdu.header['LINE'] = l_hdData['line']
            else:
                l_hdu.header['LINE'] = " "
            if 'source' in l_hdData.keys():
                l_hdu.header['OBJECT'] = l_hdData['source']            
            else:
                l_hdu.header['OBJECT'] = " "
            if 'restfreq' in l_hdData.keys():
                l_hdu.header['RESTFREQ'] = l_hdData['restfreq']
            else:
                 l_hdu.header['RESTFREQ'] = 0.0
            l_hdu.header['VELDEF']= 'RADI_LSR'
            l_hdu.header['GAINIMAG']= 0.0
            l_hdu.header['BEAMEFF']= 0.0
            l_hdu.header['FORMEFF']= 0.0
            if 'restfreq' in l_hdData.keys():
                l_hdu.header['DATE-RED'] = l_hdData['date-red']        
            else:
                l_hdu.header['DATE-RED'] = " "
            l_hdu.header['EPOCH']= 2000.0
            l_hdu.header['CRVAL'] = 0
        except KeyError as e:
            self.m_logger.error("Exception filling " + p_file_name + " header data: "+ str(e))
        # Disk writing
        try:
            #if os.path.exists(p_file_name):
            #    os.remove(p_file_name)
            l_hdu.writeto(p_file_name)
            self.m_logger.info("Wrote classfits {}".format(p_file_name))
        except Exception as e:            
            self.m_logger.error("Exception writing fitsl  file {} - {}".format(p_file_name, e))
