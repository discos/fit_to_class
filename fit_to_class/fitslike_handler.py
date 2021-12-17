# -*- coding: utf-8 -*-

import os
import re
import logging
import sys
import shutil
import pdb
import traceback
from astropy.io import fits
import numpy as np
import astropy.units as unit
from astropy.table import QTable, vstack, Column
from astropy.time import Time
import astropy.constants as const
from multiprocessing import Pool


from fit_to_class import fitslike
from fit_to_class  import awarness_fitszilla
from fit_to_class import fitslike_commons
from fit_to_class import fitslike_scangeometry
from fit_to_class import fitslike_tables


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
    _path, _filename= os.path.split(p_inpath)
    #p_logger.info("Scanning " + _filename)        
    try:
        _fits = fits.open(p_inpath)
    except Exception as e:        
        p_logger.error(f"Skipping input file {_filename} due to an excpetion: \n{str(e)}")
        return {}
    _aware= awarness_fitszilla.Awarness_fitszilla(_fits, p_inpath, _filename , p_feed)
    _aware.setOutputPath(p_outpath)
    _aware.parse()
    _repr= _aware.process()
    _errors= _aware.getErrorList()
    _fits.close()
    _fitslike= fitslike.Fitslike(_repr)
    _aware = None
    _fits =None
    _repr= _fitslike.get_inputRepr()
    _repr['file_name']= _filename
    _repr['errors']= _errors
    for er in _errors:
        p_logger.error(er)
    return _repr

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
        if not 'summary' in self.m_summary:
            return None
        return self.m_summary['summary']

    def get_geometry_definition(self) -> str:
        """Returns geometry definition as in summary if present"""        
        _summary= self.get_summary()
        if 'geometry' in _summary:
            return _summary['geometry']
        return ""

    def get_scan_type(self) -> str:
        """Returns scan type definition as in summary if present"""
        _summary= self.get_summary()
        if 'scan_type' in _summary:
            return _summary['scan_type']
        return ""

    def has_critical_error(self) -> bool:
        return self._critical_error

    def scan_data(self, p_subscan_list) -> None:
        """
        Takes data input directory and lunchs subscan data conversion
        to fitslike (intermediate rapresentation)

        Parameters
        ----------
        p_subscan_list : dict
            {
                'summary': summary_file
                'subscan': subcan_list(with inner list made after geo grouping)
            }
        """                
        # TODO VERIFCARE CHE SIA CORRETTO
        # Input files
        _summary= p_subscan_list['summary']        
        _subscans= p_subscan_list['subscan']        
        # Split parsing multiprocessing
        self.m_results=[]
        _poolSize= self.m_files_per_pool
        # Split subscan sublist in pools if presents
        _subscans_sublist= [s.tolist() for s in np.array_split(_subscans, _poolSize)]
        for _group in _subscans_sublist:
            _results=[]
            self.m_pool=Pool(_poolSize)
            for _scan in _group:
                _results.append(self.m_pool.apply_async(
                        _async_subscan,
                        [self.m_logger,
                        'fitszilla',
                        self.m_feed,
                        _scan,
                        self.m_outputPath
                        ])
                    )
            self.m_pool.close()
            self.m_pool.join()
            # Adding new data to subscan data list            
            self.m_subscans= self.m_subscans+  [x.get() for x in _results]             
        # Summary parsing on its own                
        self.m_summary= _async_subscan(self.m_logger, 'fitszilla', self.m_feed, _summary, self.m_outputPath)        
        # 
        self.m_logger.info("Subscan numbers " + str(len(self.m_subscans)))

    def geometry_group(self):
        """
        Group scans by supposed geometry
        First group is made by file sequence
        Geometry has fixed number of file for every 
        """
        _loop_len= fitslike_scangeometry.get_geo_loop_len(self.m_geometry)
        # Subscan filename sorting without the summary file
        self.m_subscans_geo= sorted( ( f for f in self.m_subscans if 'summary' not in f['file_name']),\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))
        # Ceck if number of scans is multiple of geo file len
        if len(self.m_subscans_geo) % _loop_len != 0:
            self.m_logger.error("we have {} subscans and a geometry loop of {}, MCM criteria not satisfied!".format(len(self.m_subscans_geo), _loop_len))
            self.critical_error= True
            return
        # Group sorted file by geo loop len        
        num_geo_groups= len(self.m_subscans_geo) / _loop_len
        self.m_logger.info("Grouping subscans in {} repetitions".format(num_geo_groups))
        self.m_geometry_scan_list= [s.tolist() for s in np.array_split(self.m_subscans_geo, _loop_len) ]
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
            # subscan here is:
            #             
            self.m_subscans= sorted(self.m_subscans,\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))
            subscan_to_work= self.m_subscans
                
                
        # Recognition and grouping for ON,OFF,CAL
        for subscan in subscan_to_work:
            _file_name= ''
            if 'file_name' in subscan:
                _file_name= subscan['file_name']
                self.m_logger.info(_file_name)           
            # Analyzing every subscan
            for _feed in subscan:
                # Work only selected feed
                if self.m_feed:
                    if not self.m_feed in str(_feed):
                        continue
                # Prepare on off group keys hirearchy
                if _feed not in self.m_group_on_off_cal.keys():
                    self.m_group_on_off_cal[_feed]={}
                # Section (_section) navigation
                for _section in subscan[_feed]:
                    if 'ch_' not in _section:
                        continue
                    _chObj= subscan[_feed][_section]
                    # on off group keys hierarchy
                    if _section not in self.m_group_on_off_cal[_feed].keys():
                        # replicate section data to group on off cal to retrieve them easely
                        self.m_group_on_off_cal[_feed][_section]= _chObj
                    _feed= _chObj['frontend']['feed']
                    if 'pol_tables_dict' not in _chObj.keys():
                        continue
                    _po_tables_dict= _chObj['pol_tables_dict']
                    if 'pols' not in self.m_group_on_off_cal[_feed][_section].keys():
                        self.m_group_on_off_cal[_feed][_section]['pols']={}
                    _pols_dict= self.m_group_on_off_cal[_feed][_section]['pols']
                    for _pol in _po_tables_dict.keys():
                        if _pol not in _pols_dict.keys():
                            _pols_dict[_pol]= {'signal': [],'reference': [],'cal_on':[],'cal_off': []}
                        try:
                            _po_table= QTable.read(_po_tables_dict[_pol], memmap= True)
                            # ID data type by columns, creating a list made by
                            if 'cal_on' in _po_table.colnames:
                                _pols_dict[_pol]['cal_on'].append(_po_tables_dict[_pol])
                            elif 'cal_off' in _po_table.colnames:
                                _pols_dict[_pol]['cal_off'].append(_po_tables_dict[_pol])
                            elif 'signal' in _po_table.colnames:
                                _pols_dict[_pol]['signal'].append(_po_tables_dict[_pol])
                            elif 'reference' in _po_table.colnames:
                                _pols_dict[_pol]['reference'].append(_po_tables_dict[_pol])
                        except Exception as e:
                            self.m_logger.error("{}-{}-{}-reading disk table exception : {}".format(_feed, _section, _pol,str(e)))

        # Feed base work dir cleaning
        self.m_group_path= os.path.join(self.m_outputPath, 'fits_groups')
        if os.path.exists(self.m_group_path):
            shutil.rmtree(self.m_group_path)
        os.mkdir(self.m_group_path)
        # Join tables from on off group
        for _feed in self.m_group_on_off_cal.keys():
            for _section in self.m_group_on_off_cal[_feed].keys():
                if 'ch_' not in _section:
                    continue
                 # Work dir creation per feed
                _this_group_path= os.path.join(self.m_group_path, 'group_feed_{}'.format(_feed))
                if not os.path.exists(_this_group_path):
                    os.mkdir(_this_group_path)
                #
                #pdb.set_trace()
                #
                for _pol in self.m_group_on_off_cal[_feed][_section]['pols'].keys():
                    # Join by pol and cal signal
                    for _type in self.m_group_on_off_cal[_feed][_section]['pols'][_pol].keys():
                        try:
                            _table_files= self.m_group_on_off_cal[_feed][_section]['pols'][_pol][_type]
                            if not _table_files:
                                continue
                            if len(_table_files)==0:
                                continue
                            _opened_tables= [QTable.read(t, memmap= True) for t in _table_files]
                            joined= vstack(_opened_tables)
                            _fname= "{}_{}_{}_{}.fits".format(_feed, _section, _pol, _type)
                            _gr_table_path= os.path.join(_this_group_path, _fname)
                            self.m_group_on_off_cal[_feed][_section]['pols'][_pol][_type]= _gr_table_path
                            joined.write(_gr_table_path)
                            # Deleting input tables
                            # for input_file in _table_files:
                            #     try:
                            #         os.remove(input_file)
                            #     except IOError as e:
                            #         self.m_logger.error("Can't delete intermediate file {} : {}".format(input_file, e))
                            self.m_logger.info("{}_{}_{}_{}".format(_feed, _section, _pol, _type))
                            self.m_logger.info("Grouping table into {}".format(_gr_table_path))
                        except Exception as e:
                            self.m_logger.error("{}-{}-{}-{} error writing group table on disk : {}".format(_feed, _section, _pol, _type, str(e)))
        # Deleting intermediate folder
                

    def normalize(self, p_type:str ="cal") -> None:
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

        p_type: str
        Calculation type, counts or kelvin
        "cal": kelvin
        "counts": counts

        If cal_on/cal_off are missing or off are missing simply it skips calculations,
        and provides only signal data table to:
            - output_path/fits_normalized
        Data output will be grouped by feed_sectiìon_pol

        """

        # Norm output dir creation
        self.m_norm_path= os.path.join(self.m_outputPath, 'fits_normalized')
        if os.path.exists(self.m_norm_path):
            shutil.rmtree(self.m_norm_path)
        os.mkdir(self.m_norm_path)
        # Feed traverse
        for _feed in self.m_group_on_off_cal.keys():
            # Work only selected feed, or All
            if self.m_feed:
                if self.m_feed not in str(_feed):
                    continue
            
            if not str(_feed).isdigit():
                continue
            # Feed folder
            _feed_path= os.path.join(self.m_norm_path, 'normalized_feed_{}'.format(_feed))
            if not os.path.exists(_feed_path):
                os.mkdir(_feed_path)
            # Section traverse
            for ch in self.m_group_on_off_cal[_feed].keys():
                _section= self.m_group_on_off_cal[_feed][ch]
                _section_pols= self.m_group_on_off_cal[_feed][ch]['pols']
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    _calMarkTemp= None
                    if pol not in _section_pols.keys(): continue
                    # Get calibration mark temperature if available
                    try:
                        _calMarkTemp= _section['frontend']['cal_mark_temp']
                    except:
                        self.m_logger.warning("[{}][{}][{}] no cal mark temp available".format(_feed, ch, pol))

                    # Here i have feed - section - pol - on/off/cals
                    # Read tables on off cal_on cal_off
                    _opened_tables={
                        'signal': None,
                        'reference' : None,
                        'cal_on': None,
                        'cal_off': None
                        }
                    _data={
                        'signal': QTable(),
                        'reference' : QTable(),
                        'cal_on': QTable(),
                        'cal_off': QTable()
                        }                    
                    for on_off_cal in _opened_tables.keys():
                        if on_off_cal not in _section_pols[pol].keys():
                            continue
                        try:                            
                            _table_path= _section_pols[pol][on_off_cal]
                            if not len(_table_path):
                                continue                            
                            _opened_tables[on_off_cal]= QTable.read(_table_path, memmap= True)
                        except Exception as e:
                            self.m_logger.warning("[{}][{}][{}][{}]".format(_feed, ch, pol, on_off_cal))
                            self.m_logger.error("Error reading goruped table {} - {}".format(_table_path, e))

                    # Now i've loaded all tables needed to calibrate one section
                    # Processing reference mean
                    try:                        
                        # Signal whole data set
                        _data['signal']= _opened_tables['signal']['data']
                        if not _opened_tables['signal']:                            
                            continue
                        if not len(_opened_tables['signal']):
                            continue                        
                        if _opened_tables['reference'] != None:
                            if len(_opened_tables['reference']) != 0:
                                _opened_tables['reference']= _opened_tables['reference'].group_by(['pol']).groups.aggregate(np.mean)
                                _data['reference']= _opened_tables['reference']['data']
                        # Processing cal_on mean norm
                        if _opened_tables['cal_on'] != None and _opened_tables['reference'] != None:
                            if len(_opened_tables['cal_on']) != 0 and len(_opened_tables['reference']) != 0 :
                                _opened_tables['cal_on']= _opened_tables['cal_on'].group_by(['pol']).groups.aggregate(np.mean)
                                _data['cal_on']= _opened_tables['cal_on']['data']
                                _data['cal_on']= (_data['cal_on'] - _data['signal'][0]) / _data['signal'][0]
                        # Processing cal_off mean norm
                        if _opened_tables['cal_off'] != None and _opened_tables['reference'] != None:
                            if len(_opened_tables['cal_off']) != 0 and len(_opened_tables['reference']) != 0 :
                                _opened_tables['cal_off']= _opened_tables['cal_off'].group_by(['pol']).groups.aggregate(np.mean)
                                _data['cal_off']= _opened_tables['cal_off']['data']                                            
                                _data['cal_off']= (_data['cal_off'] - _data['reference']) / _data['reference']
                    except KeyError as e:
                        self.m_logger.warning("[{}][{}][{}]".format(_feed,ch,pol))
                        self.m_logger.error("Missing mandatory keyword {}".format(e))
                    except Exception as e:
                        self.m_logger.warning("[{}][{}][{}]".format(_feed,ch,pol))
                        self.m_logger.error("Exception averaging grouped table {}".format(e))

                    # Have we data?
                    _data_raw= False
                    _on_off_possible = False
                    _calibration_possible= False
                    if len(_data['signal']) != 0:
                        _data_raw= True
                        self.m_logger.info("[{}][{}][{}] Signal data are present".format(_feed,ch,pol))
                        # Sub calculations, on - off
                        if len(_data['reference']) and len(_data['signal']):                            
                            # Check if calculations are possible
                            # TODO implementare solo conteggi senzi cal
                            _on_off_possible= True
                            self.m_logger.info("[{}][{}][{}] Reference data are present".format(_feed,ch,pol))
                            _calibration_possible= False
                            if _calMarkTemp:                                
                                if len(_data['cal_on']) or len(_data['cal_off']):
                                    _calibration_possible= True
                                    self.m_logger.info("[{}][{}][{}] Calibration data are present".format(_feed,ch,pol))
                                else:
                                    self.m_logger.warning("[{}][{}][{}] No cal_on and cal_of available?".format(_feed,ch,pol))    
                            else:
                                self.m_logger.warning("[{}][{}][{}] No calibration mark temperature available".format(_feed,ch,pol))
                    #
                    if not _data_raw:
                        self.m_logger.warning("[{}][{}][{}] No input data, skipping".format(_feed,ch,pol))
                        continue
                    self.m_logger.warning("[{}][{}][{}] Applying calibration".format(_feed,ch,pol))
                    # TODO subs signal table with calibrated data
                    # Calc with units..review
                    # Keep data as small as possible (numpy.float32)
                    try:                        
                        if _calibration_possible and p_type == "cal":                            
                            # Complete calibration
                            # On Off
                            self.m_logger.info("[{}][{}][{}] Trying Kelvin calculations".format(_feed,ch,pol))
                            _signal= _data['signal']
                            _reference= _data['reference']                                                        
                            _data['on_off']= (_signal -_reference) / _reference                            
                            # Calibration factor
                            # Np arrays, assuming data from QTable has shape[[data]]
                            _calon= np.array([])
                            try:
                                _calon= _data['cal_on'].data[0]
                            except:
                                pass
                            _caloff= np.array([])
                            try:
                                _caloff= _data['cal_off'].data[0]
                            except:
                                pass
                            _cal = np.concatenate((_calon,_caloff))
                            if len(_cal) == 0:
                                continue
                            good =  ~np.isnan(_cal) & ~np.isinf(_cal)
                            _cal = _cal[good]                                                        
                            if len(_cal) > 0:
                                meancal = np.median(_cal) if len(_cal) > 30 else np.mean(_cal)                                                                
                                calibration_factor = 1 / meancal * _calMarkTemp                                             
                                calibration_factor= np.float32(calibration_factor.value)* calibration_factor.unit
                            else:
                                continue                              
                            _data['calibrated']= _data['on_off'].astype(np.float32) * calibration_factor                    
                            # Replace data on table                             
                            _opened_tables['signal']['data']= _data['calibrated']
                            # Remove useless columns
                            if 'signal' in _opened_tables['signal'].colnames:
                                del _opened_tables['signal']['signal']                            
                            # Rewrite table on disk, overwriting signal table
                            _file_name= "{}_{}_{}_calibrated.fits".format(_feed, ch, pol)
                            _file_path= os.path.join(_feed_path , _file_name)                                   
                            _section_pols[pol]['calibrated']= _file_path
                            _opened_tables['signal'].write(_file_path, overwrite= True)                                                        
                        elif _on_off_possible:
                            self.m_logger.info("[{}][{}][{}] Trying Counts calculations".format(_feed,ch,pol))
                            # Only counts
                            _signal= _data['signal']
                            _reference= _data['reference']
                            _data['on_off']= (_signal -_reference) / _reference                            
                            _opened_tables['signal']['data']= _data['on_off']
                            # Remove useless columns
                            if 'signal' in _opened_tables['signal'].colnames:
                                del _opened_tables['signal']['signal']                
                            # Writing to disk 
                            _file_name= "{}_{}_{}_counts.fits".format(_feed, ch, pol)
                            _file_path= os.path.join(_feed_path , _file_name)                            
                            _section_pols[pol]['counts']= _file_path
                            _opened_tables['signal'].write(_file_path, overwrite= True)    
                        elif len(_opened_tables['signal']): 
                            # Signal with no calculations                                                                                
                            # Remove useless columns
                            if 'signal' in _opened_tables['signal'].colnames:
                                del _opened_tables['signal']['signal']                
                            # Writing to disk 
                            _file_name= "{}_{}_{}_signal.fits".format(_feed, ch, pol)
                            _file_path= os.path.join(_feed_path , _file_name)                            
                            _section_pols[pol]['signal']= _file_path
                            _opened_tables['signal'].write(_file_path, overwrite= True)    
                        else:
                            # Nothing to do..
                            self.m_no_cal= True
                            self.m_logger.error("[{}][{}][{}] No appropriate data to normalize this dataset".format(_feed, ch, pol))   
                        # Write data
                    except IOError as e:
                        self.m_no_cal= True                        
                        self.m_logger.error("[{}][{}][{}] Exception applying calibration to this data set : {}".format(_feed,ch,pol,e))                           
                    except Exception as e:
                        self.m_no_cal= True            
                        self.m_logger.error("[{}][{}][{}] Exception applying calibration to this data set : {}".format(_feed,ch,pol,e))                   
                    
        # Delete group data
        # if os.path.exists(self.m_group_path):
        #     shutil.rmtree(self.m_group_path, ignore_errors= True)
                        
                    
    def ClassFitsAdaptations(self, p_type:str ="cal")->None:
        """
        p_type: str
        Conversione dati kelvin oppure altro tipo
        "cal": dati kelvin
        "counts": conteggi

        Raw data groups stored as :
            self.m_group_on_off_cal[_feed][_section]['pols'][LL LR RL RR][on off cal_on cal off]

        On disk we have a folder relative to output path:
            output_path/fits_groups, folder where table data (disk part) are grouped by feed section pol on off cal
            output_path/fits_norms, folder where grouped data produce a calibrated measures (if cal is present, counts otherwise)
            output_path/classfits, single, calibrated classfits file boxing every feed and pol

        Conversion tips:
            desired cooord in az, el o ra, dec
            observed coord in crdelt2,3
        """

        def getQTableColWithUnit(p_table, p_col, p_unit):
            """ Get QTable[col_u_unit] data adding unit to returned data """
            _col_names= p_table.colnames
            for col_name in _col_names:
                if p_col in col_name:
                    _field_unit= col_name.split('_u_')
                    if len(_field_unit) == 2:
                        data= p_table[col_name] * unit.Unit(_field_unit[1])
                        data= data.to(p_unit).value
                        return data

        # Work folder
        _folder_name= 'classfits_{}'.format(p_type)
        self.m_class_path= os.path.join(self.m_outputPath, _folder_name)
        if os.path.exists(self.m_class_path):
            shutil.rmtree(self.m_class_path)
        os.mkdir(self.m_class_path)
        # Feed traverse
        for _feed in self.m_group_on_off_cal:
            # Feed filter added by user ?
            if self.m_feed:
                # Check if selected feed is the feed we want to work
                if self.m_feed not in str(_feed):
                    continue
            # Only feed, avoid extra data on this level
            if not str(_feed).isdigit():
                continue            
            # Feed folder
            _feed_path= os.path.join(self.m_class_path, 'class_feed_{}'.format(_feed))
            if not os.path.exists(_feed_path):
                os.mkdir(_feed_path)
            # Sections traverse            
            for _section in self.m_group_on_off_cal[_feed]:
                _section_data= self.m_group_on_off_cal[_feed][_section]
                #pdb.set_trace()
                _section_pols= self.m_group_on_off_cal[_feed][_section]['pols']
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    if pol not in _section_pols.keys(): continue                    
                    # Check work on data on only ( Raw Mode )                    
                    calibration_type= ''
                    _polarization_table= None                                          
                    # ONLY 1 table per pol!
                    try:                                                             
                        # Data retreiving one level up from 'on' data
                        #pdb.set_trace()
                        if p_type == "cal":
                            if 'calibrated' in _section_pols[pol].keys():
                                calibrated_data_path= _section_pols[pol]['calibrated']                                
                                self.m_logger.info("Loading normalized {}-{}-{} path {}".format(_feed, _section, pol,calibrated_data_path))
                                if calibrated_data_path: 
                                    _polarization_table= QTable.read(calibrated_data_path, memmap= True)
                                    calibration_type= 'calibrated'
                                else:
                                    self.m_logger.info("calibrated for {}-{}-{} is empty".format(_feed, _section, pol))                                    
                                
                        if p_type == "counts":
                            if  not _polarization_table and 'counts' in _section_pols[pol].keys():
                                on_off_data_table_path= _section_pols[pol]['counts']
                                self.m_logger.info("Loading counts {}-{}-{} path {}".format(_feed, _section, pol,on_off_data_table_path))
                                if on_off_data_table_path:
                                    _polarization_table= QTable.read(on_off_data_table_path, memmap= True)
                                    calibration_type= 'counts'                                
                                else:
                                    self.m_logger.info("on_off data for {}-{}-{} is empty".format(_feed, _section, pol))    
                                
                        if not _polarization_table and 'signal' in _section_pols[pol].keys():
                            signa_data_path= _section_pols[pol]['signal']
                            self.m_logger.info("Loading signal {}-{}-{} path {}".format(_feed, _section, pol,signa_data_path))
                            if signa_data_path:
                                _polarization_table= QTable.read(signa_data_path, memmap= True)
                                calibration_type= 'signal'
                            else:
                                self.m_logger.info("signal data for {}-{}-{} is empty".format(_feed, _section, pol))                                    
                                
                        if  not _polarization_table and 'reference' in _section_pols[pol].keys():
                            reference_data_path= _section_pols[pol]['reference']
                            self.m_logger.info("Loading reference {}-{}-{} path {}".format(_feed, _section, pol,reference_data_path))
                            if reference_data_path:
                                _polarization_table= QTable.read(reference_data_path, memmap= True)
                                calibration_type= 'reference'
                            else:
                                self.m_logger.info("reference data for {}-{}-{} is empty".format(_feed, _section, pol))                                                                    
                                
                        # No available input data options..
                        if not _polarization_table:
                            self.m_logger.warning("{}_{}_{} No data suitable for class conversion".format(_feed, _section, pol))
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
                    self.m_obs_genera_data={}
                    self.m_obs_genera_data['ra']= _section_data['scheduled']['ra'].to(unit.deg).value
                    self.m_obs_genera_data['dec']= _section_data['scheduled']['dec'].to(unit.deg).value
                    self.m_obs_genera_data['source']= _section_data['scheduled']['source']
                    self.m_obs_genera_data['date-red']= Time.now().to_datetime().strftime('%d/%m/%y')
                    # classfits new dict filling process
                    _class_data= {}
                    # ch by ch
                    try:
                        # Data first to get shapes                            
                        _class_data['SPECTRUM']= _polarization_table['data']
                        # data shape they must be equals ( on_off cal on )
                        data_shape= _class_data['SPECTRUM'].shape
                        _class_data['CRPIX1']=  _section_data['backend']['bins'] // 2 + 1
                        # Lavoro con i dati integrati
                        # ut
                        _tMjd= _polarization_table['data_time_mjd']
                        _timeMjd= Time(_tMjd, format='mjd', scale='utc')                        
                        _class_data['UT']= ( _tMjd - np.floor(_tMjd)) * 86400                        
                        # date
                        _class_data['DATE-OBS']= _timeMjd.strftime('%Y-%m-%dT00:00:00.000')                        
                        # lsts
                        _lsts= _timeMjd.sidereal_time('apparent', \
                                  fitslike_commons.Fitslike_commons.\
                                      get_site_location(_section_data['scheduled']['antenna']).lon)
                        _lsts= _lsts.value * unit.hr
                        # infos
                        _class_data['OBJECT']= _section_data['scheduled']['source']
                        _class_data['LINE']= "F{}-{:3.3f}-MHz"\
                            .format(_section_data['frontend']['feed'], _section_data['backend']['bandwidth'])
                        self.m_obs_genera_data['line']= _class_data['LINE']
                        try:
                            #pdb.set_trace()
                            _class_data['TELESCOP']=\
                                self.m_commons.class_telescope_name(_section_data,self.m_summary['summary']['backend_name'], pol, _section)
                        except ValueError as e:
                            self.m_logger.error("Missing summary !? TELESCOP set to empty value \n Excpetion: "+ str(e))
                            _class_data['TELESCOP']= ""
                        _mH2O= _polarization_table['weather']
                        _class_data['MH2O']= _mH2O
                        # temp
                        _class_data['TSYS']= 1.0
                        _class_data['CALTEMP']= _section_data['frontend']['cal_mark_temp'].value
                        # time
                        _class_data['LST'] = _lsts.to('s').value
                        # CDELT
                        _class_data['CDELT1']= (_section_data['frontend']['bandwidth'] / _section_data['backend']['bins']).to('Hz')
                        # freq and velocity
                        #pdb.set_trace()
                        _class_data['RESTFREQ']= self.m_summary['summary']['restfreq'].to(unit.Hz).value                        
                        # Caso semplificato (<1) se la rest freq è a 0 uso la frequenza del canale bins/2+1
                        if _class_data['RESTFREQ'] != None:
                            if _class_data['RESTFREQ'] < 1:
                                _class_data['RESTFREQ']= _section_data['frontend']['local_oscillator'].to('Hz') +\
                                                    _class_data['CDELT1'] * (_section_data['backend']['bins']/2 +1)
                                _class_data['RESTFREQ']= _class_data['RESTFREQ'].value
                        self.m_obs_genera_data['restfreq']= _class_data['RESTFREQ']
                        # Versione chiave singola VELO-FRAME
                        #  self.m_summary['summary']['velo_frame']
                        # Dato velocità da summary.VRAD
                        # se manca sul summary : _section_data['scheduled']['vlsr']
                        _velo_keyword="VELO-{}".format(self.m_summary['summary']['velo_frame'])
                        _class_data['VELOCITY']= _velo_keyword
                        self.m_logger.info("Velo column name {}".format(_velo_keyword))
                        if self.m_summary['summary']['velo_rad'] != .0:
                            self.m_logger.warn("Velocity taken from summary.velo_rad {}".format(self.m_summary['summary']['velo_rad']))
                            _class_data[_velo_keyword]= self.m_summary['summary']['velo_rad'].to("m/s").value                        
                        else:
                            self.m_logger.warn("Velocity taken from observed.vlsr {}".format(_section_data['scheduled']['vlsr']))
                            _class_data[_velo_keyword]= _section_data['scheduled']['vlsr'].to("m/s").value
                        _df= (_section_data['backend']['bandwidth'] / _section_data['backend']['bins']).to('Hz')
                        _class_data['CDELT1']= _df.value
                        self.m_obs_genera_data['cdelt1']= _class_data['CDELT1']
                        _deltav= - _df/ _class_data['RESTFREQ'] * const.c
                        _class_data['DELTAV']= _deltav.value
                        # LOG test
                        #self.m_logger.warn("RESTFREQ {}".format(_class_data['RESTFREQ']))
                        # Objects Coordinates   
                        # CRVAL 2,3: target object coordiantes                        
                        _class_data['CRVAL2']= self.m_summary['summary']['target_ra'].to(unit.deg).value
                        _class_data['CRVAL3']= self.m_summary['summary']['target_dec'].to(unit.deg).value
                        # CDELT 2,3: observer - target
                        _observed_ra= getQTableColWithUnit(_polarization_table,'data_ra', unit.deg)                        
                        _offset_ra= (_observed_ra - _class_data['CRVAL2'])* np.cos(_class_data['CRVAL3'])
                        _class_data['CDELT2']= _offset_ra
                        _observed_dec= getQTableColWithUnit(_polarization_table,'data_dec', unit.deg)                        
                        _offset_dec= (_observed_dec - _class_data['CRVAL3'])
                        _class_data['CDELT3']= _offset_dec
                        # az el deg
                        _class_data['AZIMUTH']= getQTableColWithUnit(_polarization_table,'data_az', unit.deg)
                        _class_data['ELEVATIO']= getQTableColWithUnit(_polarization_table,'data_el', unit.deg)
                        # data
                        _class_data['OBSTIME'] = _section_data['backend']['integration_time'].to('s').value
                        _class_data['MAXIS1'] = _section_data['backend']['bins']
                        self.m_obs_genera_data['maxis1']= _class_data['MAXIS1']
                        # we have to shape classfits data properly according to data shape
                        # es data shape is 12, 16384 we have to replicate data this shape
                        rows= data_shape[0]
                        for k in _class_data.keys():
                            _value= _class_data[k]
                            if "SPECTRUM" in k:
                                continue
                            try:
                                # this prevent already shaped (rows,) to be warped
                                _class_data[k]= np.full((rows,), _value)
                            except Exception  as e:
                                self.m_logger.error("Error reshaping classfits data: " +str(e))
                        # Transform dictionary based converted table in QTable on disk
                        _class_table= QTable()
                        for col in _class_data.keys():
                            try:
                                # Skip empty columns
                                _class_table.add_column( Column(data= _class_data[col], name= col) )
                            except ValueError as e:
                                self.m_logger.error("Exception adding column [{}] to classfit converted table {}".format(col, e))
                        # Table name generation                                                
                        _class_name= '{}_{}_{}_{}_class.fits'.format(_feed, _section, pol, calibration_type)       
                        _class_path= os.path.join(_feed_path, _class_name)
                        self.m_logger.info("class table path: {}".format(_class_path))
                        try:                            
                            self.classfitsWrite(_class_path, _class_table, p_type)
                        except ValueError as e:                            
                            self.m_logger.error("Error writing small classfit table to disk {}".format(e))                       
                    except Exception as e:                        
                        self.m_logger.error("[Exception] Error preparing class data: {}".format(e))                    
                        self.m_logger.error("\nEXCEPTION\n-------------------------------------")
                        traceback.print_exc()
                        self.m_logger.error("\n-------------------------------------")                                                                    

        # Delete norm data
        # try:
        #     if os.path.exists(self.m_norm_path) and self.m_norm_path:
        #         shutil.rmtree(self.m_norm_path, ignore_errors= True)
        # except Exception as e:
        #         self.m_logger.error("Error cleaning normalized data folder: {}\n {}\n".format(self.m_norm_path, e))


    def classfitsWrite(self, p_file_name, p_class_table, p_cal_type:str ="cal"):
        """
        Convert Qtable classfits column to fits file format adding appropriate
        header data

        Parameters
        ---------
        p_file_name: string
            Output file name
        p_class_table: string
            QTable with clasfits data            
        p_cal_type: string
            Calibratin type
        """              
        _newCols=[]
        # for every column expressed in classfits definition ..
        for classCol in self.m_commons.getClassfitsColumnsZip():            
            # [ name, form, unit ] column by column data building 
            # fill one column looking into every feed[on], and builds column data 
            _columnFound = False
            # conditionals, some fields needs dedicated approach            
            # converted fits data matches with classfits columns? 
            # some columns need special care
            _inferredCol= classCol[0]
            if _inferredCol in p_class_table.colnames:                
                # found match, add data to column data 
                _columnFound= True                
            try:
                # adding column to classfits if fitszilla representation matches it
                if _columnFound:                   
                    # some fields needs dedicated approach 
                    if classCol[0] == "SPECTRUM":
                        # checking calibration type, then adjust SPECTRUM's unit
                        _unit= classCol[2]
                        if p_cal_type == "counts":
                            _unit= "Counts"  
                        _rows= p_class_table[_inferredCol][0].shape[0]
                        _newCols.append(fits.Column(array= p_class_table[_inferredCol],name= classCol[0],format= "{}D".format(_rows),unit= _unit))
                    elif classCol[0] == "VELOCITY":
                        _velo_name= p_class_table[classCol[0]][0]                                                 
                        self.m_logger.warning("Velo classfit column name {}".format(_velo_name))                            
                        _newCols.append(fits.Column(array= p_class_table[_velo_name],name= _velo_name,format= classCol[1],unit= classCol[2]) )
                    else:
                        _newCols.append(fits.Column(array= p_class_table[_inferredCol],name= classCol[0],format= classCol[1],unit= classCol[2]) )

            except Exception as e:
                self.m_logger.error("Classfits column creation exception: {} - {}".format(classCol, e))    
                traceback.print_exc()            

        _hdData= self.m_obs_genera_data
        # header 
        _hdu= fits.PrimaryHDU()        
        # data 
        try:
            # TEST 
            #for col in _newCols:
            #    print(col.name + " " + str(col.array.shape))
            _cdefs= fits.ColDefs(_newCols)            
            _hdu= fits.BinTableHDU().from_columns(_cdefs)
        except Exception as e:            
            self.m_logger.error("Exception creating classfits model file {} - {}".format(p_file_name, e))            
            return
        #header 
        try:
            # TODO Rivedere i campi DELT, PIX ...
            _hdu.header['EXTNAME'] = "MATRIX"
            _hdu.header['EXTVER'] = 1
            _hdu.header['MAXIS'] = 4
            if 'maxis1' in _hdData.keys():
                _hdu.header['MAXIS1'] = _hdData['maxis1']
            else:
                _hdu.header['MAXIS1'] = 0.0
            _hdu.header['MAXIS2'] = 4
            _hdu.header['MAXIS3'] = 4
            _hdu.header['MAXIS4'] = 4            
            _hdu.header['CTYPE1']= "FREQ"
            _hdu.header['CRVAL1']= 0.0
            if 'cdelt1' in _hdData.keys():
                _hdu.header['CDELT1'] = _hdData['cdelt1']
            else:
                _hdu.header['CDELT1'] = 0.0
            _hdu.header['CTYPE2']= "RA"
            if 'ra' in _hdData.keys():
                _hdu.header['CRVAL2']= _hdData['ra']
            else:
                _hdu.header['CRVAL2']= 0.0
            _hdu.header['CDELT2']= 0.0
            _hdu.header['CRPIX2']= 0.0
            _hdu.header['CTYPE3']= "DEC"
            if 'dec' in _hdData.keys():
                _hdu.header['CRVAL3']= _hdData['dec']
            else:
                _hdu.header['CRVAL3']= 0.0            
            _hdu.header['CDELT3']= 0.0
            _hdu.header['CRPIX3']= 0.0
            _hdu.header['CTYPE4']= "STOKES"
            _hdu.header['CRVAL4']= 0.0
            _hdu.header['CDELT4']= 0.0
            _hdu.header['CRPIX4']= 0.0
            _hdu.header['SUBSCAN']= 1
            if 'line' in _hdData.keys():
                _hdu.header['LINE'] = _hdData['line']
            else:
                _hdu.header['LINE'] = " "
            if 'source' in _hdData.keys():
                _hdu.header['OBJECT'] = _hdData['source']            
            else:
                _hdu.header['OBJECT'] = " "
            if 'restfreq' in _hdData.keys():
                _hdu.header['RESTFREQ'] = _hdData['restfreq']
            else:
                 _hdu.header['RESTFREQ'] = 0.0
            # TODO Controllare la composizione (ora fissa LSR)
            #_hdu.header['VELDEF']= self.m_summary['summary']['velo_def'] + '-' +
            _hdu.header['GAINIMAG']= 0.0
            _hdu.header['BEAMEFF']= 0.0
            _hdu.header['FORMEFF']= 0.0
            if 'restfreq' in _hdData.keys():
                _hdu.header['DATE-RED'] = _hdData['date-red']        
            else:
                _hdu.header['DATE-RED'] = " "
            _hdu.header['EPOCH']= 2000.0
            _hdu.header['CRVAL'] = 0
        except KeyError as e:
            self.m_logger.error("Exception filling " + p_file_name + " header data: "+ str(e))
        # Disk writing
        try:
            #if os.path.exists(p_file_name):
            #    os.remove(p_file_name)
            _hdu.writeto(p_file_name)
            self.m_logger.info("Wrote classfits {}".format(p_file_name))
        except Exception as e:            
            self.m_logger.error("Exception writing fitsl  file {} - {}".format(p_file_name, e))
