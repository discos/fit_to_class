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
from collections import defaultdict
from multiprocessing import Pool
import sys
import os
import shutil
import traceback
import pdb

from fit_to_class import fitslike
from fit_to_class  import awarness_fitszilla
from fit_to_class.fitslike_commons import keywords as kws
from fit_to_class import fitslike_commons


def _envelope_subscan(p_logger, p_ftype, p_feed, p_path) -> dict:
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
    p_logger.info("Scanning " + p_path)
    l_path, l_filename= os.path.split(p_path)
    p_logger.info("fitszilla parsing: " + p_path)
    try:
        l_fits = fits.open(p_path)
    except Exception as e:
        p_logger.error("Skipping input file {} due to an excpetion: \n{}".format(p_path, e))
        return {}
    l_aware= awarness_fitszilla.Awarness_fitszilla(l_fits, p_path, p_feed)
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

    def __init__(self, p_raw, p_scanType, p_feed, p_files_per_scan= 3):
        """
        Internal struct definition

        Parameters
        ----------

        p_raw : bool
            Avoid grouping data by on off cal
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
        self.m_raw= p_raw
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

    def scan_data(self, p_dataDir):
        """
        Takes data input directory and launchs subscan data conversion
        to fitslike

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
        for dirpath, dirnames, filenames in os.walk(self.m_dataDir):
            l_inputFiles.extend(filenames)
        l_inputFiles= [f for f in l_inputFiles if re.match(l_filt,f)]
        # Split parsing multiprocessing
        self.m_results=[]
        l_poolSize= self.m_files_per_pool
        for l_group in chunker(l_inputFiles, l_poolSize ):
            l_results=[]
            self.m_pool=Pool(l_poolSize)
            for l_fPath in l_group:
                l_results.append( self.m_pool.apply_async(
                        _envelope_subscan,
                        [self.m_logger,
                        'fitszilla',
                        self.m_feed,
                        p_dataDir +'/'+ l_fPath
                        ])
                    )
            self.m_pool.close()
            self.m_pool.join()
            self.m_subscans= self.m_subscans + [x.get() for x in l_results]
        self.m_logger.info("subscan numbers " + str(len(self.m_subscans)))

    def group_on_off_cal(self):
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

        def _getGroupPol( p_group):
            """ Get Polarization flag """
            return p_group['pol'][0]

   

        " ordino le subscan in base al file name "
        self.m_subscans= sorted(self.m_subscans,\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))

        # Work dir
        os.chdir(self.m_outputPath)
        # Feed base work dir cleaning
        l_base_path= 'tmp_feed_*'
        if os.path.exists(l_base_path):
            shutil.rmtree(l_base_path)
        # Subscan data grouping
        for l_subscan in self.m_subscans:
            l_file_name= ''
            if 'file_name' in l_subscan:
                l_file_name= l_subscan['file_name']
                self.m_logger.info(l_file_name)
            # Summary file doesn't group
            if 'summary' in l_subscan.keys():
                self.m_summary= l_subscan
                self.m_logger.info("summary.fits excluded from on off cal")
                continue
            # Analyzing every subscan
            for l_feed in l_subscan:
                # Work only selected feed
                if self.m_feed:
                    if not self.m_feed in str(l_feed):
                        continue
                # Prepare on off group keys hirearchy
                if l_feed not in self.m_group_on_off_cal.keys():
                    self.m_group_on_off_cal[l_feed]={}
                # Work dir creation per feed
                l_base_path= 'tmp_feed_{}/'.format(l_feed)
                if not os.path.exists(l_base_path):
                    os.mkdir(l_base_path)
                # Section (l_chx) navigation
                for l_chx in l_subscan[l_feed]:
                    if 'ch_' not in l_chx:
                        continue
                    l_chObj= l_subscan[l_feed][l_chx]
                    # on off group keys hirearchy
                    if l_chx not in self.m_group_on_off_cal[l_feed].keys():
                        # replicate section data to group on off cal to retrieve them easely
                        self.m_group_on_off_cal[l_feed][l_chx]= l_chObj
                    # subscan type [not a map]
                    #if self.m_scanType == 'on_off' or self.m_scanType == 'nod':
                    l_feed= l_chObj['frontend']['feed']
                    # split by polarization so we have
                    # subscan - feed - ch_x(section) - pol ( L R Q U )
                    # pol navigation
                    if 'groups' not in l_chObj.keys():
                        continue
                    l_grouped_pol= l_chObj['groups'].group_by(['pol'])
                    for group in l_grouped_pol.groups:
                        # Find pol
                        l_pol= _getGroupPol(group)
                        # grouping or not ( Raw Mode, every entr is on )
                        if self.m_raw:
                            # on off group keys hirearchy
                            # list per feed ch_x setup
                            # Astropy fits table on disk creation, this table is a temp table to work
                            # subscan data with the goal of reducing RAM usage and improve computation
                            # stability with big input data
                            #l_onOffdict= {
                            #    'on': [], 'off': [],'cal_on':[], 'cal_off': []
                            #}
                            
                            l_onOffdict= {
                                'on': l_base_path + "section_{}_pol_{}_on.fits".format(l_chx, l_pol),
                                'off': l_base_path + "section_{}_pol_{}_off.fits".format(l_chx, l_pol),
                                'cal_on':l_base_path + "section_{}_pol_{}_cal_on.fits".format(l_chx, l_pol),
                                'cal_off': l_base_path + "section_{}_pol_{}_cal_off.fits".format(l_chx, l_pol)
                            }
                            
                            #pdb.set_trace()
                            if l_pol not in self.m_group_on_off_cal[l_feed][l_chx].keys():
                                self.m_group_on_off_cal[l_feed][l_chx][l_pol]= l_onOffdict
                            try:
                                # Converting qtable to fits table and join existing table
                                # Reading on disk QTable and append to astropy fits table
                                # through conversion                                
                                #l_stack_path= l_base_path + "section_{}_pol_{}_{}.fits".format(l_chx, l_pol, l_file_name)
 #                               self.m_group_on_off_cal[l_feed][l_chx][l_pol]['on'].append(l_stack_path)                                
                                l_disk_table= QTable()
                                if os.path.exists(l_stack_path):
                                    
                                    # TODO 
                                    #l_disk_table= QTable.read(l_stack_path)
                                    #l_disk_table= vstack([l_disk_table, group])
                                    #l_disk_table.write(l_stack_path, overwrite= True)
                                    self.m_logger.info('Overw. existing table {}'.format(l_stack_path))
                                else:
                                    self.m_logger.info('Creating new table {}'.format(l_stack_path))
                                    l_disk_table= QTable(group)
                                    l_disk_table.write(l_stack_path)
                            except Exception as e:
                                self.m_logger.error("{}-{}-{}-on table stacking exception : {}".format(l_feed, l_chx, l_pol,str(e)))
                        else:
                            l_onOffdict= {
                                'on': QTable(), 'off':QTable(),'cal_on':QTable(), 'cal_off':QTable()
                                }
                            # scanning joined table group
                            # cal on off for table group
                            l_isCal= _is_cal(l_chObj, group)
                            l_isOn= _is_on(l_chObj)
                            # on off group keys hirearchy
                            if l_pol not in self.m_group_on_off_cal[l_feed][l_chx].keys():
                                self.m_group_on_off_cal[l_feed][l_chx][l_pol]= l_onOffdict
                            # Grouping feed sub scans on off cal on off
                            if l_isOn and not l_isCal: # ON group
                                l_on= self.m_group_on_off_cal[l_feed][l_chx][l_pol]['on']
                                try:
                                    self.m_group_on_off_cal[l_feed][l_chx][l_pol]['on']= vstack([l_on, group])
                                except Exception as e:
                                    self.m_logger.error("{}-{}-{}-on table stacking exception : {}".format(l_feed, l_chx, l_pol,str(e)))
                                self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is on ')

                            elif not l_isOn and not l_isCal:  # OFF group
                                l_off= self.m_group_on_off_cal[l_feed][l_chx][l_pol]['off']
                                try:
                                    self.m_group_on_off_cal[l_feed][l_chx][l_pol]['off']= vstack([l_off, group])
                                except Exception as e:
                                    self.m_logger.error("{}-{}-{}-off table stacking exception : {}".format(l_feed, l_chx, l_pol,str(e)))
                                self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is off ')
                            elif l_isOn and  l_isCal: # CAL ON group
                                l_calon= self.m_group_on_off_cal[l_feed][l_chx][l_pol]['cal_on']
                                try:
                                    self.m_group_on_off_cal[l_feed][l_chx][l_pol]['cal_on']= vstack([l_calon, group])
                                except Exception as e:
                                    self.m_logger.error("{}-{}-{}-cal_on table stacking exception : {}".format(l_feed, l_chx, l_pol,str(e)))
                                self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is cal on ')
                            elif not l_isOn and  l_isCal: # CAL OFF group
                                l_caloff= self.m_group_on_off_cal[l_feed][l_chx][l_pol]['cal_off']
                                try:
                                    self.m_group_on_off_cal[l_feed][l_chx][l_pol]['cal_off']= vstack([l_caloff, group])
                                except Exception as e:
                                    self.m_logger.error("{}-{}-{}-cal_off table stacking exception : {}".format(l_feed, l_chx, l_pol,str(e)))
                                self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is cal off ')


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

        """

        def _tableToDict(p_table) -> dict:
            """
            Disgregate :) and build a single row dictionary
            Gets unit from column names _u_unit
            """
            "  "
            l_out={}
            for field in p_table.colnames:
                splits= field.split('_u_')
                if len(splits) > 1:
                    l_out[splits[0]]= p_table[field].data.data
                    l_out[splits[0]] *= unit.Unit(splits[1])
                else:
                    try:
                        l_out[field]= p_table[field].data.data
                    except TypeError:
                        " cannot mean strings"
                        l_out[field]= p_table[field][0]
            return l_out


        for l_feed in self.m_group_on_off_cal.keys():
            # Work only selected feed, or All
            if self.m_feed:
                if self.m_feed not in str(l_feed):
                    continue
            for ch in self.m_group_on_off_cal[l_feed].keys():
                l_section= self.m_group_on_off_cal[l_feed][ch]
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    if pol not in l_section.keys(): continue
                    " get calibration mark temperature if available "
                    try:
                        l_calMarkTemp= l_section['frontend']['cal_mark_temp']
                    except:
                        self.m_logger.warning("[{}][{}][{}]['on'] is empty (why?)".format(l_feed,ch,pol))
                        continue
                    " Here i have feed - section - pol - on/off/cals "
                    " i want to decompose table group to a dictionary "
                    for group in l_section[pol].keys():
                        self.m_group_on_off_cal[l_feed][ch][pol][group]= _tableToDict(l_section[pol][group])
                    l_polarization= self.m_group_on_off_cal[l_feed][ch][pol]
                    " flag can calibrate "
                    can_calibrate= True
                    " Avg off "
                    #l_offAvgData= []
                    try:
                        l_offAvg= np.mean(l_polarization['off']['data'], axis= 0)
                    except KeyError:
                        self.m_logger.warning(" off not present in {}-{}-{}".format(l_feed,ch,pol))
                        can_calibrate= False
                    " Avg call on "
                    l_calOnAvg= None
                    try:
                        if l_polarization['cal_on']['data'].shape[0] :
                            l_calOnAvg= np.mean(l_polarization['cal_on']['data'], axis= 0)
                            l_calOnAvg= (l_calOnAvg - l_offAvg) / l_offAvg
                    except KeyError :
                        self.m_logger.warning(" cal_on not present in {}-{}-{}".format(l_feed,ch,pol))
                        can_calibrate= False
                    " Avg cal off "
                    l_calOffAvg= None
                    try:
                        if l_polarization['cal_off']['data'].shape[0]:
                            l_calOffAvg= np.mean(l_polarization['cal_off']['data'], axis= 0)
                            l_calOffAvg = (l_calOffAvg - l_offAvg) / l_offAvg
                    except KeyError:
                        self.m_logger.warning(" cal_off not present in {}-{}-{}".format(l_feed,ch,pol))
                        can_calibrate= False
                    " On - Off "
                    on_off= []
                    try:
                        for elOn in l_polarization['on']['data']:
                            on_off.append( (elOn - l_offAvg)/l_offAvg )
                    except:
                        self.m_logger.warning(" ON data not present in {}-{}-{}".format(l_feed,ch,pol))
                        can_calibrate= False
                    " adding on_off list to section "
                    l_polarization['on_off']= on_off
                    self.m_group_on_off_cal[l_feed][ch][pol]['on_off']= np.array(on_off)
                    " Rescale with cal mark temp "
                    " average non lienarity from receiver at differents power input levels "
                    if not can_calibrate:
                        self.m_logger.error("Cannot calibrate without cal_on, cal_off data ")
                        continue
                    cal = np.concatenate((l_calOnAvg, l_calOffAvg))
                    try:
                        if len(cal) == 0:
                            continue
                        good =  ~np.isnan(cal) & ~np.isinf(cal)
                        cal = cal[good]
                        if len(cal) > 0:
                            meancal = np.median(cal) if len(cal) > 30 else np.mean(cal)
                            calibration_factor = 1 / meancal * l_calMarkTemp
                        else:
                            continue
                        " Calibrated spectrum added to chx "
                        calibrated=[]
                        for elOnOff in l_polarization['on_off']:
                            calibrated.append(elOnOff * calibration_factor)
                        " adding calibrated data to section"
                        self.m_group_on_off_cal[l_feed][ch][pol]['calibrated']= np.array(calibrated)
                    except Exception as e:
                        traceback.print_exc()
                        self.m_no_cal= True
                        self.m_logger.error("Cannot apply calibration to this data set")
                        self.m_logger.error(str(e))

    def ClassFitsAdaptations(self):
        """
        Generazione struttura dati secondo la definizione del classfist

        data in group on off cal sono divisi per
        self.m_group_on_off_cal['ch_0']['on'][0].keys()

        info base

        coordinate comandate in az, el o ra, dec
        coordinate osservate in crdelt2,3
        spettri separati per polarizzazione e per feed
        un file per ogni uno

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



        for l_feed in self.m_group_on_off_cal:
            # Feed filter added by user ?
            if self.m_feed:
                # Check if selected feed is the feed we want to work
                if self.m_feed not in str(l_feed):
                    continue
            classfits= []
            for l_ch in self.m_group_on_off_cal[l_feed]:
                " navigating polarizations "
                l_chx= self.m_group_on_off_cal[l_feed][l_ch]
                for pol in ['LL', 'RR', 'LR', 'RL']:
                    if pol not in self.m_group_on_off_cal[l_feed][l_ch].keys(): continue
                    # We work starting from 'on' data set
                    l_table_list= self.m_group_on_off_cal[l_feed][l_ch][pol]['on']
                    # For every pol it reads related disk Table
                    for l_pol_table_path in l_table_list:
                        l_polarization_table= QTable.read(l_pol_table_path)
                        # Check work on data on only ( Raw Mode )
                        if self.m_raw:
                            calibrated_data= None
                            on_off_data= None
                        else:
                            # Data retreiving one level up from 'on' data
                            if 'calibrated' in self.m_group_on_off_cal[l_feed][l_ch][pol].keys():
                                calibrated_data= self.m_group_on_off_cal[l_feed][l_ch][pol]['calibrated']
                            else:
                                calibrated_data= None
                            if 'on_off' in self.m_group_on_off_cal[l_feed][l_ch][pol].keys():
                                on_off_data= self.m_group_on_off_cal[l_feed][l_ch][pol]['on_off']
                            else:
                                on_off_data= None
                        #pdb.set_trace()
                        # Generic observation data copy to dedicated dict, more copies
                        # below during calculations
                        self.m_obs_general_data={}
                        self.m_obs_general_data['ra']= l_chx['scheduled']['ra'].to(unit.deg).value
                        self.m_obs_general_data['dec']= l_chx['scheduled']['dec'].to(unit.deg).value
                        self.m_obs_general_data['source']= l_chx['scheduled']['source']
                        self.m_obs_general_data['date-red']= Time.now().to_datetime().strftime('%d/%m/%y')
                        # classfits new dict filling process
                        l_class_data= {}
                        # ch by ch
                        try:
                            # Data first to get shapes
                            l_class_data['SPECTRUM_RAW']= l_polarization_table['data']
                            l_class_data['SPECTRUM_ON_OFF']= on_off_data
                            l_class_data['SPECTRUM_CAL']= calibrated_data
                            # data shape they must be equals ( on_off cal on )
                            data_shape= l_class_data['SPECTRUM_RAW'].shape
                            l_class_data['CRPIX1']=  l_chx['backend']['bins'] // 2 + 1
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
                                          get_site_location(l_chx['scheduled']['antenna']).lon)
                            l_lsts= l_lsts.value * unit.hr
                            # infos
                            l_class_data['OBJECT']= l_chx['scheduled']['source']
                            l_class_data['LINE']= "F{}-{:3.3f}-MHz"\
                                .format(l_chx['frontend']['feed'], l_chx['backend']['bandwidth'])
                            self.m_obs_general_data['line']= l_class_data['LINE']
                            try:
                                #pdb.set_trace()
                                l_class_data['TELESCOP']=\
                                    self.m_commons.class_telescope_name(l_chx,self.m_summary['summary']['backend_name'], pol, l_ch)
                            except ValueError as e:
                                self.m_logger.error(str(e))
                            l_mH2O= l_polarization_table['weather']
                            l_class_data['MH2O']= l_mH2O
                            # temp
                            l_class_data['TSYS']= 1.0
                            l_class_data['CALTEMP']= l_chx['frontend']['cal_mark_temp'].value
                            # time
                            l_class_data['LST'] = l_lsts.to('s').value
                            # CDELT
                            l_class_data['CDELT1']= (l_chx['frontend']['bandwidth'] / l_chx['backend']['bins']).to('Hz')
                            # freq and velocity
                            l_class_data['RESTFREQ']= self.m_summary['summary']['restfreq'].to(unit.Hz).value
                            self.m_obs_general_data['restfreq']= l_class_data['RESTFREQ']
                            l_class_data['VELOCITY']= l_chx['scheduled']['vlsr'].to("m/s").value
                            l_df= (l_chx['backend']['bandwidth'] / l_chx['backend']['bins']).to('Hz')
                            l_class_data['CDELT1']= l_df.value
                            self.m_obs_general_data['cdelt1']= l_class_data['CDELT1']
                            l_deltav= - l_df/ l_class_data['RESTFREQ'] * const.c
                            l_class_data['DELTAV']= l_deltav.value
                            # LOG test
                            #self.m_logger.warn("RESTFREQ {}".format(l_class_data['RESTFREQ']))
                            # Objects Coordinates
                            l_class_data['CDELT2'] = l_chx['scheduled']['ra_offset'].to(unit.deg).value
                            l_class_data['CDELT3'] = l_chx['scheduled']['dec_offset'].to(unit.deg).value

                            l_class_data['AZIMUTH']= getQTableColWithUnit(l_polarization_table,'data_az', unit.deg)
                            l_class_data['ELEVATIO']= getQTableColWithUnit(l_polarization_table,'data_el', unit.deg)
                            l_class_data['CRVAL2']= getQTableColWithUnit(l_polarization_table,'data_ra', unit.deg)
                            l_class_data['CRVAL3']= getQTableColWithUnit(l_polarization_table,'data_dec', unit.deg)
                            # data
                            l_class_data['OBSTIME'] = l_chx['backend']['integration_time']
                            l_class_data['MAXIS1'] = l_chx['backend']['bins']
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
                                except Exception  as e :
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
                            l_splitted_path= os.path.splitext(l_pol_table_path)
                            l_class_path= l_splitted_path[0].replace('.','_')+'_class.fits'
                            self.m_logger.info("input table path: {}".format(l_pol_table_path))
                            self.m_logger.info("class table path: {}".format(l_class_path))
                            try:                                
                                l_class_table.write(l_class_path)
                            except ValueError as e:                                
                                self.m_logger.error("-------- ERROR begin--------")
                                traceback.print_stack()
                                self.m_logger.error("Error writing small classfit table to disk {}".format(e))
                                self.m_logger.error("-------- ERROR end--------")                                
                            classfits.append(l_class_path)
                        except TypeError as e:
                            traceback.print_exc()
                            self.m_logger.error("Error preparing class data: " +str(e))
                        except Exception as e:
                            traceback.print_exc()
                            self.m_logger.error("Error preparing class data: " +str(e))

            # merging classfits data dicts by polarization
            # classfits is a list of classfits data dict
            # classList= defaultdict(list)
            # for d in classfits:
            #     for k,v in d.items():
            #         try:
            #             if not len(v): continue
            #             for n in range(v.shape[0]):
            #                 classList[k].append( v[n] )
            #         except Exception as e:
            #             self.m_logger.error("Exception on appending data classfits data pot: " +str(e))

            # self.m_group_on_off_cal[l_feed]['classlist']= classList


    def classfitsWrite(self, p_on_what):
        """
        Scrittura file con calcolo header
        header prende i dati dai dati generici ricavati dalla scansione scansione
        astropy fits works per column, i have to transpose all data while
        generating new fits cols ( input data are per row basis

        Parameters
        ---------

        p_group: string
            which group, on ,off , cal_on, cal_off
        p_on_what: string
            wich kind of normalized data to use in case of 'on' group
            it assumes 'on', 'on_off', 'cal'

        """
        " for every feed "
        for l_feed in self.m_group_on_off_cal:
            # TODO remove or evaluate
            if  l_feed not in range(0,16):continue
            # Work on ly selected feed
            if self.m_feed:
                if self.m_feed not in str(l_feed):
                    continue
            l_outFileName= self.m_outputPath+ "feed_{}_{}.fits".format(l_feed,p_on_what)
            self.m_logger.info("Preparing classfits file : " + l_outFileName)
            l_newCols=[]
            " for every column expressed in classfits definition .."
            for classCol in self.m_commons.getClassfitsColumnsZip():
                " [ name, form, unit ] column by column data building "
                " fill one column looking into every feed[on], and builds column data "
                l_columnFound = False
                " conditionals, some fields needs dedicated approach"
                l_classList= self.m_group_on_off_cal[l_feed]['classlist']
                " converted fits data matches with classfits columns? "
                " some columns need special care"
                l_inferredCol= classCol[0]
                "  we can choose between on on-off and calibrated spectrum for 'on' group "
                if classCol[0] == "SPECTRUM":
                    l_inferredCol += '_'+ p_on_what.upper()
                " "
                if l_inferredCol in l_classList.keys():
                    " found match, add data to column data "
                    l_columnFound= True
                try:
                    " adding column to classfits if fitszilla representation matches it"
                    if l_columnFound:
                        " some fields needs dedicated approach "
                        if classCol[0] == "SPECTRUM":
                            l_rows= l_classList[l_inferredCol][0].shape[0]
                            l_newCols.append(fits.Column(array= l_classList[l_inferredCol],name= classCol[0],format= "{}D".format(l_rows),unit= classCol[2]))
                        else:
                            l_newCols.append(fits.Column(array= l_classList[l_inferredCol],name= classCol[0],format= classCol[1],unit= classCol[2]) )

                except Exception as e:
                    self.m_logger.error("classfits column creation exception: "+ str(e))
                    self.m_logger.error("column: " +str(classCol))

            l_hdData= self.m_obs_general_data
            " header "
            l_hdu= fits.PrimaryHDU()
            try:
                l_hdu.header['CTYPE1']= "FREQ"
                l_hdu.header['CRVAL1']= 0
                l_hdu.header['CRVAL2']= l_hdData['ra']
                l_hdu.header['CRVAL3']= l_hdData['dec']
                l_hdu.header['OBJECT'] = l_hdData['source']
                l_hdu.header['SOURCE'] = l_hdData['source']
                l_hdu.header['DATE-RED'] = l_hdData['date-red']
                l_hdu.header['LINE'] = l_hdData['line']
                l_hdu.header['CDELT1'] = l_hdData['cdelt1']
                l_hdu.header['RESTFREQ'] = l_hdData['restfreq']
                l_hdu.header['MAXIS1'] = l_hdData['maxis1']
            except KeyError as e:
                self.m_logger.error("Exception filling " + l_outFileName + " header data: "+ str(e))
                " data "
            try:
                " TEST "
                for col in l_newCols:
                    print(col.name + " " + str(col.array.shape))
                l_cdefs= fits.ColDefs(l_newCols)
                l_hdu= fits.BinTableHDU().from_columns(l_cdefs)
            except Exception as e:
                traceback.print_exc()
                self.m_logger.error("Exception creating classfits model file")
                self.m_logger.error("classfits file: " + l_outFileName)
                self.m_logger.error(str(e))
                return
            try:
                if os.path.exists(l_outFileName):
                    os.remove(l_outFileName)
                l_hdu.writeto(l_outFileName)
            except Exception as e:
                self.m_logger.error("Exception writings file")
                self.m_logger.error("classfits file: " + l_outFileName)
                self.m_logger.error(str(e))
