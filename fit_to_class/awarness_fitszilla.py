# -*- coding: utf-8 -*-
"""Fits like parser

    Fits format awarness, it handles fits data toward fits like representation
"""

import numpy as np
import logging
import os
import shutil
import sys
import pdb

from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.io.fits as fits
import astropy.units as unit
from astropy.time import Time
from astropy.table import  QTable, vstack, Column
from memory_profiler import profile

from fit_to_class.fitslike_commons import keywords as kws
from fit_to_class import fitslike_keywords
from fit_to_class import fitslike_commons
from commons import wrapkeys
from fit_to_class import awarness_fitszilla_coord as coord


class Awarness_fitszilla():
    """fitszilla data parser"""

    def __init__(self, l_fitszilla, p_path, p_file_name, p_feed= None):
        """
        Store fitzilla file

        Parameters
        ----------
        l_fitszilla : astropy fits
            Fitszilla to be handled before it gets represented as fitslike

        p_path: string
            path of opened fits

        p_feed: int
            Number of feed to parse (None parse all)

        """
        self.m_errors= [] # errori interni, stringa composta, cchi come
        l_path, self.m_fileName = os.path.split(p_path)
        self.m_path= l_path
        self.m_file_name= p_file_name
        self.m_outputPath= self.m_path
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_jsonRepr = fitslike_keywords.Keyword_json('fitszilla')
        self.m_components = self.m_jsonRepr.fitslike_components()
        self.m_parsingDicts = {}
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        for l_component in self.m_components:
            self.m_parsingDicts[l_component] = \
                self.m_jsonRepr.parser(l_component)
        self.m_fitszilla = l_fitszilla
        self.m_feed= p_feed

    def parse(self):
        """
        Fitszilla parsing

        Parse fitzilla input file following parsing dictionaried
        It only extract keywords from fitzilla without processing them.

        Returns
        -------
        A dictionary = {
            fitslike_key : fitszilla_value
            ... : ...
        }
        """
        self.m_intermediate = {}
        for l_component in self.m_parsingDicts.keys():
            for l_key in self.m_parsingDicts[l_component].keys():
                self.m_intermediate[l_key] = None
                l_keyParserEntry = self.m_parsingDicts[l_component][l_key]
                # Decoding entry
                l_tableName, l_headerOrData, l_inputKeyword  =\
                    fitslike_keywords.Keyword_helper.\
                    decode_entries(l_keyParserEntry)
                if l_tableName is None:
                    self.m_logger.error("empty table name %s, %s, %s",\
                                 l_tableName, l_headerOrData, l_inputKeyword)
                    return
                # todo semplificare la gestione tabella
                l_table = fitslike_keywords.Keyword_helper.\
                    get_table_index(self.m_jsonRepr.input_tables_dict(),
                                     l_tableName)
                if l_table is None:
                    self.m_logger.error("table not found %s", l_tableName)
                    return
                try:
                    l_fitsTable = self.m_fitszilla[l_table]
                except:
                    continue
                # is in header or in data ?
                l_headerOrData = fitslike_keywords.Keyword_helper.\
                    is_header_or_data(l_headerOrData)
                # looking for keyword or data
                l_fitsTableContent = None
                if l_headerOrData is True:
                    l_fitsTableContent = l_fitsTable.header
                else:
                    l_fitsTableContent = l_fitsTable.data
                # fitszilla keyword look up
                try:
                    self.m_intermediate[l_key] = l_fitsTableContent[l_inputKeyword]
                except:
                    self.m_intermediate[l_key]= None
                    #self.m_logger.error("Missing [table, keyword] %s: %s",
                    #                    l_key, l_inputKeyword)
        return self.m_intermediate

    def setOutputPath(self, p_path):
        """set output path to write intermediate processing file"""
        self.m_outputPath= p_path

    def getErrorList(self):
        """
        getter errors
        """
        return self.m_errors

    def process(self):
        """
        Parsed keyword processing

        This class knows how to prepare data to fit fitslike.
        Those data are store into a processed representation.
        Processing order is mandatory, coordinates should be processed after
        spectrum.

        This function creates dictionary:
            ch_X:{
                'frontend': {}
                'backend': {}
                'spectrum': {}
                'coordinates': {}
                'exrtas':{}
            }
        for every channel from fitszilla, hence data are packed on feed-data
        basis.

        _process_spectrum defines :
            'frontend'
            'backend'
            'spectrum'

        _process_coordinates defines:
            'coordinates'

        Returns
        -------
        Processed data representation.

        """
        self.m_processedRepr = {}
        print (self.m_fileName)
        self.m_scheduled= {}
        if self.m_fileName.lower().startswith('sum') :
            self._process_summary()
        else:
            self._process_observation()
            self._process_spectrum()
            self._process_coordinates()
            self._process_extras()
            self._group_data_and_coords()
        return self.m_processedRepr

    def _process_summary(self):
        """
        Keywords from summary.fits
        """        
        l_keys= [ 'geometry', 'scan_type','restfreq', 'backend_name', 'target_ra', 'target_dec', 'velo_def', 'velo_frame', 'velo_rad' ]
        # TODO Aggiustare la lettura delle keyword da file, ora forzate per la mappa
        _geometry= wrapkeys.get_value(self.m_intermediate,'geometry', False)         
        _scan_type= wrapkeys.get_value(self.m_intermediate,'scan_type')         
        l_restFreq= wrapkeys.get_value(self.m_intermediate, 'sum_restfreq') * unit.MHz
        l_target_ra= wrapkeys.get_value(self.m_intermediate,'target_ra') * unit.rad
        l_target_dec= wrapkeys.get_value(self.m_intermediate,'target_dec') * unit.rad        
        if self.m_intermediate['sum_backend_name']== 0.0:
            self.m_intermediate['sum_backend_name']= 'UNKNOWN'
            self._errorFromMissingKeyword('scheduled', 'obs_backend_name')
        _velo_def= wrapkeys.get_value(self.m_intermediate,'velo_def')
        # TODO completare match velo def
        _velo_def_match={
            'RD': 'RADI',
            'OP': 'OPTI',
            'Z': 'RELA',
            'DEFAULT': 'RADI'      
        }        
        if _velo_def in _velo_def_match.keys():
            l_velo_def= _velo_def_match[_velo_def]
        else:
            self.m_logger.warning("VELOCITY DEFINITON NOT HANDLED {}".format(_velo_def))
            l_velo_def= _velo_def_match['DEFAULT']
        if l_velo_def != 'RADI':
            self.m_logger.warning("VELOCITY DEFINITON DIFFERENT FROM \'RADIO\' IS PROBABLY NOT SUPPORTED BY GILDAS {}".format(_velo_def))
        _velo_frame= wrapkeys.get_value(self.m_intermediate,'velo_frame')
        _velo_frame_match={
            'LSRD': 'LSR',
            'LSRK': 'LSR',
            'BARY': 'HEL',
            'OBS': 'OBS',
            'EAR': 'EAR',
            'TOPOCEN': 'TOP',
            'DEFAULT': 'LSR'
        }
        if _velo_frame in _velo_frame_match.keys():
            l_velo_frame= _velo_frame_match[_velo_frame]
        else:
            self.m_logger.warning("VELOCITY FRAME NOT HANDLED {}".format(_velo_frame))
            l_velo_frame= _velo_frame_match['DEFAULT']
        l_velo_rad= wrapkeys.get_value(self.m_intermediate,'velo_rad') * unit.Unit("km/s")
        # Zip dict
        l_values= [_geometry, _scan_type, l_restFreq, self.m_intermediate['sum_backend_name'], l_target_ra, l_target_dec, l_velo_def, l_velo_frame, l_velo_rad]                
        self.m_processedRepr['summary']= dict(zip(l_keys, l_values))        

    def _process_observation(self):
        """
        General observation data review
        Cope with apporpiated phys. units

        feed independent
        """
        try:
            self.m_intermediate['obs_ra'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_ra') * unit.rad
            self.m_intermediate['obs_dec'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_dec')* unit.rad
            self.m_intermediate['obs_ra_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_ra_offset') *unit.rad
            self.m_intermediate['obs_dec_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_dec_offset')* unit.rad
            self.m_intermediate['obs_az_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_az_offset')* unit.rad
            self.m_intermediate['obs_el_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_el_offset')*unit.rad
            self.m_intermediate['obs_gal_lat_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_gal_lat_offset')*unit.rad
            self.m_intermediate['obs_gal_lon_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_gal_lon_offset')*unit.rad                
            self.m_intermediate['obs_user_lat_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_user_lat_offset')*unit.rad
            self.m_intermediate['obs_user_lon_offset'] = \
                wrapkeys.get_value(self.m_intermediate,'obs_user_lon_offset')*unit.rad   
            self.m_intermediate['file_name']= self.m_fileName
            self.m_intermediate['obs_vlsr'] = wrapkeys.get_value(self.m_intermediate,'obs_vlsr')* unit.Unit("km/s")        
        except Exception as e:
            self.m_logger.error(f"{str(e)}")
            _func_name= sys._getframe().f_code.co_name
            raise Exception(f"Mandatory observation keywords cannot have issues!\n Where: {type(self).__name__}::{_func_name}")
        "todo : trasportare le coordinate  per ogni feed?"
        l_scheduled= {}
        try:
            l_scheduled['source']= self.m_intermediate['obs_source']
            l_scheduled['receiver_code']= self.m_intermediate['obs_receiver_code']
            l_scheduled['antenna']= self.m_intermediate['obs_site']
            l_scheduled['date']= self.m_intermediate['obs_date']
            l_scheduled['ra']= self.m_intermediate['obs_ra']
            l_scheduled['dec']= self.m_intermediate['obs_dec']
            l_scheduled['ra_offset']= self.m_intermediate['obs_ra_offset']
            l_scheduled['dec_offset']= self.m_intermediate['obs_dec_offset']
            l_scheduled['az_offset']= self.m_intermediate['obs_az_offset']
            l_scheduled['el_offset']= self.m_intermediate['obs_el_offset']
            l_scheduled['user_lat_offset']= self.m_intermediate['obs_user_lat_offset']
            l_scheduled['user_lon_offset']= self.m_intermediate['obs_user_lon_offset']
            l_scheduled['user_offset_frame']= self.m_intermediate['obs_user_offset_frame']
            l_scheduled['signal']= self.m_intermediate['obs_signal']
            l_scheduled['scan_type']= self.m_intermediate['obs_scantype']
            l_scheduled['file_name']= self.m_intermediate['file_name']
            l_scheduled['vlsr']= self.m_intermediate['obs_vlsr']
        except KeyError as e:
            self.m_logger.error("Key exception : " + str(e))
            self._errorFromMissingKeyword('scheduled', str(e))
        self.m_scheduled= l_scheduled.copy()
        # Check offset frame type
        if self.m_scheduled['user_offset_frame'] == 'GAL':
            self.m_logger.error("Galactic frame conversion not supported! cannot continue");

    def _process_spectrum(self):
        """
        Spectrum keyword processing.
        Fills spectrum final representation

        frontend[backend_id]{
            "fe_feeds"
			"fe_if"
			"fe_polarizations",
			"fe_be_id"
			"fe_frequency"
			"fe_bandwith"
			"fe_local_oscillator"
			"fe_cal_mark_temp"
            "fe_flag_cal"
        }

        backend[frontend_id]{
            "be_id"
			"be_bins"
			"be_sample_rate"
			"be_bandwith"
			"be_frequency"
			"be_integration"
			"be_data_type"
        }

        ch_X: {
            scheduled: scheduled observation data
            frontend: frontend[be_id]
            backend: backend[fe_id]
            spectrum

        }

        Returns
        -------
        None.

        """

        def _add_unit_to_fe(p_feDict):
            """
            Adding units to front end dictionary fields

            Parameters
            ----------
            p_feDict : TYPE
                DESCRIPTION.

            Returns
            -------
            None.

            """
            p_feDict['frequency'] *= unit.MHz
            p_feDict['bandwidth'] *= unit.MHz
            p_feDict['local_oscillator'] *= unit.MHz
            p_feDict['cal_mark_temp'] *= unit.K


        " todo portare fuori el definizioni dei dizionari"
        # Front end dict keys
        l_feDictKeys= [
            'be_id', 'feed', 'if', 'polarizations',
            'frequency', 'bandwidth',
            'local_oscillator', 'cal_mark_temp'
            ]
        # Back end dic keys
        l_beDictKeys= [
            'id', 'bins', 'sample_rate',
            'bandwidth', 'frequency', 'data_type', 'integration_time'
            ]
        # zip front end
        l_frontEnds= {}
        l_zipFrontEnds= {}
        try:
            l_zipFrontEnds = zip(
                        self.m_intermediate['fe_be_id'],
                        self.m_intermediate['fe_feeds'],
                        self.m_intermediate['fe_if'],
                        self.m_intermediate['fe_polarizations'],
                        self.m_intermediate['fe_frequency'],
                        self.m_intermediate['fe_bandwidth'],
                        self.m_intermediate['fe_local_oscillator'],
                        self.m_intermediate['fe_cal_mark_temp'],
                        )
        except Exception as e:
            self.m_logger.error("frontend end zip error: " +str(e))            
        # create dict[backend_id]= front end
        for l_zipFe in l_zipFrontEnds:
            l_feDict= dict(zip(l_feDictKeys, l_zipFe))
            l_feDict['fe_tsys']= 1.0
            # Adding units
            _add_unit_to_fe(l_feDict)
            l_frontEnds[l_feDict['be_id']]= l_feDict.copy()
        """ Pre check backend columns, old MED backend sw lacks of this 2 columns on on section table
           We got bw and freq from rf input table, freq is first spectrum channel if freq and in this case
           needs to be calculated from sky frep and lo
         """
        "@todo rivedere perchè non distinguo correttamente tra lista e None "
        try:
            len(self.m_intermediate['be_frequency'])
        except:
            self.m_intermediate['be_frequency'] =  self.m_intermediate['fe_frequency'] \
                                                - self.m_intermediate['fe_local_oscillator']
            self.m_logger.warning("section table frequency column not found ")
        try:
            len(self.m_intermediate['be_bandwidth'])
        except:
            self.m_intermediate['be_bandwidth'] = self.m_intermediate['fe_bandwidth']
            self.m_logger.warning("section table bandwidth column not found ")
        #  zip backend
        l_backEnds= {}
        l_zipBackend= {}
        try:
            l_zipBackend= zip(
                        self.m_intermediate['be_id'],
                        self.m_intermediate['be_bins'],
                        self.m_intermediate['be_sample_rate'],
                        self.m_intermediate['be_bandwidth'],
                        self.m_intermediate['be_frequency'],
                        self.m_intermediate['be_data_type']
                        )
        except Exception as e:
            self.m_logger.error("back end zip error: " +str(e))        

        # create dict[backend_id]= back end
        for l_zipBe in l_zipBackend:
            l_beDict= dict(zip(l_beDictKeys, l_zipBe))
            l_beDict['integration_time']= self.m_intermediate['be_integration']
            l_beDict['integration_time'] = l_beDict['integration_time'] * unit.ms
            l_beDict['bandwidth'] *= unit.Unit('MHz')
            l_backEnds[l_beDict['id']]= l_beDict.copy()

        # Creates chX_feed_pol: frontend, backend, spectrum
        for l_elBe in l_backEnds.keys():
            l_innerDict= {}
            l_innerDict['scheduled']= self.m_scheduled.copy()
            l_innerDict['frontend']= l_frontEnds[l_elBe]
            l_innerDict['backend']= l_backEnds[l_elBe]
            l_innerDict['spectrum']= {}
            " flag cal data separation before integration "
            l_innerDict['spectrum']['data']={}
            l_innerDict['spectrum']['data']= np.asarray(
                    self.m_intermediate['ch'+str(l_elBe)]
                    )
            l_innerDict['spectrum']['flag_cal']= np.asarray(self.m_intermediate['data_flag_cal'])
            l_feed= l_innerDict['frontend']['feed']
            # Grouping by feed: ch_x
            if l_feed not in self.m_processedRepr.keys():
                self.m_processedRepr[l_feed]={}
            # Work only selected feed or all feeds
            if not self.m_feed:
                self.m_processedRepr[l_feed]['ch_'+str(l_elBe)] = l_innerDict.copy()
            else:
                if self.m_feed in str(l_feed):
                    self.m_processedRepr[l_feed]['ch_'+str(l_elBe)] = l_innerDict.copy()


    def _process_coordinates(self):
        """
        Coordinate data processing

        Every data table entry fitszilla coordinates refers to central feed.

        feed offsets:
            "fe_x_offset"
			"fe_y_offset"

        coordinates : [commons to every data table entry]
            "data_time"
			"data_ra"
			"data_dec"
			"data_az"
			"data_el"
			"data_par_angle"
			"data_derot_angle"

            - Converts radians to arcsec ?
            - Apply offset and rotation angle to every feed in az, el
            - Convert final az, el in ra, dec

        Replicate coordinates for every ch_x (feed only)
        Apply feed offset (az, el)
        Apply derot_angle to every feed (az, el)
        Infer ra, dec for every feed

        Returns
        -------
        None.

        """

        # Process coordinates for every table entries
        l_coordinatesDict= {
            'data_time': np.asarray(self.m_intermediate['data_time']),
            'data_az': np.asarray(self.m_intermediate['data_az']) * unit.rad,
            'data_el': np.asarray(self.m_intermediate['data_el']) * unit.rad,
            'data_derot_angle': np.asarray(
                self.m_intermediate['data_derot_angle']
                )* unit.rad
            }
        # reduce feeds removing duplicate (left/right)
        l_feedCoordinatesDict= dict.fromkeys(self.m_intermediate['fe_feeds'])
        # copy usefull coord. for every feed
        for l_feeds in l_feedCoordinatesDict.keys():
            if self.m_feed:
                if not self.m_feed in str(l_feeds):
                    continue
            l_feedCoordinatesDict[l_feeds]= l_coordinatesDict.copy()
        #pdb.set_trace()
        # Apply offset + rotation adjust for every feed        
        # rest angle for every feed
        l_feedXOffsets = self.m_intermediate['fe_x_offset']* unit.rad
        l_feedYOffsets = self.m_intermediate['fe_y_offset']* unit.rad
        l_feedsRestAngles = coord._feeds_rest_angle(l_feedXOffsets, l_feedYOffsets)
        # for every feed..
        # decor feed - coordinates dict with feed rest angle
        # update observing angle on feed basis
        # offset correction at observation angle
        # calculate ra dec
        # copy the coordinates to processed repr on feed basis
        for l_feed in l_feedCoordinatesDict.keys():
            if self.m_feed:
                if not self.m_feed in str(l_feed):
                    continue
            l_feedCoord = l_feedCoordinatesDict[l_feed]
            l_feedCoord['rest_angle']= \
                l_feedsRestAngles[l_feed]
            l_feedCoord['fe_x_offset']= l_feedXOffsets[l_feed]
            l_feedCoord['fe_y_offset']= l_feedYOffsets[l_feed]
            l_feedObsAngle = coord._observing_angle(
                l_feedsRestAngles[l_feed],
                l_feedCoord['data_derot_angle'])
            # with feed 0 skip correction
            if coord._offset_needs_correction(l_feedCoord):
                l_correctedXoff, l_correctedYoff = coord._offset_corrections(
                    l_feedObsAngle,
                    l_feedCoord['fe_x_offset'],
                    l_feedCoord['fe_y_offset']
                    )
            else:
                l_correctedXoff = l_feedCoord['fe_x_offset']
                l_correctedYoff = l_feedCoord['fe_y_offset']
            l_obstime = Time(l_feedCoord['data_time'] * unit.day,
                             format= 'mjd',
                             scale= 'utc')
            l_feedCoord['time_mjd'] = np.asarray(l_obstime)
            l_location = self.m_commons.get_site_location(
                self.m_intermediate['site'].lower()
                )
            # Final coordinates ra, dec after calculations
            l_feedCoord['data_ra'], l_feedCoord['data_dec'] = \
                coord._azel_to_radec(self.m_scheduled,
                                    l_obstime,
                                    l_feedCoord['data_el'],
                                    l_feedCoord['data_az'],
                                    l_correctedXoff,
                                    l_correctedYoff,
                                    l_location)
            # coordinates dict storage
            for feed in self.m_processedRepr.keys():
                if self.m_feed:
                    if not self.m_feed in str(feed):
                        continue
                for chx in self.m_processedRepr[feed]:
                    _section = self.m_processedRepr[feed][chx]
                    if _section['frontend']['feed'] == l_feed:
                        _section['coordinates']= l_feedCoord.copy()

    def _process_extras(self):
        """
        Extra data processing
        Filling not present keywords
        """
        for feed in self.m_processedRepr.keys():
            if self.m_feed:
                if not self.m_feed in str(feed):
                    continue
            for chx in self.m_processedRepr[feed]:
                _section = self.m_processedRepr[feed][chx]
                _section['extras']= {}
                " weather "
                l_weather= self.m_intermediate['ex_weather']
                _section['extras']['weather']= []
                for el in l_weather:
                    _section['extras']['weather'].append(\
                        self.m_commons.calculate_weather(el[1] + 273.15, el[0]))



    def _group_data_and_coords(self):
        """
        Groups data table entry data, those data are composed by multiple entries
        that needs to be reduced / group for example by flag_cal
        Grouped data are arranged similarly to fitszilla data table row entries
        After grouping data by flag cal spectrum data field is deleted

        Applying stokes differentiation
        Adding polarization infos

        At the end of operation data table has integrated data and grouped table

        Stage 2:

            Check on off cal..
            grouping by on off cal also

        @todo Attenzione, non riesco ad aggregare la tabella se inserisco i dati cone unit !!
        """

        def _is_cal(p_subscan, p_group):
            """
            Genera il flag calibrazione attiva

            Parameters
            ----------
            p_sub : astropy table group
                group with flag_cal
            p_subscan : subscan with meta data
                group with flag_cal

            Returns
            -------
            true calibrazione attiva
            """
            l_signal= p_subscan['scheduled']['signal']
            if l_signal in kws['keys_cal_on']:
                return True
            # TODO Check this part..
            #elif np.any(p_group['flag_cal']):
            #    return True
            # TODO Non chiaro...
            #if p_group['flag_cal'] > 0:
            #    return True
            return False

        def _is_on(p_subscan, p_feed):
            """
            @TODO Check according to appropriate keyword
            The main signal keyword refers to central feed (feed 0)
            if an appropriate feed reference is not provided through
            a dedicated fitszilla keyword.            

            Parameters
            ----------
            p_subscan : dict
                subscan ch_x
            p_feed : subscan with meta data
                group with flag_cal

            Returns
            -------
            None.
            """                        
            if 'scheduled' not in p_subscan:
                raise Exception("Cannot continue without mandatory keywords [scheduled]")
            # Signal value
            if 'signal' not in p_subscan['scheduled']:                
                raise Exception("Cannot continue without mandatory keywords [scheduled][signal]")
            _signal_value= p_subscan['scheduled']['signal']
            if _signal_value == None :
                self.m_logger.warning("Signal keyword not present using offset to set is_on flag!")
                try:
                    _az_offset= p_subscan['ccordinates']['az_offset']
                    return _az_offset > 1e-4 * unit.rad        
                except KeyError as e:
                    raise Exception("Cannot continue without mandatory keywords ['ccordinates']['az_offset']")
            # Reference feed, this keyword is not mandatory
            _reference_feed= 0
            if 'reference_feed' not in p_subscan['scheduled']:
                self.m_logger.warn("Reference feed number relating to SIGNAL keyword is missing [scheduled][reference_feed]")
            else:
                _reference_feed= int(self.p_subscan['scheduled']['reference_feed'])            
            #  SIGNAL keywords checking
            if _signal_value in kws['keys_on']:
                # Signal key is related to 'On source' context
                if p_feed == 0:
                    return True
                else:
                    return False
            else:
                # Signal key is NOT related to 'On source' context
                if p_feed == 0:
                    return False
                else:
                    return True            

        # Start function

        for feed in self.m_processedRepr.keys():
            # Work only selected feed
            if self.m_feed:
                if not self.m_feed in str(feed):
                    continue
            for chx in self.m_processedRepr[feed]:
                _section= self.m_processedRepr[feed][chx]
                l_coo= _section['coordinates']
                l_spec= _section['spectrum']
                # Check, some scan split data feed one per file
                # Expecting numpy ndarray for spectrum data
                if 'data' not in l_spec.keys():
                    continue
                shape = l_spec['data'].shape
                if not shape : continue
                if shape[0]== 0 : continue
                # Stokes ?
                if _section['backend']['data_type']== 'stokes':
                    " split stokes, repeating data "
                    data= l_spec['data']
                    l_bins= _section['backend']['bins']
                    self.m_logger.info("section {} : data shape {} - bins {}".format(chx, shape, l_bins))
                    try:
                        L= (data[:,:l_bins], "LL")
                        R= (data[:, l_bins: 2*l_bins], "RR")
                        Q= (data[:, 2*l_bins: 3*l_bins],"LR")
                        U= (data[:, -l_bins: ], "RL")
                        self.m_logger.info("shapes: L {} - R {} - Q {} - U {}".format(L[0].shape, R[0].shape, Q[0].shape, U[0].shape))
                    except Exception as e:
                        self.m_logger.error("Splitting stokes data by pol raised an excp: {}".format(e))
                        self.m_logger.error("section {} : data {} -  bins {}".format(chx, data,l_bins))
                        break
                    l_tGroupPol=[]
                    for obj in [L, R, Q, U]:
                        poldata= obj[0].astype(np.int32)
                        polLabel= obj[1]
                        polLabel= np.full((poldata.shape[0],), polLabel)
                        l_pol_table= QTable()
                        try:
                            l_keys= ["data_time_mjd", "pol", "data_az",
                                    "data_el", "data_derot_anngle", "data_ra",
                                    "data_dec", "weather", "data", "flag_cal"]
                            l_data= [l_coo["data_time"], polLabel, l_coo["data_az"],
                                     l_coo["data_el"],l_coo["data_derot_angle"],l_coo["data_ra"],
                                     l_coo["data_dec"], _section['extras']['weather'], poldata,
                                     l_spec["flag_cal"]]
                            for n, v in zip(l_keys, l_data):
                                # If i have scalar without quantity
                                try:
                                    name_unit= n +"_u_" + str(v.unit)
                                    col=  Column(v.value, name= name_unit )
                                except:
                                    col=  Column(v, name= n )
                                l_pol_table.add_column(col)

                        except Exception as e:
                            self.m_logger.error("Exception creating data tables for stokes data : {}".format(e) )                            
                            continue
                        # Add this pol table to table list
                        l_tGroupPol.append(l_pol_table)
                    # Group and aggregation
                    if not l_tGroupPol:
                        _section['groups']= QTable()
                        continue                                      
                    # Regroup data from varius pol into one table
                    l_oneTable= vstack(l_tGroupPol)
                    
                else: # SPECTRUM OR SINGLE POL
                    " manage pwr spectrum "
                    l_oneTable = QTable()
                    l_shape= l_coo["time_mjd"].shape
                    l_pol= _section['frontend']['polarizations']
                    " uniform polarization strings "
                    if l_pol== 'LCP': l_pol= 'LL'
                    if l_pol== 'RCP': l_pol= 'RR'
                    l_polShaped= np.full(l_shape, l_pol)
                    l_keys= [l_coo["data_time"],l_polShaped,l_coo["data_az"],
                             l_coo["data_el"],l_coo["data_derot_angle"], l_coo["data_ra"],
                             l_coo["data_dec"], _section['extras']['weather'], l_spec['data'],
                             l_spec["flag_cal"]]
                    l_names= ["data_time_mjd", "pol", "data_az",
                            "data_el", "data_derot_angle", "data_ra",
                            "data_dec", "weather", "data", "flag_cal"]
                    for n, v in zip(l_names, l_keys):
                        # If i have scalar without quantity
                        try:
                            name_unit= n +"_u_" + str(v.unit)
                            col=  Column(v.value, name= name_unit )
                        except:
                            col=  Column(v, name= n )
                        l_oneTable.add_column(col)                                     
                
                # Remove data already present in groups
                del l_coo["time_mjd"]
                del l_coo["data_time"]
                del l_coo["data_az"]
                del l_coo["data_el"]
                del l_coo["data_derot_angle"]
                del l_coo["data_ra"]
                del l_coo["data_dec"]
                del _section['extras']['weather']
                del l_spec["data"]
                del l_spec["flag_cal"]                              

                # One table to rule'em all   
                
                # Adding ON OFF column to the whole table
                _isOn= _is_on(_section, feed)
                if _isOn:
                    on_col_data= [1]*len(l_oneTable)
                else:
                    on_col_data= [0]*len(l_oneTable)
                on_col= Column(on_col_data, 'is_on')
                l_oneTable.add_column(on_col)                
                # Adding on off cal column
                l_temporary_tables= []                
            
                for group in  l_oneTable.group_by(['pol']).groups:
                    # Cal on cal off
                    _isCal= _is_cal(_section, group)
                    #Cal on
                    if _isCal and np.any(group['is_on']):
                        cal_on_col= Column([1]*len(group),'cal_on')
                        group.add_column(cal_on_col)
                        self.m_logger.info("{}--{}_{}_{} is cal_on".format(self.m_file_name, feed, chx, group['pol'][0]))
                    #Cal off
                    if _isCal and not np.any(group['is_on']):
                        cal_off_col= Column([1]*len(group),'cal_off')
                        group.add_column(cal_off_col)
                        self.m_logger.info("{}--{}_{}_{} is cal_off".format(self.m_file_name, feed, chx, group['pol'][0]))
                    # Signal
                    if not _isCal and np.any(group['is_on']):
                        signal_col= Column([1]*len(group),'signal')
                        group.add_column(signal_col)           
                        self.m_logger.info("{}--{}_{}_{} is signal".format(self.m_file_name, feed, chx, group['pol'][0]))
                    # Off
                    if not _isCal and not np.any(group['is_on']):
                        off_col= Column([1]*len(group),'reference')
                        group.add_column(off_col)
                        self.m_logger.info("{}--{}_{}_{} is reference".format(self.m_file_name, feed, chx, group['pol'][0]))
                    # Removing unusefull data from table group 
                    del group['is_on']
                    del group['flag_cal']
                    # Save group with new columns
                    l_temporary_tables.append(group)
                # Stack and group with on/off and cal on/cal off
                l_stacked= vstack(l_temporary_tables)
                _section['groups_by_pol']= l_stacked.group_by(['pol']).groups
                # Write groups to disk
                l_pols={}
                for group in _section['groups_by_pol']:
                    # Write table to disk
                    if 'pol' in group.colnames:
                        l_pols[group['pol'][0]]= self._writeTableToFile(group, feed, chx, group['pol'][0])
                    else:
                        self.m_logger.error("Missing polarization label! skip table")
                # Keep trace of on disk tables
                _section['pol_tables_dict']= l_pols                
                # Remove unesfull data on fits representation
                del _section['groups_by_pol']                    
                    


    def _writeTableToFile(self, p_table, p_feed, p_section, p_pol) -> str:
        """
        Writes table group to disk grouped by feed, section, pol
        Doesn't overwrite output folder, it creates the folder when missing
        (one folder per feed)

        Parameters
        ----------
        p_table : QTable
            table to write on disk.
        p_feed : int
            feed number
        p_section : string
            section number
        p_pol : string
            polarization
        
        Returns:
        --------
        table on disk path
        """
        try:
            l_fname= os.path.join(self.m_outputPath, "fits")
            # Create fits/feed_{}
            if not os.path.exists(l_fname):
                os.mkdir(l_fname)
            l_fname= os.path.join(l_fname,'feed_{}'.format(p_feed))
            if not os.path.exists(l_fname):             
                os.mkdir(l_fname)
            # If file il is already present remove it
            l_fname= os.path.join(l_fname, "{}_{}_{}.fits".format(p_section, p_pol, self.m_fileName))
            if os.path.exists(l_fname):
                os.remove(l_fname)
            p_table.write(l_fname)
            return l_fname
        except IOError as e:
            self.m_logger.error("IO Exception writing intermediate table to disk: \n {} \n {}".format(l_fname, e))
            return ''
        except Exception as e:
            self.m_logger.error("Exception writing intermediate table to disk: \n {} \n {}".format(l_fname, e))
            return ''


    def _errorFromMissingKeyword(self, p_section, p_key):
        """
        Set  error for process upon missing keyws from fitszilla

        Parameters
        ----------
        p_key : string
            keyword section
        p_key : string
            missing keyword

        Returns
        -------
        None.
        """
        l_mandatory={}
        l_mandatory['scheduled']=['obs_backend_name','obs_site']
        if p_key in l_mandatory[p_section]:
            self.m_errors.append("key error: {}:{}".format(p_section, p_key))

