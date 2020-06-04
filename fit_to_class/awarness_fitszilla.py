# -*- coding: utf-8 -*-
"""Fits like parser

    Fits format awarness, it handles fits data toward fits like representation
"""
from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.io.fits as fits
import astropy.units as unit
from astropy.time import Time
from astropy.table import  QTable, vstack, Column
import numpy as np
import copy
import logging
import os
import traceback
import pdb
from memory_profiler import profile

from fit_to_class import fitslike_keywords
from fit_to_class import fitslike_commons


class Awarness_fitszilla():
    """fitszilla data parser"""

    def __init__(self, l_fitszilla, p_path):
        """
        Store fitzilla file

        Parameters
        ----------
        l_fitszilla : astropy fits
        Fitszilla to be handled before it gets represented as fitslike

        p_path: string
        path of opened fits

        Returns
        -------
        None
        """
        self.m_errors= [] # errori interni, stringa composta, cchi come
        l_path, self.m_fileName = os.path.split(p_path)        
        self.m_commons = fitslike_commons.Fitslike_commons() 
        self.m_jsonRepr = fitslike_keywords.Keyword_json('fitszilla')
        self.m_components = self.m_jsonRepr.fitslike_components()        
        self.m_parsingDicts = {}
        self.m_logger = logging.getLogger(self.m_commons.logger_name())        
        for l_component in self.m_components:            
            self.m_parsingDicts[l_component] = \
                self.m_jsonRepr.parser(l_component)        
        self.m_fitszilla = l_fitszilla        

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
        self.m_scheduled= {}
        if 'summary' in self.m_fileName :
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
        l_keys= [ 'restfreq', 'backend_name' ]        
        l_restFreq= self.m_intermediate['sum_restfreq'] * unit.MHz
        if self.m_intermediate['sum_backend_name']== 0.0:
            self.m_intermediate['sum_backend_name']= 'UNKNOWN'
            self._errorFromMissingKeyword('scheduled', 'obs_backend_name')
        l_values= [l_restFreq, self.m_intermediate['sum_backend_name']]        
        self.m_processedRepr['summary']= dict(zip(l_keys, l_values))
        print(self.m_processedRepr['summary'])
            
    def _process_observation(self):
        """
        General observation data review
        Cope with apporpiated phys. units
        """        
        self.m_intermediate['obs_ra'] = \
            self.m_intermediate['obs_ra'] * unit.rad
        self.m_intermediate['obs_dec'] = \
            self.m_intermediate['obs_dec']* unit.rad
        self.m_intermediate['obs_ra_offset'] = \
            self.m_intermediate['obs_ra_offset'] *unit.rad
        self.m_intermediate['obs_dec_offset'] = \
            self.m_intermediate['obs_dec_offset']* unit.rad
        self.m_intermediate['obs_az_offset'] = \
            self.m_intermediate['obs_az_offset']* unit.rad
        self.m_intermediate['obs_el_offset'] = \
            self.m_intermediate['obs_el_offset']*unit.rad
        self.m_intermediate['file_name']= self.m_fileName
        self.m_intermediate['obs_vlsr'] *=  unit.Unit("km/s")
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
            l_scheduled['signal']= self.m_intermediate['obs_signal']
            l_scheduled['scan_type']= self.m_intermediate['obs_scantype']
            l_scheduled['file_name']= self.m_intermediate['file_name']            
            l_scheduled['vlsr']= self.m_intermediate['obs_vlsr']
        except KeyError as e:
            self.m_logger.error("Key exception : " + str(e))
            self._errorFromMissingKeyword('scheduled', str(e))
        self.m_scheduled= l_scheduled.copy()
       
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
            traceback.print_exc()
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
            traceback.print_exc()                    
        
        # create dict[backend_id]= back end
        for l_zipBe in l_zipBackend:
            l_beDict= dict(zip(l_beDictKeys, l_zipBe))
            l_beDict['integration_time']= self.m_intermediate['be_integration']            
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
            " Grouping by feed: ch_x"
            if l_feed not in self.m_processedRepr.keys():
                self.m_processedRepr[l_feed]={}                
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
        def _coordinate_feeds_rest_angle(p_xoffsets, p_yoffsets):
            """
            Rest angle generation for every feed:
                - derotaion angle in resting condition                    

            Parameters
            ----------
            p_xoffsets : list
                feeds X offset list
            p_yoffsets : list
                feeds y offset list

            Returns
            -------
            todo

            """            
            if (len(p_xoffsets)) <= 2:
                return np.array([0]*len(p_xoffsets))
            l_npXOffsets= np.asarray(p_xoffsets)
            l_npYOffsets= np.asarray(p_yoffsets)
            l_num_lat_feeds= len(p_xoffsets) -1
            l_angle_range= np.arange(1, 0, -1/l_num_lat_feeds)
            l_rest_angle_def = l_angle_range * 2 * np.pi * unit.rad
            l_w0= np.where((l_npXOffsets[1:] > 0) & (l_npYOffsets[1:] == 0.))[0][0]
            return np.concatenate(([0], np.roll(l_rest_angle_def.to(unit.rad).value, l_w0))) * unit.rad
    
        def _coordinates_observing_angle(p_rest_angle, p_derot_angle):
            """
            Observign angle calculation for one feed

            Parameters
            ----------
            p_rest_angle : double, rad
                feed rest angle.
            p_derot_angle : array, rad
                actual derotation angle 

            Returns
            -------
            Feed observing angle
            rad

            """
            if not hasattr(p_rest_angle, 'unit'):
                p_rest_angle *= unit.rad
            if not hasattr(p_derot_angle, 'unit'):
                p_derot_angle *= unit.rad
            #pdb.set_trace()
            return p_rest_angle + (2 * np.pi * unit.rad - p_derot_angle)

        def _coordinates_offset_needs_correction(p_coord):
            """
            Check this feed offset amount.
            Nearly no offset feeds don't need correction

            Parameters
            ----------
            p_coord : Dictionary
                Coordinate infos for on feed         

            Returns
            -------
            True needs correction

            """               
            if np.abs(p_coord['fe_x_offset'] ) < \
                np.radians(0.001 / 60.) * unit.rad and \
                np.abs(p_coord['fe_y_offset'] ) < \
                    np.radians(0.001 / 60.)* unit.rad :
                    return False
            return True
        
        def _coordinates_offset_corrections(p_obs_angle, p_xoffset, p_yoffset):
            """
            Feed offset correction at observation angle

            Parameters
            ----------
            p_obs_angle : float rad
                Observation angle
            p_xoffset : float rad
                feed x offset
            p_xoffset : float rad
                feed y offset

            Returns
            -------
            Corrected feed offset float rad.

            """
            l_sep = np.sqrt(p_xoffset**2. + p_yoffset**2.)
            l_corr_xoff = l_sep * np.cos(p_obs_angle)
            l_corr_yoff = l_sep * np.sin(p_obs_angle)
            return l_corr_xoff, l_corr_yoff
            
        def _coordinates_azel_to_radec(p_obstimes,
                                       p_el, p_az, 
                                       p_xoffs, p_yoffs,
                                       p_location):
            """
            Converion from alt-az to ra-dec.
            Offset must be correccted based on observation time.
            
            Returns:
            --------
            Actual ra dec lists
            """
            # Calculate observing angle            
            p_yoffs = p_yoffs
            p_xoffs = p_xoffs
            l_el = copy.deepcopy(p_el)
            l_az = copy.deepcopy(p_az)
            l_el += p_yoffs
            l_az += p_xoffs / np.cos(l_el)
            l_coordsAltAz = AltAz(az=Angle(l_az), 
                                  alt=Angle(l_el),
                                  location= p_location,
                                  obstime= p_obstimes)            
            # According to line_profiler, coords.icrs is *by far* the longest
            # operation in this function, taking between 80 and 90% of the
            # execution time. Need to study a way to avoid this.
            l_coords_deg = l_coordsAltAz.transform_to(ICRS)
            l_ra = np.radians(l_coords_deg.ra)
            l_dec = np.radians(l_coords_deg.dec)
            
            return l_ra, l_dec        

        # Process coordinates for every table entries    
        l_coordinatesDict= {
            'data_time': np.asarray(self.m_intermediate['data_time']),
            'data_az': np.asarray(self.m_intermediate['data_az']) \
                * unit.rad,
            'data_el': np.asarray(self.m_intermediate['data_el']) \
                * unit.rad,
            'data_derot_angle': np.asarray(
                self.m_intermediate['data_derot_angle']
                )* unit.rad
            }
        # reduce feeds removing duplicate (left/right)
        l_feedCoordinatesDict= dict.fromkeys(self.m_intermediate['fe_feeds'])
        # copy usefull coord. for every feed
        for l_feeds in l_feedCoordinatesDict.keys():
            l_feedCoordinatesDict[l_feeds]= l_coordinatesDict.copy()
        #pdb.set_trace()
        # Apply offset + rotation adjust for every feed
        # rest angle for every feed            
        l_feedXOffsets = self.m_intermediate['fe_x_offset']* unit.rad
        l_feedYOffsets = self.m_intermediate['fe_y_offset']* unit.rad
        l_feedsRestAngles = _coordinate_feeds_rest_angle(l_feedXOffsets, l_feedYOffsets)
        # for every feed..
        # decor feed - coordinates dict with feed rest angle        
        # update observing angle on feed basis
        # offset correction at observation angle
        # calculate ra dec
        # copy the coordinates to processed repr on feed basis        
        for l_feed in l_feedCoordinatesDict.keys():
            l_feedCoord = l_feedCoordinatesDict[l_feed]
            l_feedCoord['rest_angle']= \
                l_feedsRestAngles[l_feed]
            l_feedCoord['fe_x_offset']= l_feedXOffsets[l_feed]
            l_feedCoord['fe_y_offset']= l_feedYOffsets[l_feed]
            l_feedObsAngle = _coordinates_observing_angle(
                l_feedsRestAngles[l_feed],
                l_feedCoord['data_derot_angle'])
            # with feed 0 skip correction
            if _coordinates_offset_needs_correction(l_feedCoord):                
                l_correctedXoff, l_correctedYoff = _coordinates_offset_corrections(
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
                _coordinates_azel_to_radec(l_obstime,
                                           l_feedCoord['data_el'],
                                           l_feedCoord['data_az'],
                                           l_correctedXoff,
                                           l_correctedYoff,
                                           l_location)  
                
            # coordinates dict storage
            for feed in self.m_processedRepr.keys():
                for chx in self.m_processedRepr[feed]:
                    l_chx = self.m_processedRepr[feed][chx]
                    if l_chx['frontend']['feed'] == l_feed:
                        l_chx['coordinates']= l_feedCoord.copy()
     
    def _process_extras(self):
        """
        Extra data processing
        Filling not present keywords    
        """                        
        for feed in self.m_processedRepr.keys():
            for chx in self.m_processedRepr[feed]:
                l_chx = self.m_processedRepr[feed][chx]
                l_chx['extras']= {}
                " weather "
                l_weather= self.m_intermediate['ex_weather']
                l_chx['extras']['weather']= []
                for el in l_weather:
                    l_chx['extras']['weather'].append(\
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
        
        @todo Attenzione, non riesco ad aggregare la tabella se inserisco i dati cone unit !!
        """
        for feed in self.m_processedRepr.keys():
                for chx in self.m_processedRepr[feed]:
                    l_chx= self.m_processedRepr[feed][chx]                    
                    l_coo= l_chx['coordinates']
                    l_spec= l_chx['spectrum']                                       
                    if l_chx['backend']['data_type']== 'stokes':
                        " split stokes, repeating data "
                        data= l_spec['data']
                        l_bins= l_chx['backend']['bins']                        
                        L= (data[:,:l_bins], "LL")
                        R= (data[:, l_bins: 2*l_bins], "RR")
                        Q= (data[:, 2*l_bins: 3*l_bins],"LR")
                        U= (data[:, -l_bins: ], "RL")
                        l_tGroup=[]                        
                        for obj in [L, R, Q, U]:
                            poldata= obj[0].astype(np.int32)
                            polLabel= obj[1]
                            polLabel= np.full((poldata.shape[0],), polLabel)                            
                            l_oneTable= QTable()                            
                            try:                                
                                l_keys= [l_coo["data_time"], polLabel, l_coo["data_az"],
                                         l_coo["data_el"],l_coo["data_derot_angle"],l_coo["data_ra"],
                                         l_coo["data_dec"], l_chx['extras']['weather'], poldata,
                                         l_spec["flag_cal"]]
                                l_names= ["data_time_mjd", "pol", "data_az",
                                        "data_el", "data_derot_anngle", "data_ra",
                                        "data_dec", "weather", "data", "flag_cal"]
                                for n, v in zip(l_names, l_keys):
                                    # If i have scalar without quantity 
                                    try:                                
                                        name_unit= n +"_u_" + str(v.unit)
                                        col=  Column(v.value, name= name_unit )                                        
                                    except:                                                    
                                        col=  Column(v, name= n )
                                    l_oneTable.add_column(col)  
                                
                            except Exception as e:
                                traceback.print_exc()
                                self.m_logger.error("Exception creating data tables fro stokes pols : " +str(e) )
                                pdb.set_trace()
                                continue                            
                            " add this pol table to table list "
                            l_tGroup.append(l_oneTable)
                        " group and aggregation " 
                        if l_tGroup:
                            l_oneTable= vstack(l_tGroup)                            
                            l_oneTable= l_oneTable.group_by(['pol', 'flag_cal'])                            
                            l_oneTableAggregated= l_oneTable.groups.aggregate(np.mean)
                            " @todo for an unkown reason i cannot aggregate data column..doing it manually "                            
                            if 'data' not in l_oneTableAggregated.colnames:
                                self.m_logger.warning("Aggretating data column manually")
                                l_spectrum= []
                                for gr in l_oneTable.groups:
                                    l_data= gr['data']
                                    l_spectrum.append(np.mean(l_data, axis= 0))
                                l_dataCol= Column(l_spectrum, name= 'data')
                                l_oneTableAggregated.add_column(l_dataCol)
                                l_oneTableAggregated= l_oneTableAggregated.group_by(['pol', 'flag_cal'])
                            l_chx['groups']= l_oneTableAggregated
                            
                        else:
                            l_chx['groups']= QTable()
                    else: # SPECTRUM OR SINGLE POL 
                        " manage pwr spectrum "
                        l_oneTable = QTable()
                        l_shape= l_coo["time_mjd"].shape        
                        l_pol= l_chx['frontend']['polarizations']                                        
                        " uniform polarization strings "
                        if l_pol== 'LCP': l_pol= 'LL'
                        if l_pol== 'RCP': l_pol= 'RR'    
                        l_polShaped= np.full(l_shape, l_pol)                        
                        l_keys= [l_coo["data_time"],l_polShaped,l_coo["data_az"],
                                 l_coo["data_el"],l_coo["data_derot_angle"], l_coo["data_ra"],
                                 l_coo["data_dec"], l_chx['extras']['weather'], l_spec['data'],
                                 l_spec["flag_cal"]]
                        l_names= ["data_time_mjd", "pol", "data_az",
                                "data_el", "data_derot_angle", "data_ra",
                                "data_dec", "weather", "data", "flag_cal"]
                        for n, v in zip(l_names, l_keys):
                            " If i have scalar without quantity "
                            try:                                
                                name_unit= n +"_u_" + str(v.unit)
                                col=  Column(v.value, name= name_unit )
                            except:                                                    
                                col=  Column(v, name= n )
                            l_oneTable.add_column(col)  
                        " group and aggregate "
                        l_integrationTime = len(l_oneTable) * l_chx['backend']['integration_time']
                        l_chx['backend']['integration_time']= l_integrationTime
                        l_oneTable= l_oneTable.group_by(['pol','flag_cal'])             
                        l_oneTableAggregated= l_oneTable.groups.aggregate(np.mean)
                        " @todo for an unkown reason i cannot aggregate data column..doing it manually "
                        if 'data' not in l_oneTableAggregated.colnames:
                            self.m_logger.warning("Aggretating data column manually")
                            l_spectrum= []
                            for gr in l_oneTable.groups:
                                l_data= gr['data']
                                l_spectrum.append(np.mean(l_data, axis= 0))
                            l_dataCol= Column(l_spectrum, name= 'data')
                            l_oneTableAggregated.add_column(l_dataCol)
                            l_oneTableAggregated= l_oneTableAggregated.group_by(['pol', 'flag_cal'])
                        l_chx['groups']= l_oneTableAggregated
                                                
                    " remove data already present in groups "                        
                    del l_coo["time_mjd"]
                    del l_coo["data_time"]
                    del l_coo["data_az"]
                    del l_coo["data_el"]
                    del l_coo["data_derot_angle"]
                    del l_coo["data_ra"]
                    del l_coo["data_dec"]
                    del l_chx['extras']['weather']
                    del l_spec["data"]
                    del l_spec["flag_cal"]                     
                                
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
            
    