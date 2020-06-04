# -*- coding: utf-8 -*-
"""Fitslike commons helper"""


from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.units as unit
from astropy.io import fits
import numpy as np
import re
import pdb

keywords= {
    "key_on":"SIGNAL",
    "key_off": "REFERENCE",
    "key_off_cal": "REFCAL",
    "key_on_cal": "REFSIG",
    "key_sig_cal": "SIGCAL",
    "keys_on": ["SIGNAL", "SIGCAL", "REFSIG"],
    "keys_off": ["REFCAL", "REFERENCE"],
    "keys_cal_on": ["REFCAL", "SIGCAL", "REFSIG"],
    "keys_cal_off": ["REFERENCE", "SIGNAL"]
    
}

class Fitslike_commons():
    """Helper class"""

    @staticmethod
    def logger_name():
        """
        Getter logger name

        Returns
        -------
        String with logger name.

        """
        return 'root'
    
    @staticmethod
    def build_attribute_dict(p_lstAttributes):
        """Build fitlslike component attribute dictionary
        Parameters
        ----------
        p_lstAttributes : string list
            List with dictionary entries

        Returns
        -------
        Attribute dictionary with 'value' and  'type' entries
        """
        l_outDict = {}
        for l_attribute in p_lstAttributes:
            l_outDict[l_attribute] = {'value': '', 'type': ''}
        return l_outDict.copy()
    
    @staticmethod
    def dump_attribute(p_key, p_attributeValue):
        """fitslike component value dumping
        component's dict :
            {'key' : {'value' :'', 'type':''}}
        It prints key and value type part
        """
        print(p_key + " " + str(p_attributeValue))
        
    @staticmethod
    def get_site_location(p_site):
        """
        Getter site earth location
        """
        locations = {
            'srt': EarthLocation(4865182.7660, 791922.6890, 4035137.1740, unit=unit.m),
             'medicina': EarthLocation(Angle("11:38:49", unit.deg),
                                       Angle("44:31:15", unit.deg),
                                       25 * unit.meter),
             'greenwich': EarthLocation(lat=51.477*unit.deg, lon=0*unit.deg)}
        return locations[p_site.lower()]

    @staticmethod
    def filter_input_files(p_type):
        """
        Input file extension filter based         
        
        Parameters
        ----------
        p_type : string
            input data file type

        Returns
        -------
        None.

        """    
        l_inputDict={
            'fitszilla': r"[\s\S]*\.fits"
            }
        try:
            return l_inputDict[p_type]
        except KeyError:
            return ''
        return ''
    
    @staticmethod
    def telescope_label_from_channel_name(p_chName):
        """
        Generation telescope name for feed
        """            
        
        def interpret_chan_name(chan_name):
            """Get feed, polarization and baseband info from chan name.
        
            Examples
            >>> feed, polar, baseband = interpret_chan_name('blablabal')
            >>> feed  # None
            >>> polar  # None
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Ch0')
            >>> feed
            0
            >>> polar  # None
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Feed1_LCP')
            >>> feed
            1
            >>> polar
            'LCP'
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Feed2_LCP_3')
            >>> feed
            2
            >>> polar
            'LCP'
            >>> baseband
            3
            """
            chan_re = re.compile(r'^Ch([0-9]+)$'
                         r'|^Feed([0-9]+)_([a-zA-Z]+)$'
                         r'|^Feed([0-9]+)_([a-zA-Z]+)_([0-9]+)$')
            
            matchobj = chan_re.match(chan_name)
            if not matchobj:
                return None, None, None
        
            matches = [matchobj.group(i) for i in range(7)]
            polar, baseband = None, None
            if matches[6] is not None:
                baseband = int(matchobj.group(6))
                polar = matchobj.group(5)
                feed = int(matchobj.group(4))
            elif matches[3] is not None:
                polar = matchobj.group(3)
                feed = int(matchobj.group(2))
            else:
                feed = int(matchobj.group(1))
    
            return feed, polar, baseband
        
        _, polar, _ = interpret_chan_name(p_chName)
    
        if polar.startswith('L'):
            return 'LL'
        elif polar.startswith('R'):
            return 'RR'
        elif polar.startswith('Q'):
            return 'LR'
        elif polar.startswith('U'):
            return 'RL'
        else:
            raise ValueError('Unrecognized polarization')
    
    @staticmethod 
    def class_telescope_name(p_section, p_bk, p_pol, p_section_name):
        """
        Generazione strigna telescope classfits        
        """                         
        try:                
            return '{}-{}-{}-{}'.format(
                        p_section['scheduled']['antenna'][:3],
                        p_section['scheduled']['receiver_code'][0],
                        p_bk[:3],                        
                        p_pol
                        )
        except Exception as e:
            raise ValueError("can't build TELESCOP string [err]: " + str(e))

    @staticmethod
    def calculate_weather(p_tmp, p_u):
        """Get the meters of H2O, using the formula from old converter.
        Unsure if this is correct in all cases"""
        RS = 8.314472
        mw = 0.018015
        md = 0.0289644
        eps0 = mw / md
        k1 = 77.60
        k2 = 70.4
        k3 = 3.739E5    
        p_tmp = p_tmp - 273.15
        H = (np.log10(p_u) - 2.0) / 0.4343 + (17.62 * p_tmp) / (243.12 + p_tmp)    
        DPT = 243.12 * H / (17.62 - H)    
        p_tmp = p_tmp + 273.15    
        Tm = 0.673 * p_tmp + 83.0
        C = 1E6 * mw / (k2 - k1 * eps0 + k3 / Tm) / RS
        e0 = np.exp(1.81 + 17.27 * DPT / (DPT + 237.5))
        ZWDS = 0.002277 * (0.005 + 1255 / p_tmp) * e0            
        return ZWDS * C * 100.  
    
    @staticmethod
    def getClassfitsColumnsZip():
        " returns name, format, unit for every classfits entry"    
        LIST_TTYPE = \
            ["MJD",
             "MAXIS1", "SCAN", "TELESCOP", "TSYS",
             "IMAGFREQ", "DELTAV", "TAU-ATM", "MH2O",
             "TOUTSIDE", "PRESSURE", "CRVAL2", "CRVAL3",
             "ELEVATIO", "AZIMUTH", "DATE-OBS", "UT",
             "LST", "OBSTIME", "CRPIX1", "RESTFREQ",
             "OBJECT", "VELOCITY", "CDELT1", "CDELT2",
             "CDELT3", "LINE", "SIGNAL", "CAL_IS_ON",
             "CALTEMP", "SPECTRUM"]
        
        LIST_TFORM = \
            ["1D",
             "1J", "1J", "12A", "1E",
             "1E", "1E", "1E", "1E",
             "1E", "1E", "1E", "1E",
             "1E", "1E", "23A ", "1D",
             "1D", "1E", "1E", "1D",
             "12A", "1E", "1E", "1D",
             "1D", "12A", "1J", "1J",
             "1D" , "ND" ]
        
        LIST_TUNIT = \
            ["d",
             " ", "", "", "K",
             "Hz", "m.s-1", "neper", "mm",
             "K", "hPa", "deg", "deg",
             "deg", "deg", "", "s",
             "s", "s", "", "Hz",
             "", "m.s-1", "Hz", "deg",
             "deg", "",  "", "", "K", "K"]
            
        return list(zip(LIST_TTYPE, LIST_TFORM, LIST_TUNIT))      
    
    @staticmethod
    def classfitsColumnsModel(additional_columns=None, **kwargs):
        """classfits structure model"""
        
        model_primary_header = """
        SIMPLE  =                    T
        BITPIX  =                    8
        NAXIS   =                    0
        EXTEND  =                    T
        BLOCKED =                    T
        ORIGIN  = 'SRT'
        CREATOR = '      '
        END
        """

        model_header = """
        XTENSION= 'BINTABLE'
        BITPIX  =                    8         / Always 8.
        NAXIS   =                    2         / Always 2: tables have 2 dimensions.
        NAXIS1  =                 7275         / Number of bytes per row.
        NAXIS2  =                    4         / Number of rows.
        PCOUNT  =                    0         / Usually 0.
        GCOUNT  =                    1         / Always 1.
        TFIELDS =                   18         / Number of columns.
        EXTNAME = 'MATRIX  '                   / Just a name, not important.
        EXTVER  =                    1         / Always 1.
        MAXIS   =                    4         / Number of dimensions in the data.
        MAXIS1  =                 1793         / Dummy number of channels (see TTYPE1).
        MAXIS2  =                    1         /
        MAXIS3  =                    1         /
        MAXIS4  =                    1         /
        CTYPE1  = 'FREQ    '                   / Dim1: freq => MAXIS1 = Nb channels.
        CRVAL1  =  0.0000000000000E+00         / Frequency offset, always 0.
        CDELT1  =  0.0000000000000E+00         / Frequency resolution [Hz].
        CRPIX1  =  0.0000000000000E+00         / Dummy reference channel (see TTYPE18).
        CTYPE2  = 'RA      '
        CRVAL2  =  0.0000000000000E+00
        CDELT2  =  0.0000000000000E+00
        CRPIX2  =  0.0000000000000E+00
        CTYPE3  = 'DEC     '
        CRVAL3  =  0.0000000000000E+00
        CDELT3  =  0.0000000000000E+00
        CRPIX3  =  0.0000000000000E+00
        CTYPE4  = 'STOKES  '
        CRVAL4  =  0.0000000000000E+00
        CDELT4  =  0.0000000000000E+00
        CRPIX4  =  0.0000000000000E+00
        SUBSCAN =                    1         / Subscan number.  Often 1.
        LINE    = '            '               / Name of your line, up to 12 chars.
        OBJECT  = '            '               / Name of your source, up to 12 chars.
        RESTFREQ=  0.0000000000000E+00         / Rest (signal) frequency at ref chan.
        VELO-HEL=  0.0000000000000E+00         / Velocity at ref.  chan [m.s-1].
        VELDEF  = 'RADI-LSR'                   / Type of velocity.
        GAINIMAG=  0.0000000000000E+00         / Ratio Image/Signal.
        BEAMEFF =  0.0000000000000E+00         / Beam efficiency.
        FORWEFF =  0.0000000000000E+00         / Forward efficiency.
        EPOCH   =  2.0000000000000E+03         / Epoch of coordinates.
        DATE-RED= '15/07/97'                   / Date of reduction.
        """
        
        LIST_TTYPE = \
            ["MJD",
             "MAXIS1", "SCAN", "TELESCOP", "TSYS",
             "IMAGFREQ", "DELTAV", "TAU-ATM", "MH2O",
             "TOUTSIDE", "PRESSURE", "CRVAL2", "CRVAL3",
             "ELEVATIO", "AZIMUTH", "DATE-OBS", "UT",
             "LST", "OBSTIME", "CRPIX1", "RESTFREQ",
             "OBJECT", "VELOCITY", "CDELT1", "CDELT2",
             "CDELT3", "LINE", "SIGNAL", "CAL_IS_ON",
             "CALTEMP"]
        
        LIST_TFORM = \
            ["1D",
             "1J", "1J", "12A", "1E",
             "1E", "1E", "1E", "1E",
             "1E", "1E", "1E", "1E",
             "1E", "1E", "23A ", "1D",
             "1D", "1E", "1E", "1D",
             "12A", "1E", "1E", "1D",
             "1D", "12A", "1J", "1J",
             "1D"]
        
        LIST_TUNIT = \
            ["d",
             " ", "", "", "K",
             "Hz", "m.s-1", "neper", "mm",
             "K", "hPa", "deg", "deg",
             "deg", "deg", "", "s",
             "s", "s", "", "Hz",
             "", "m.s-1", "Hz", "deg",
             "deg", "",  "", "", "K"]
                         
        cols = []
        list_ttype = LIST_TTYPE
        list_tform = LIST_TFORM
        list_tunit = LIST_TUNIT
    
        for ttype, tform, tunit in zip(list_ttype, list_tform, list_tunit):
            newcol = fits.Column(name=ttype, format=tform, unit=tunit)
            cols.append(newcol)
        coldefs = fits.ColDefs(cols)
        if additional_columns is not None:
            coldefs += fits.ColDefs(additional_columns)
    
        hdu = fits.BinTableHDU.from_columns(
            coldefs, header=fits.Header.fromstring(model_header, sep='\n'),
            name='MATRIX', **kwargs)
    
        primary_hdu = fits.PrimaryHDU(
            header=fits.Header.fromstring(model_primary_header, sep='\n'))
        return fits.HDUList([primary_hdu, hdu])