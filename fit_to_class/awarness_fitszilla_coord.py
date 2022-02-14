import numpy as np
import astropy.units as unit
from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import copy

def _feeds_rest_angle(p_xoffsets, p_yoffsets):
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

def _observing_angle(p_rest_angle, p_derot_angle):
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

def _offset_needs_correction(p_coord):
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

def _offset_corrections(p_obs_angle, p_xoffset, p_yoffset):
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

def _azel_to_radec(p_obs_data,
                    p_obstimes,
                    p_el, p_az,
                    p_xoffs, p_yoffs,
                    p_location):
    """
    Conversion from alt-az to ra-dec.
    Offset must be correccted based on observation time.
    Returns:
    --------
    Actual ra dec lists            
    
    """            
    l_el = copy.deepcopy(p_el)
    l_az = copy.deepcopy(p_az)            
    
    # Basic version, header fitszilla version < V1.21                                     
    # TODO è corretto fare la conversione standard con gal? non si tengono in considerazione gli offset..
    if p_obs_data['user_offset_frame']== 'GAL':
        print("ERROR: Cannot convert coordinates! skupping..data will result as incomplete?")
        # feed offset (corrected)relates to az, el they need recalc to add az
        l_el += p_yoffs
        l_az += p_xoffs / np.cos(l_el)
        # TODO Conversion to adapt orign to az, el allowing data to fit ICRS !?
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
                
    # Pre Apply offset for HOR
    if p_obs_data['user_offset_frame']== 'HOR':        
        l_el = l_el + p_yoffs - p_obs_data['user_lat_offset']
        l_az = l_az +( p_xoffs - p_obs_data['user_lon_offset'] ) / np.cos(l_el)                
    
    # Frame change
    l_coordsAltAz = AltAz(az=Angle(l_az), alt=Angle(l_el), location= p_location, obstime= p_obstimes)
    l_coords_deg = l_coordsAltAz.transform_to(ICRS)
    
    # Post apply offset for
    if p_obs_data['user_offset_frame']== 'EQ':        
        l_dec = np.radians(l_coords_deg.dec) - p_obs_data['user_lat_offset']
        l_ra = np.radians(l_coords_deg.ra) - p_obs_data['user_lon_offset'] / np.cos(l_dec)
    else:
        l_ra = np.radians(l_coords_deg.ra)
        l_dec = np.radians(l_coords_deg.dec)
                
    return l_ra, l_dec

