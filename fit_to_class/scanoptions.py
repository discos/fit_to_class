import os
import logging
from fit_to_class import fitslike_scangeometry
from fit_to_class import fitslike_scantype


class ScanOptions:
    """
    Scanning parameters 
    it takes global script variable lile scan_geometry and scan_types    
    """

    # Scan and type structures
    scan_geometry= fitslike_scangeometry.ScanGeometry()
    available_scan_geo_elements= scan_geometry.get_types()
    scan_geo_rules= scan_geometry.get_rules()
    scan_types= fitslike_scantype.ScanType()
    available_scan_types= scan_types.get_types()
    output_fname_suffix= '_classfits'

    def __init__(self, p_logger):
        self.folder="."
        self._output_path=""
        self.scan_conf_path=""
        self.type= None
        self.feed= 0
        self.raw=False
        self.geometry= None
        self.parallel= 4
        self._errors= False
        self._logger= p_logger                


    def __str__(self):
        return "SCAN OPTIONS:\n" + \
                "Folder: " + self.folder +"\n"\
                "Output Path " + self._output_path +"\n"\
                "Scan option conf path" + self.scan_conf_path + "\n"\
                "Scan Type: " + ScanOptions.scan_types.get_str_type(self.type) +"\n"\
                "Geometry: " + str(self.geometry.get_geometry()) +"\n"\
                "Feed: " + str(self.feed) +"\n"\
                "parallelism: " + str(self.parallel) + "\n"

    # OUTPUT PATH

    def get_output_path(self) -> str:
        return self._output_path

    def get_errors(self) -> bool:
        return self._errors

    # FOLDER

    @property
    def folder(self):
        return self._folder

    @folder.setter
    def folder(self, p_folder):
        if type(p_folder) is not str:            
            raise ValueError("Folder argument shoud be of string type")            
        if not os.path.isdir(p_folder):            
            raise ValueError("Folder argument shoud be a valid path")            
        self._folder= p_folder.rstrip('/')              
        self._output_path= os.path.join(self._folder,ScanOptions.output_fname_suffix +'/')          


    # SCAN PIPELINE CONF PATH

    @property
    def scan_conf_path(self):
        return self._scan_conf_path

    @scan_conf_path.setter
    def scan_conf_path(self, p_path):
        if type(p_path) is not str:            
            raise ValueError("Scan conf file path argument shoud be of string type")       
        if not p_path:
            self._scan_conf_path= p_path
            return
        if not os.path.isfile(p_path):            
            raise ValueError("Scan conf file path argument shoud be a valid path")            
        self._scan_conf_path= p_path        

    # SCAN TYPE

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, p_type_str):
        self._type= ScanOptions.scan_types.get_enum_type(p_type_str)
    
    # SCAN GEOMETRY

    @property
    def geometry(self):
        return self._geometry
    
    @geometry.setter
    def geometry(self, p_geometry_str):
        if not ScanOptions.scan_geometry.parse(p_geometry_str):
            raise ValueError("Invalid geometry format or values")            
        self._geometry= ScanOptions.scan_geometry

    # RAW SCAN

    @property
    def raw(self):
        return self._raw

    @raw.setter
    def raw(self, p_value):
        if type(p_value) is not bool:
            raise ValueError("Raw scan option nees BOOL type")
        self._raw= p_value

    # SPECIFIC FEED 

    @property
    def feed(self):
        return self._feed

    @feed.setter
    def feed(self, p_feed):
        if p_feed is None:
            self._feed= None
            return
        if type(p_feed) is not int:
            raise ValueError("Feed specification need INT type")
        self._feed= p_feed

    
    # PARALLELISM

    @property
    def parallel(self):
        return self._parallel
    
    @parallel.setter
    def parallel(self, p_par):
        if type(p_par) is not int:
            raise ValueError("Parallel computation requires INT value")
        self._parallel= p_par
    
