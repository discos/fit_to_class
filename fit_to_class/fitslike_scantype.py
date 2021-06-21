from enum import Enum

class ScanTypes(Enum):
    NOTVALID=0
    ONOFF= 1
    NODDING= 2
    OTF= 3
    RASTER = 4

class ScanType:
    """
    Scan types
    """    

    def __init__(self):
        self.__types= {}
        self.__types["ONOFF"] = ScanTypes.ONOFF
        self.__types["NODDING"] = ScanTypes.NODDING
        self.__types["OTF"] = ScanTypes.OTF
        self.__types["RASTER"]= ScanTypes.RASTER
    
    def get_types(self) -> list:
        """ Getter available string values"""
        valid_types= []
        for key,value in self.__types.items():
            valid_types.append(key)
        return valid_types

    def get_enum_type(self, p_type_str):
        try:
            return self.__types[p_type_str]
        except:
            return ScanTypes.NOTVALID

    def get_str_type(self, p_type) -> str:
        for key, value in self.__types.items():
            if p_type == key:
                return key
        return "Unknown"

    