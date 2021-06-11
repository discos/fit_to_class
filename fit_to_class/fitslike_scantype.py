class ScanType:
    """
    Scan types
    """
    NOTVALID=0
    ONOFF= 1
    NODDING= 2
    OTF= 3
    RASTER = 4

    def __init__(self):
        self.__types= {}
        self.__types["ONOFF"] = ScanType.ONOFF
        self.__types["NODDING"] = ScanType.NODDING
        self.__types["OTF"] = ScanType.OTF
        self.__types["RASTER"]= ScanType.RASTER
    
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
            return ScanType.NOTVALID

    