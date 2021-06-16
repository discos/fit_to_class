
import re


def get_geo_loop_len(p_geo_list) ->int:
    """
    Len calc from geo tuples list
    """
    sum= 0        
    for tup in p_geo_list:
        try:
            sum= sum + tup[0]
        except:
            pass
    return sum

class ScanGeometry:
    """ 
    Scan Geometry property definition
    """
    ON=0
    OFF=1
    CALON=2
    CALOFF=3

    def __init__(self) -> None:
        self.__geometry= []
        self.__types={}
        self.__types[ScanGeometry.ON]= "ON"
        self.__types[ScanGeometry.OFF]= "OFF"
        self.__types[ScanGeometry.CALON]= "CALON"
        self.__types[ScanGeometry.CALOFF]= "CALOFF"
        self.__elementPattern= "(\d+)([A-Z]+)"

    def get_geometry(self) -> list:
        return self.__geometry.copy()

    def get_types(self) -> list:
        """
        Getter vaild string types
        """
        valid_types= []
        for key,value in self.__types.items():
            valid_types.append(value)
        return valid_types

    def get_rules(self) -> str:
        """
        Getter geo rules
        """
        return " Number+Element comma separated ex: 1ON,1OFF,1CALON. Empty is a valid geometry"

    def parse(self, p_geo_str) -> bool:
        """
        Input geo string parsing
        Format:
            ON
            OFF
            CAL_ON
            CAL_OFF

            ex: 2ON,2OFF,1CAL_ON
            ex: 1ON

        Criteria:
            Empty string is valid
            Not empty not parsable string in invalid
            On non empty string we need at least 1 type (ex: 1ON)

        Return:
            True: input format valid            
            
        """        
        if p_geo_str is None:
            return True
        if type(p_geo_str) is not str:
            return False
        if not p_geo_str:
            return True
        splitted= p_geo_str.split(",")
        filtered= list(filter(None,splitted))
        if not filtered:
            return False
        self.__geometry.clear()
        for el in filtered:
            self.__geometry.append(self._validateGeoElement(self.__elementPattern, el))
        for tup in self.__geometry:
            if 0 in tup:
                return False
            if None in tup:
                return False
        return True

    def get_loop_lenght(self) -> int:
        """
        Getter geo loop file len
        It's valid only for non otf scans
        """
        sum= 0        
        for tup in self.__geometry:
            try:
                sum= sum + tup[0]
            except:
                pass
        return sum

    def _validateGeoElement(self, p_pattern, p_element) -> tuple:
            """
            Single geo element validation
            Assuming is string
            NUMBER+TYPE
            ex: 2OFF 
            """
            match= re.search(p_pattern, p_element)
            if not match:
                return (0,None)            
            return (match.group(1), match.group(2))
