import astropy.units as unit

def _tableToDict(p_table) -> dict:
            """
            Disgregate :) and build a single row dictionary
            Gets unit from column names _u_unit
            """
            "  "
            _out= {}
            if len(p_table) == 0:
                print("Error: Empty table..why?")
                return _out
            for field in p_table.colnames:
                splits= field.split('_u_')
                if len(splits) > 1:
                    _out[splits[0]]= p_table[field].data.data
                    _out[splits[0]] *= unit.Unit(splits[1])
                else:
                    try:
                        _out[field]= p_table[field].data.data
                    except TypeError:
                        # Cannot average strings
                        _out[field]= p_table[field][0]
            return _out
