

def set_value(p_dict, p_key, p_value, p_raise= True ):
    """
    Dict key wrapper
    """
    if p_key not in p_dict:
        if p_raise:
            raise Exception(f"Missing key {p_key}")
    else:
        p_dict[p_key]= p_value

def get_value(p_dict, p_key, p_raise= True):
    """
    Get value from dict
    """
    if p_key not in p_dict:
        if p_raise:            
            raise Exception(f"Missing key {p_key}")
        else:
            return "KEY_UNAVAILABLE"
    # Raise with None value
    if p_dict[p_key]== None:
        raise Exception(f"Key with unacceptable None value: {p_key}")
    return p_dict[p_key]