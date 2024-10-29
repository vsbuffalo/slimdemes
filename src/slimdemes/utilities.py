def safe_float(value):
    """
    Return None if value is None/"None" otherwise
    convert to float.
    """
    if isinstance(value, str) and value.lower() == "none":
        return None
    if value is None:
        return None
    return float(value)
