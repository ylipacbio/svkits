VERSION = (0, 1, 0)

def get_version():
    return ".".join([str(i) for i in VERSION])

__version__ = get_version()
