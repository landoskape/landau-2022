import os


def get_src_dir():
    return os.path.abspath(os.path.dirname(__file__))

def get_save_dir():
    return os.path.join(os.path.dirname(get_src_dir()), "saved_data")
