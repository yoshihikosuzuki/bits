import configparser
import pickle


def load_config(config_fname):
    config = configparser.ConfigParser()
    config.read(config_fname)
    return config


def load_pickle(pkl_fname):
    with open(pkl_fname, 'rb') as f:
        return pickle.load(f)


def save_pickle(obj, pkl_fname):
    with open(pkl_fname, 'wb') as f:
        pickle.dump(obj, f, protocol=4)

