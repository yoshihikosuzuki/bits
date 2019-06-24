import pickle


def load_pickle(pkl_fname):
    with open(pkl_fname, 'rb') as f:
        return pickle.load(f)


def save_pickle(obj, pkl_fname):
    with open(pkl_fname, 'wb') as f:
        pickle.dump(obj, f, protocol=4)

