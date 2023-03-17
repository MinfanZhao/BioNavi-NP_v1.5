import pickle
def prepare_rxn_fps(filename):
    with open(filename, 'rb') as f:
        rxn_fps = pickle.load(f)
    return rxn_fps
    #return [pair[1] for pair in rxn_fps]