import pickle

def read_pickle_file(file):
    pickle_data = pickle.load(file)
    return pickle_data

def write_pickle_file(a, file):
	with open(file, 'wb') as handle:
    	pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)