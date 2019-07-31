import h5py

def export_h5_dataset(export_fname, label, data):
    data_file = h5py.File(export_fname, 'a')
    data_file.create_dataset(label, data=data)
    data_file.close()

def import_h5_dataset(import_fname, label):
    data_file = h5py.File(import_fname, 'r')
    data = data_file[label][...]
    data_file.close()
    return data