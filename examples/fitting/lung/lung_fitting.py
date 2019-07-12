import h5py
import numpy as np

import mesh_tools

def load_data():
    filename = '/home/psam012/opt/cvt-data/2018-10-30 pre CT lung/PointCloudAlignedToLaserLinesPreCT.hdf5'
    coords_file = h5py.File(filename, 'r')
    geometric_coords = np.array((coords_file[coords_file.keys()[0]]))
    coords_file.close()
    return geometric_coords

nodeFilename = '/home/psam012/opt/cvt-data/2018-11-9 Post microCT/MeshFitting/data/left_surface.exnode'
elementFilename = '/home/psam012/opt/cvt-data/2018-11-9 Post microCT/MeshFitting/data/left_surface.exelem'

mesh_tools.exfile_to_OpenCMISS(nodeFilename, elementFilename, coordinateField, basis, region, meshUserNumber,
                      dimension=2, interpolation='linear'):