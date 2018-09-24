import numpy as np
from opencmiss.iron import iron
import fields

def interpolate_field(field, element_ids=[], num_values=4,dimension=3, derivative_number=1, xi=None, elems=None):

    if xi is None:
        XiNd = fields.generate_xi_grid_fem(num_points=num_values)

        num_elem_values = XiNd.shape[0]
        num_Xe = len(element_ids)
        total_num_values = num_Xe * num_elem_values
        values = np.zeros((num_Xe, num_elem_values, dimension))
        xi = np.zeros((num_Xe, num_elem_values, dimension))
        elements = np.zeros((num_Xe, num_elem_values, 1))

        for elem_idx, element_id in enumerate(element_ids):
            for point_idx in range(num_elem_values):
                single_xi = XiNd[point_idx,:]
                values[elem_idx, point_idx, :] = field.ParameterSetInterpolateSingleXiDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, derivative_number, element_id, single_xi, dimension)
            xi[elem_idx, :, :] = XiNd
            elements[elem_idx, :] = element_id

        values = np.reshape(values, (total_num_values, dimension))
        xi = np.reshape(xi, (total_num_values, dimension))
        elements = np.reshape(elements, (total_num_values))
        return values, xi, elements
    else:
        num_values = xi.shape[0]
        values = np.zeros((num_values, dimension))
        for point_idx in range(xi.shape[0]):
            element_id = elems[point_idx]
            single_xi = xi[point_idx,:]
            values[point_idx, :] = field.ParameterSetInterpolateSingleXiDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, derivative_number, int(element_id), single_xi, dimension)
        return values