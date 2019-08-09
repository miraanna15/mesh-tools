import numpy as np
import mesh_tools
import morphic

def exfile_to_morphic(nodeFilename, elementFilename, coordinateField,
                      dimension=2, interpolation='linear',elements=[]):
    """Convert an exnode and exelem files to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    nodeFilename -- exnode filename
    elementFilename -- exelem filename
    coordinateField -- the field to read in
    dimension -- dimension of mesh to read in
    """

    # Create morphic mesh
    mesh = morphic.Mesh()

    # Load exfiles
    exnode = mesh_tools.Exnode(nodeFilename)
    exelem = mesh_tools.Exelem(elementFilename, dimension)

    if elements == []:
        ex_elems = exelem.elements
    else:
        ex_elems = []
        for elem in exelem.elements:
            if elem.number in elements:
                ex_elems.append(elem)

    node_ids = []
    for elem_idx, elem in enumerate(ex_elems):
        node_ids.append(elem.nodes)
    node_ids = np.unique(node_ids)

    # Add nodes
    if interpolation == 'hermite':
        derivatives = range(1,9)
    else:
        derivatives = [1]

    for node_num in node_ids:
        coordinates = []
        for component in range(1, 4):
            component_name = ["x", "y", "z"][component - 1]
            componentValues = []
            for derivative_idx, derivative in enumerate(derivatives):
                componentValues.append(exnode.node_value(coordinateField,
                                                     component_name, node_num,
                                                     derivative))
            coordinates.append(componentValues)
        if interpolation != 'hermite':
            coordinates = np.hstack(coordinates)

        mesh.add_stdnode(node_num, coordinates, group='_default')
        #print('Morphic node added', node_num, coordinates)

    if dimension == 2:
        if interpolation == 'linear':
            element_interpolation = ['L1', 'L1']
        if interpolation == 'quadratic':
            element_interpolation = ['L2', 'L2']
    elif dimension == 3:
        if interpolation == 'linear':
            element_interpolation = ['L1', 'L1', 'L1']
        if interpolation == 'quadratic':
            element_interpolation = ['L2', 'L2', 'L2']
        if interpolation == 'cubic':
            element_interpolation = ['L3', 'L3', 'L3']
        if interpolation == 'hermite':
            element_interpolation = ['H3', 'H3', 'H3']

    # Add elements
    for elem in ex_elems:
        mesh.add_element(elem.number, element_interpolation, elem.nodes)
        #print('Morphic element added', elem.number)

    # Generate the mesh
    mesh.generate(True)

    return mesh

def exfile_to_OpenCMISS(nodeFilename, elementFilename, coordinateField, basis,
                        region, meshUserNumber, dimension=2,
                        interpolation='linear', pressure_basis=None,
                        use_pressure_basis=False, elements=[]):
    """Convert an exnode and exelem files to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    nodeFilename -- exnode filename
    elementFilename -- exelem filename
    coordinateField -- the field to read in
    dimension -- dimension of mesh to read in
    interpolation -- the interpolation of the mesh to read in
    """
    from opencmiss.iron import iron

    # Load exfiles
    exnode = mesh_tools.Exnode(nodeFilename)
    exelem = mesh_tools.Exelem(elementFilename, dimension)

    if elements == []:
        ex_elems = exelem.elements
    else:
        ex_elems = []
        for elem in exelem.elements:
            if elem.number in elements:
                ex_elems.append(elem)

    totalNumberOfNodes = len(exnode.nodeids)
    totalNumberOfElements = len(ex_elems)

    # Start the creation of a manually generated mesh in the region
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, dimension)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(totalNumberOfElements)

    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, totalNumberOfNodes)
    nodes.UserNumbersAllSet(exnode.nodeids)
    nodes.CreateFinish()

    MESH_COMPONENT1 = 1
    MESH_COMPONENT2 = 2

    if (use_pressure_basis):
        mesh.NumberOfComponentsSet(2)
    else:
        mesh.NumberOfComponentsSet(1)

    elements = iron.MeshElements()
    elements.CreateStart(mesh, MESH_COMPONENT1, basis)

    elemNums = []
    # elemNums = range(1,len(ex_elems))
    for elem in ex_elems:
        elemNums.append(elem.number)

    elements.UserNumbersAllSet(elemNums)
    node_ids = []

    for elem_idx, elem in enumerate(ex_elems):
        elements.NodesSet(elem_idx+1, elem.nodes)
        node_ids.append(elem.nodes)
    elements.CreateFinish()
    node_ids = np.unique(node_ids)

    if (use_pressure_basis):
        linear_elem_node_idxs = [0, 3, 12, 15, 48, 51, 60, 63]
        pressureElements = iron.MeshElements()
        pressureElements.CreateStart(mesh, MESH_COMPONENT2, pressure_basis)
        pressureElements.UserNumbersAllSet(elemNums)
        for elem_idx, elem in enumerate(ex_elems):
            pressureElements.NodesSet(elem_idx+1, np.array(
                    elem.nodes, dtype=np.int32)[linear_elem_node_idxs])
        pressureElements.CreateFinish()

    mesh.CreateFinish()

    coordinates, node_ids = mesh_tools.extract_exfile_coordinates(nodeFilename, coordinateField,node_ids, interpolation)

    return mesh, coordinates, node_ids


def morphic_to_OpenCMISS(morphicMesh, region, basis, meshUserNumber,
                         dimension=2, interpolation='linear',
                         UsePressureBasis=False, pressureBasis=None,
                         include_derivatives=True, elements=[]):
    """Convert an exnode and exelem files to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    morphicMesh -- morphic mesh
    dimension -- dimension of mesh to read in
    include_derivatives -- whether to include derivatives when returning mesh
                           coordinates.
    """
    from opencmiss.iron import iron
    if elements == []:
        element_list = [int(cid)+1  for cid in morphicMesh.get_element_cids()]
    else:
        element_list = []
        for element in morphicMesh.elements:
            if element.id in elements:
                element_list.append(element.id+1)
    # Create mesh topology
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, dimension)
    if (UsePressureBasis):
        mesh.NumberOfComponentsSet(2)
    else:
        mesh.NumberOfComponentsSet(1)

    node_list = []
    for element_idx in element_list:
        node_list.append( morphicMesh.elements[element_idx-1].node_ids)
    node_list = np.unique(node_list)+1


    mesh.NumberOfElementsSet(len(element_list))
    nodes = iron.Nodes()
    nodes.CreateStart(region, len(node_list))
    nodes.UserNumbersAllSet((np.array(node_list)).astype('int32'))
    nodes.CreateFinish()

    MESH_COMPONENT1 = 1
    MESH_COMPONENT2 = 2
    elements = iron.MeshElements()
    elements.CreateStart(mesh, MESH_COMPONENT1, basis)
    elements.UserNumbersAllSet((np.array(element_list).astype('int32')))

    for element_idx in element_list:
        element = morphicMesh.elements[element_idx-1]
        elements.NodesSet(element_idx, np.array(element.node_ids, dtype='int32')+1)
    elements.CreateFinish()

    if (UsePressureBasis):
        linear_elem_node_idxs = [0, 3, 12, 15, 48, 51, 60, 63]
        pressureElements = iron.MeshElements()
        pressureElements.CreateStart(mesh, MESH_COMPONENT2, pressureBasis)
        pressureElements.UserNumbersAllSet(element_list)
        for element_idx in element_list:
            pressureElements.NodesSet(element_idx,
                             np.array(morphicMesh.elements[element_idx-1].node_ids, dtype='int32')[linear_elem_node_idxs]+1)
        pressureElements.CreateFinish()

    mesh.CreateFinish()

    # Add nodes
    if interpolation == 'hermite':
        derivatives = range(1,9)
    else:
        derivatives = [1]

    if include_derivatives:
        coordinates = np.zeros((len(node_list), 3,len(derivatives)))
    else:
        derivatives = [1]
        coordinates = np.zeros((len(node_list), 3))

    morphicNodes = morphicMesh.nodes
    _, morphicNodesIds = morphicMesh.get_node_ids()

    for node_indx,node_id in enumerate(node_list):

        for component_idx, component in enumerate(range(1, 4)):
            for derivative_idx, derivative in enumerate(derivatives):

                if include_derivatives:
                    coordinates[node_indx,component_idx, derivative_idx] = \
                        morphicNodes[node_id-1].values[component_idx]
                else:
                    coordinates[node_indx,component_idx] = \
                        morphicNodes[node_id-1].values[component_idx]


    return mesh, coordinates, node_list


def OpenCMISS_to_morphic(c_mesh, geometric_field,
                      dimension=2, interpolation='linear'):
    """Convert an OpenCMISS mesh to a morphic mesh.

    Only Linear lagrange elements supported.

    Keyword arguments:
    morphicMesh -- morphic mesh
    dimension -- dimension of mesh to read in
    """

    from opencmiss.iron import iron
    mesh_nodes = iron.MeshNodes()
    mesh_elements = iron.MeshElements()
    c_mesh.NodesGet(1, mesh_nodes)
    c_mesh.ElementsGet(1, mesh_elements)
    # Create morphic mesh
    mesh = morphic.Mesh()

    mesh_elements.NodesGet()

    # Load exfiles

    # Add nodes
    if interpolation == 'linear':
        derivatives = [1]
    elif interpolation == 'hermite':
        derivatives = range(1,9)
    for node_num in exnode.nodeids:
        coordinates = []
        for component in range(1, 4):
            component_name = ["x", "y", "z"][component - 1]
            componentValues = []
            for derivative_idx, derivative in enumerate(derivatives):
                componentValues.append(
                    geometric_field.ParameterSetGetNodeDP(
                            iron.FieldVariableTypes.U,
                            iron.FieldParameterSetTypes.VALUES, 1, derivative,
                            node_num, component))
            coordinates.append(componentValues)

        mesh.add_stdnode(node_num, coordinates, group='_default')
        #print('Morphic node added', node_num, coordinates)

    if dimension == 2:
        if interpolation == 'linear':
            element_interpolation = ['L1', 'L1']
        if interpolation == 'quadratic':
            element_interpolation = ['L2', 'L2']
    elif dimension == 3:
        if interpolation == 'linear':
            element_interpolation = ['L1', 'L1', 'L1']
        if interpolation == 'quadratic':
            element_interpolation = ['L2', 'L2', 'L2']
        if interpolation == 'hermite':
            element_interpolation = ['H3', 'H3', 'H3']

    # Add elements
    for elem in exelem.elements:
        mesh.add_element(elem.number, element_interpolation, elem.nodes)
        #print('Morphic element added', elem.number)

    # Generate the mesh
    mesh.generate(True)

    return mesh