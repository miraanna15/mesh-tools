import os
import mesh_tools
from morphic.utils import convert_hermite_lagrange
from opencmiss.iron import iron

# Parameters
dimension = 3
results_folder='./results/'
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# Read in ex mesh as a morphic mesh
coordinateField = 'coordinates'
nodeFilename = 'test_mesh/lung_mesh.exnode'
elementFilename = 'test_mesh/lung_mesh.exelem'
cubic_hermite_morphic_mesh = mesh_tools.exfile_to_morphic(nodeFilename, elementFilename, coordinateField,
                      dimension=dimension, interpolation='hermite')

#Convert mesh from cubic hermite to cubic Lagrange
cubic_lagrange_morphic_mesh = convert_hermite_lagrange(cubic_hermite_morphic_mesh, tol=1e-9)

# Export mesh in ex format using OpenCMISS
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

coordinateSystemUserNumber = 1
regionUserNumber = 3
basisUserNumber = 2
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1

coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(dimension)
basis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*dimension)
basis.QuadratureNumberOfGaussXiSet([4]*dimension)
basis.CreateFinish()

region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

mesh, coordinates, node_nums, element_nums = mesh_tools.morphic_to_OpenCMISS(cubic_lagrange_morphic_mesh, region, basis, meshUserNumber,
                      dimension=dimension, interpolation='cubicLagrange',UsePressureBasis=False, pressureBasis=None)

decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CreateFinish()

geometric_field = iron.Field()
geometric_field.CreateStart(geometricFieldUserNumber, region)
geometric_field.MeshDecompositionSet(decomposition)
geometric_field.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometric_field.CreateFinish()

# Update the geometric field parameters
geometric_field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES)
for node_idx, node in enumerate(node_nums):
    for component_idx, component in enumerate([1, 2, 3]):
        for derivative_idx, derivative in enumerate(
                range(1, coordinates.shape[2] + 1)):
            geometric_field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                    iron.FieldParameterSetTypes.VALUES,
                                                    1, derivative, node,
                                                    component,
                                                     coordinates[
                                                        node_idx, component_idx, derivative_idx])
geometric_field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES)

output_file = os.path.join(results_folder, 'cubic_lagrange_mesh')
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport(output_file, "FORTRAN")
fields.ElementsExport(output_file, "FORTRAN")
fields.Finalise()