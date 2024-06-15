#!/usr/bin/env python
#
# This is an example script for a linear elasticity cantilever bea  using OpenCMISS calls in python.
# By Chris Bradley
#
#

import sys

# Intialise OpenCMISS
from opencmiss.iron import iron

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------

HEIGHT = 2.0 # mm
WIDTH = 2.0 # mm
LENGTH = 10.0 # mm

YOUNGS_MODULUS = 30.0E6 # mg.mm^-1.ms^-2
POISSONS_RATIO = 0.3
THICKNESS = 1.0 # mm (for plane strain and stress)

LINEAR_LAGRANGE = 1
QUADRATIC_LAGRANGE = 2
CUBIC_LAGRANGE = 3
CUBIC_HERMITE = 4
LINEAR_SIMPLEX = 5
QUADRATIC_SIMPLEX = 6
CUBIC_SIMPLEX = 7

DIRICHLET_BCS = 1
NEUMANN_BCS = 2

# Boundary condition for 1 & 2D. Analytic for 3D
boundaryConditionType = DIRICHLET_BCS
MAX_DISPLACEMENT = -0.10*HEIGHT;
MAX_FORCE = -10.0 # N.mm^-2

if (boundaryConditionType == DIRICHLET_BCS):
    DISPLACEMENT_BC = MAX_DISPLACEMENT
elif (BOUNDARY_CONDITION_TYPE == NEUMANN_BCS):
    DISPLACEMENT_BC = MAX_FORCE
else:
    print('Invalid boundary condition type')
    exit()
   
(CONTEXT_USER_NUMBER,
 COORDINATE_SYSTEM_USER_NUMBER,
 REGION_USER_NUMBER,
 BASIS_USER_NUMBER,
 GENERATED_MESH_USER_NUMBER,
 MESH_USER_NUMBER,
 DECOMPOSITION_USER_NUMBER,
 DECOMPOSER_USER_NUMBER,
 GEOMETRIC_FIELD_USER_NUMBER,
 ELASTICITY_DEPENDENT_FIELD_USER_NUMBER,
 ELASTICITY_MATERIALS_FIELD_USER_NUMBER,
 ELASTICITY_ANALYTIC_FIELD_USER_NUMBER,
 ELASTICITY_DERIVED_FIELD_USER_NUMBER,
 ELASTICITY_EQUATIONS_SET_FIELD_USER_NUMBER,
 ELASTICITY_EQUATIONS_SET_USER_NUMBER,
 ELASTICITY_PROBLEM_USER_NUMBER) = range(1,17)

NUMBER_OF_GAUSS_XI = 4

numberOfGlobalXElements = 2
numberOfGlobalYElements = 1
numberOfGlobalZElements = 1
interpolationType = LINEAR_LAGRANGE
#interpolationType = LINEAR_SIMPLEX

# Override with command line arguments if need be
if len(sys.argv) > 1:
    if len(sys.argv) > 5:
        sys.exit('ERROR: too many arguments- currently only accepting up to 4 options: numberOfGlobalXElements numberOfGlobalYElements numberOfGlobalZElements interpolationType')
    numberOfGlobalXElements = int(sys.argv[1])
    if len(sys.argv) > 2:
        numberOfGlobalYElements = int(sys.argv[2])
    if len(sys.argv) > 3:
        numberOfGlobalZElements = int(sys.argv[3])
    if len(sys.argv) > 4:
        interpolationType = int(sys.argv[4])

if (numberOfGlobalZElements >= 0):
    if (numberOfGlobalYElements >= 0):
        if (numberOfGlobalXElements >= 0):
            if (numberOfGlobalZElements == 0):
                if(numberOfGlobalYElements == 0):
                    numberOfDimensions = 1
                else:
                    numberOfDimensions = 2
            else:
                numberOfDimensions = 3
        else:
            sys.exit('ERROR: number of global X elements must be greater than 0.')
    else:
        sys.exit('ERROR: number of global Y elements must be greater than 0.')
else:
    sys.exit('ERROR: number of global Z elements must be greater than 0.')

if (interpolationType == LINEAR_LAGRANGE):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE
    numberOfNodesXi = 2
    numberOfGaussXi = 2
elif (interpolationType == QUADRATIC_LAGRANGE):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif (interpolationType == CUBIC_LAGRANGE):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE
    numberOfNodesXi = 4
    numberOfGaussXi = 4
elif (interpolationType == CUBIC_HERMITE):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.CUBIC_HERMITE
    numberOfNodesXi = 2
    numberOfGaussXi = 4
elif (interpolationType == LINEAR_SIMPLEX):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX
    numberOfNodesXi = 2
    gaussOrder = 4
elif (interpolationType == QUADRATIC_SIMPLEX):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX
    numberOfNodesXi = 3
    gaussOrder = 4
elif (interpolationType == CUBIC_SIMPLEX):
    interpolationTypeXi = iron.BasisInterpolationSpecifications.CUBIC_SIMPLEX
    numberOfNodesXi = 4
    gaussOrder = 5
else:
    sys.exit('The interpolation type of ',interpolationType,' is invalid.')

haveHermite = (interpolationType == CUBIC_HERMITE)
haveSimplex = (interpolationType == LINEAR_SIMPLEX or interpolationType == QUADRATIC_SIMPLEX or interpolationType == CUBIC_SIMPLEX)

elementFactor = 1
if (numberOfDimensions == 2):
    if (haveSimplex):
        elementFactor = 2
    numberOfElements = numberOfGlobalXElements*numberOfGlobalYElements*elementFactor
    numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
    numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
    numberOfNodes = numberOfXNodes*numberOfYNodes            
else:
    if (haveSimplex):
        elementFactor = 6
    numberOfElements = numberOfGlobalXElements*numberOfGlobalYElements*numberOfGlobalZElements*elementFactor
    numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
    numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
    numberOfZNodes = numberOfGlobalZElements*(numberOfNodesXi-1)+1
    numberOfNodes = numberOfXNodes*numberOfYNodes*numberOfZNodes

numberOfXi = numberOfDimensions
if (not haveSimplex):
    numberOfGauss = pow(numberOfGaussXi,numberOfXi)

#-----------------------------------------------------------------------------------------------------------
# CONTEXT AND WORLD REGION
#-----------------------------------------------------------------------------------------------------------

context = iron.Context()
context.Create(CONTEXT_USER_NUMBER)

worldRegion = iron.Region()
context.WorldRegionGet(worldRegion)

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

iron.OutputSetOn("LinearCantilever")

iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"",["BoundaryConditionsVariable_NeumannIntegrate"])

# Get the computational nodes information
computationEnvironment = iron.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

worldWorkGroup = iron.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)

#-----------------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(numberOfDimensions)
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# REGION
#-----------------------------------------------------------------------------------------------------------

region = iron.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("Cantilever")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BASIS
#-----------------------------------------------------------------------------------------------------------

basis = iron.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
if (haveSimplex):
    basis.TypeSet(iron.BasisTypes.SIMPLEX)
else:
    basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(numberOfXi)
basis.InterpolationXiSet([interpolationTypeXi]*numberOfXi)
if (haveSimplex):
    basis.QuadratureOrderSet(gaussOrder)
else:
    basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
basis.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# MESH
#-----------------------------------------------------------------------------------------------------------

generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(GENERATED_MESH_USER_NUMBER,region)
generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
if (numberOfDimensions == 2):
    generatedMesh.ExtentSet([LENGTH,HEIGHT])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements])
else:
    generatedMesh.ExtentSet([LENGTH,WIDTH,HEIGHT])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements])
mesh = iron.Mesh()
generatedMesh.CreateFinish(MESH_USER_NUMBER,mesh)

#-----------------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

decomposition = iron.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = iron.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

geometricField = iron.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
if (numberOfDimensions == 3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
# EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create linear elasiticity equations set
elasticityEquationsSetField = iron.Field()
elasticityEquationsSet = iron.EquationsSet()
if (numberOfDimensions == 2):
    elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                           iron.EquationsSetTypes.LINEAR_ELASTICITY,
                                           iron.EquationsSetSubtypes.TWO_DIMENSIONAL_PLANE_STRESS]
else:
    elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                           iron.EquationsSetTypes.LINEAR_ELASTICITY,
                                           iron.EquationsSetSubtypes.THREE_DIMENSIONAL_ISOTROPIC]
elasticityEquationsSet.CreateStart(ELASTICITY_EQUATIONS_SET_USER_NUMBER,region,geometricField,
                         elasticityEquationsSetSpecification,
                         ELASTICITY_EQUATIONS_SET_FIELD_USER_NUMBER,elasticityEquationsSetField)
elasticityEquationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET DEPENDENT
#-----------------------------------------------------------------------------------------------------------

elasticityDependentField = iron.Field()
elasticityEquationsSet.DependentCreateStart(ELASTICITY_DEPENDENT_FIELD_USER_NUMBER,elasticityDependentField)
elasticityDependentField.LabelSet("ElasticityDependent")
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Displacement")
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.T,"Traction")
elasticityEquationsSet.DependentCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET MATERIALS
#-----------------------------------------------------------------------------------------------------------

elasticityMaterialsField = iron.Field()
elasticityEquationsSet.MaterialsCreateStart(ELASTICITY_MATERIALS_FIELD_USER_NUMBER,elasticityMaterialsField)
elasticityMaterialsField.LabelSet("ElasticityMaterials")
elasticityMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Materials")
elasticityEquationsSet.MaterialsCreateFinish()    
# Initialise the analytic field values
elasticityMaterialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                   1,YOUNGS_MODULUS)
elasticityMaterialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                   2,POISSONS_RATIO)
if(numberOfDimensions==2):
    elasticityMaterialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                       3,THICKNESS)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET ANALYTIC
#-----------------------------------------------------------------------------------------------------------

elasticityAnalyticField = iron.Field()
if(numberOfDimensions==3):
    elasticityEquationsSet.AnalyticCreateStart(iron.EquationsSetLinearElasticityAnalyticFunctionTypes.CANTILEVER_END_LOAD,
                                               ELASTICITY_ANALYTIC_FIELD_USER_NUMBER,elasticityAnalyticField)
    elasticityAnalyticField.LabelSet("ElasticityAnalytic")
    elasticityAnalyticField.VariableLabelSet(iron.FieldVariableTypes.U,"Analytic")
    elasticityEquationsSet.AnalyticCreateFinish()    
    # Initialise the analytic field values
    elasticityAnalyticField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      1,LENGTH)
    elasticityAnalyticField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      2,HEIGHT)
    elasticityAnalyticField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      3,WIDTH)
    elasticityAnalyticField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      4,YOUNGS_MODULUS)
    elasticityAnalyticField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      5,MAX_FORCE)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET DERIVED
#-----------------------------------------------------------------------------------------------------------

# Create a field for the derived field. Have three variables U - Small strain tensor, V - Cauchy stress, W - Elastic Work
if(numberOfDimensions==2):
    numberOfTensorComponents = 3
else:
    numberOfTensorComponents = 6
elasticityDerivedField = iron.Field()
elasticityDerivedField.CreateStart(ELASTICITY_DERIVED_FIELD_USER_NUMBER,region)
elasticityDerivedField.LabelSet("ElasticityDerived")
elasticityDerivedField.TypeSet(iron.FieldTypes.GENERAL)
elasticityDerivedField.DecompositionSet(decomposition)
elasticityDerivedField.GeometricFieldSet(geometricField)
elasticityDerivedField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
elasticityDerivedField.NumberOfVariablesSet(3)
elasticityDerivedField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.V,iron.FieldVariableTypes.W])
elasticityDerivedField.VariableLabelSet(iron.FieldVariableTypes.U,"SmallStrain")
elasticityDerivedField.VariableLabelSet(iron.FieldVariableTypes.V,"CauchyStress")
elasticityDerivedField.VariableLabelSet(iron.FieldVariableTypes.W,"ElasticWork")
elasticityDerivedField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfTensorComponents)
elasticityDerivedField.NumberOfComponentsSet(iron.FieldVariableTypes.V,numberOfTensorComponents)
elasticityDerivedField.NumberOfComponentsSet(iron.FieldVariableTypes.W,1)
for componentIdx in range(1,numberOfTensorComponents+1):
    elasticityDerivedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    elasticityDerivedField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,componentIdx,1)
elasticityDerivedField.ComponentMeshComponentSet(iron.FieldVariableTypes.W,1,1)
for componentIdx in range(1,numberOfTensorComponents+1):
    elasticityDerivedField.ComponentInterpolationSet(iron.FieldVariableTypes.U,componentIdx,iron.FieldInterpolationTypes.ELEMENT_BASED)
    elasticityDerivedField.ComponentInterpolationSet(iron.FieldVariableTypes.V,componentIdx,iron.FieldInterpolationTypes.ELEMENT_BASED)
elasticityDerivedField.ComponentInterpolationSet(iron.FieldVariableTypes.W,1,iron.FieldInterpolationTypes.ELEMENT_BASED)
elasticityDerivedField.CreateFinish()

# Create the derived equations set fields
elasticityEquationsSet.DerivedCreateStart(ELASTICITY_DERIVED_FIELD_USER_NUMBER,elasticityDerivedField)
elasticityEquationsSet.DerivedVariableSet(iron.EquationsSetDerivedTensorTypes.SMALL_STRAIN,iron.FieldVariableTypes.U)
elasticityEquationsSet.DerivedVariableSet(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,iron.FieldVariableTypes.V)
elasticityEquationsSet.DerivedVariableSet(iron.EquationsSetDerivedTensorTypes.ELASTIC_WORK,iron.FieldVariableTypes.W)
elasticityEquationsSet.DerivedCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

elasticityEquations = iron.Equations()
elasticityEquationsSet.EquationsCreateStart(elasticityEquations)
#elasticityEquations.SparsityTypeSet(iron.EquationsSparsityTypes.FULL)
elasticityEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
#elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
#elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.TIMING)
#elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.MATRIX)
elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.ELEMENT_MATRIX)
elasticityEquationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# PROBLEM
#-----------------------------------------------------------------------------------------------------------

elasticityProblem = iron.Problem()
elasticityProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                                  iron.ProblemTypes.LINEAR_ELASTICITY,
                                  iron.ProblemSubtypes.NONE]
elasticityProblem.CreateStart(ELASTICITY_PROBLEM_USER_NUMBER,context,elasticityProblemSpecification)
elasticityProblem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# CONTROL LOOPS
#-----------------------------------------------------------------------------------------------------------

elasticityProblem.ControlLoopCreateStart()
elasticityProblem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVER
#-----------------------------------------------------------------------------------------------------------

# Create problem solver
elasticitySolver = iron.Solver()
elasticityProblem.SolversCreateStart()
elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,elasticitySolver)
#elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
#elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.MONITOR)
#elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
#elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.TIMING)
#elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.SOLVER)
elasticitySolver.OutputTypeSet(iron.SolverOutputTypes.MATRIX)
elasticitySolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
elasticityProblem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create solver equations and add equations set to solver equations
elasticitySolver = iron.Solver()
elasticitySolverEquations = iron.SolverEquations()
elasticityProblem.SolverEquationsCreateStart()
elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,elasticitySolver)
elasticitySolver.SolverEquationsGet(elasticitySolverEquations)
#elasticitySolverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.FULL)
elasticitySolverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
elasticityEquationsSetIndex = elasticitySolverEquations.EquationsSetAdd(elasticityEquationsSet)
elasticityProblem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

elasticityBoundaryConditions = iron.BoundaryConditions()
elasticitySolverEquations.BoundaryConditionsCreateStart(elasticityBoundaryConditions)

if (numberOfDimensions == 2):
    # Set left hand edge to be built in. 
    for yNodeIdx in range(0,numberOfYNodes):
        nodeNumber = yNodeIdx*numberOfXNodes+1
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                 iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                 iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
            if (haveHermite):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)                
        if (boundaryConditionType == DIRICHLET_BCS):
            #Set downward displacement on the right hand edge 
            nodeNumber = numberOfNodes
            nodeDomain = decomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,MAX_DISPLACEMENT)
        else:
            #Set downward force on the right hand edge
            if (numberOfDimensions == 2):
                nodeNumber = numberOfNodes
                nodeDomain = decomposition.NodeDomainGet(nodeNumber,1)
                if (nodeDomain == computationalNodeNumber):
                    elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.T,1,
                                                         iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                         iron.BoundaryConditionsTypes.FIXED,0.0)
                    elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.T,1,
                                                         iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                         iron.BoundaryConditionsTypes.FIXED,MAX_FORCE)
else:
    #3D - Use analytic
    elasticitySolverEquations.BoundaryConditionsAnalytic()

elasticitySolverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVE
#-----------------------------------------------------------------------------------------------------------

elasticityProblem.Solve()

# Calculate the derived fields
elasticityEquationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTensorTypes.SMALL_STRAIN)
elasticityEquationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS)
elasticityEquationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTensorTypes.ELASTIC_WORK)

#-----------------------------------------------------------------------------------------------------------
# OUTPUT
#-----------------------------------------------------------------------------------------------------------

if(numberOfDimensions == 3):
    iron.AnalyticAnalysis_Output(elasticityDependentField,'CantileverEndLoad')

fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("LinearCantilever","FORTRAN")
fields.ElementsExport("LinearCantilever","FORTRAN")
fields.Finalise()

# Finalise OpenCMISS
iron.Finalise()
