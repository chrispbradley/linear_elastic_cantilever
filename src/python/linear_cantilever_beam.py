#!/usr/bin/env python
#
# This is an example script for a linear elasticity cantilever bea  using OpenCMISS calls in python.
# By Chris Bradley
#
#

import sys

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

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
    interpolationTypeXi = oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE
    numberOfNodesXi = 2
    numberOfGaussXi = 2
elif (interpolationType == QUADRATIC_LAGRANGE):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif (interpolationType == CUBIC_LAGRANGE):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE
    numberOfNodesXi = 4
    numberOfGaussXi = 4
elif (interpolationType == CUBIC_HERMITE):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.CUBIC_HERMITE
    numberOfNodesXi = 2
    numberOfGaussXi = 4
elif (interpolationType == LINEAR_SIMPLEX):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.LINEAR_SIMPLEX
    numberOfNodesXi = 2
    gaussOrder = 4
elif (interpolationType == QUADRATIC_SIMPLEX):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX
    numberOfNodesXi = 3
    gaussOrder = 4
elif (interpolationType == CUBIC_SIMPLEX):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.CUBIC_SIMPLEX
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

context = oc.Context()
context.Create(CONTEXT_USER_NUMBER)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

oc.OutputSetOn("LinearCantilever")

oc.DiagnosticsSetOn(oc.DiagnosticTypes.IN,[1,2,3,4,5],"",["BoundaryConditionsVariable_NeumannIntegrate"])

# Get the computational nodes information
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)

#-----------------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(numberOfDimensions)
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# REGION
#-----------------------------------------------------------------------------------------------------------

region = oc.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("Cantilever")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BASIS
#-----------------------------------------------------------------------------------------------------------

basis = oc.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
if (haveSimplex):
    basis.TypeSet(oc.BasisTypes.SIMPLEX)
else:
    basis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
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

generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(GENERATED_MESH_USER_NUMBER,region)
generatedMesh.TypeSet(oc.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
if (numberOfDimensions == 2):
    generatedMesh.ExtentSet([LENGTH,HEIGHT])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements])
else:
    generatedMesh.ExtentSet([LENGTH,WIDTH,HEIGHT])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements])
mesh = oc.Mesh()
generatedMesh.CreateFinish(MESH_USER_NUMBER,mesh)

#-----------------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

decomposition = oc.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.TypeSet(oc.DecompositionTypes.CALCULATED)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = oc.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

geometricField = oc.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(oc.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(oc.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
if (numberOfDimensions == 3):
    geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
# EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create linear elasiticity equations set
elasticityEquationsSetField = oc.Field()
elasticityEquationsSet = oc.EquationsSet()
if (numberOfDimensions == 2):
    elasticityEquationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
                                           oc.EquationsSetTypes.LINEAR_ELASTICITY,
                                           oc.EquationsSetSubtypes.TWO_DIMENSIONAL_PLANE_STRESS]
else:
    elasticityEquationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
                                           oc.EquationsSetTypes.LINEAR_ELASTICITY,
                                           oc.EquationsSetSubtypes.THREE_DIMENSIONAL_ISOTROPIC]
elasticityEquationsSet.CreateStart(ELASTICITY_EQUATIONS_SET_USER_NUMBER,region,geometricField,
                         elasticityEquationsSetSpecification,
                         ELASTICITY_EQUATIONS_SET_FIELD_USER_NUMBER,elasticityEquationsSetField)
elasticityEquationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET DEPENDENT
#-----------------------------------------------------------------------------------------------------------

elasticityDependentField = oc.Field()
elasticityEquationsSet.DependentCreateStart(ELASTICITY_DEPENDENT_FIELD_USER_NUMBER,elasticityDependentField)
elasticityDependentField.LabelSet("ElasticityDependent")
elasticityDependentField.VariableLabelSet(oc.FieldVariableTypes.U,"Displacement")
elasticityDependentField.VariableLabelSet(oc.FieldVariableTypes.T,"Traction")
elasticityEquationsSet.DependentCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET MATERIALS
#-----------------------------------------------------------------------------------------------------------

elasticityMaterialsField = oc.Field()
elasticityEquationsSet.MaterialsCreateStart(ELASTICITY_MATERIALS_FIELD_USER_NUMBER,elasticityMaterialsField)
elasticityMaterialsField.LabelSet("ElasticityMaterials")
elasticityMaterialsField.VariableLabelSet(oc.FieldVariableTypes.U,"Materials")
elasticityEquationsSet.MaterialsCreateFinish()    
# Initialise the analytic field values
elasticityMaterialsField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                   1,YOUNGS_MODULUS)
elasticityMaterialsField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                   2,POISSONS_RATIO)
if(numberOfDimensions==2):
    elasticityMaterialsField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                       3,THICKNESS)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET ANALYTIC
#-----------------------------------------------------------------------------------------------------------

elasticityAnalyticField = oc.Field()
if(numberOfDimensions==3):
    elasticityEquationsSet.AnalyticCreateStart(oc.EquationsSetLinearElasticityAnalyticFunctionTypes.CANTILEVER_END_LOAD,
                                               ELASTICITY_ANALYTIC_FIELD_USER_NUMBER,elasticityAnalyticField)
    elasticityAnalyticField.LabelSet("ElasticityAnalytic")
    elasticityAnalyticField.VariableLabelSet(oc.FieldVariableTypes.U,"Analytic")
    elasticityEquationsSet.AnalyticCreateFinish()    
    # Initialise the analytic field values
    elasticityAnalyticField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                      1,LENGTH)
    elasticityAnalyticField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                      2,HEIGHT)
    elasticityAnalyticField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                      3,WIDTH)
    elasticityAnalyticField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                      4,YOUNGS_MODULUS)
    elasticityAnalyticField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                      5,MAX_FORCE)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS SET DERIVED
#-----------------------------------------------------------------------------------------------------------

# Create a field for the derived field. Have three variables U - Small strain tensor, V - Cauchy stress, W - Elastic Work
if(numberOfDimensions==2):
    numberOfTensorComponents = 3
else:
    numberOfTensorComponents = 6
elasticityDerivedField = oc.Field()
elasticityDerivedField.CreateStart(ELASTICITY_DERIVED_FIELD_USER_NUMBER,region)
elasticityDerivedField.LabelSet("ElasticityDerived")
elasticityDerivedField.TypeSet(oc.FieldTypes.GENERAL)
elasticityDerivedField.DecompositionSet(decomposition)
elasticityDerivedField.GeometricFieldSet(geometricField)
elasticityDerivedField.DependentTypeSet(oc.FieldDependentTypes.DEPENDENT)
elasticityDerivedField.NumberOfVariablesSet(3)
elasticityDerivedField.VariableTypesSet([oc.FieldVariableTypes.U,oc.FieldVariableTypes.V,oc.FieldVariableTypes.W])
elasticityDerivedField.VariableLabelSet(oc.FieldVariableTypes.U,"SmallStrain")
elasticityDerivedField.VariableLabelSet(oc.FieldVariableTypes.V,"CauchyStress")
elasticityDerivedField.VariableLabelSet(oc.FieldVariableTypes.W,"ElasticWork")
elasticityDerivedField.NumberOfComponentsSet(oc.FieldVariableTypes.U,numberOfTensorComponents)
elasticityDerivedField.NumberOfComponentsSet(oc.FieldVariableTypes.V,numberOfTensorComponents)
elasticityDerivedField.NumberOfComponentsSet(oc.FieldVariableTypes.W,1)
for componentIdx in range(1,numberOfTensorComponents+1):
    elasticityDerivedField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)
    elasticityDerivedField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,componentIdx,1)
elasticityDerivedField.ComponentMeshComponentSet(oc.FieldVariableTypes.W,1,1)
for componentIdx in range(1,numberOfTensorComponents+1):
    elasticityDerivedField.ComponentInterpolationSet(oc.FieldVariableTypes.U,componentIdx,oc.FieldInterpolationTypes.ELEMENT_BASED)
    elasticityDerivedField.ComponentInterpolationSet(oc.FieldVariableTypes.V,componentIdx,oc.FieldInterpolationTypes.ELEMENT_BASED)
elasticityDerivedField.ComponentInterpolationSet(oc.FieldVariableTypes.W,1,oc.FieldInterpolationTypes.ELEMENT_BASED)
elasticityDerivedField.CreateFinish()

# Create the derived equations set fields
elasticityEquationsSet.DerivedCreateStart(ELASTICITY_DERIVED_FIELD_USER_NUMBER,elasticityDerivedField)
elasticityEquationsSet.DerivedVariableSet(oc.EquationsSetDerivedTensorTypes.SMALL_STRAIN,oc.FieldVariableTypes.U)
elasticityEquationsSet.DerivedVariableSet(oc.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,oc.FieldVariableTypes.V)
elasticityEquationsSet.DerivedVariableSet(oc.EquationsSetDerivedTensorTypes.ELASTIC_WORK,oc.FieldVariableTypes.W)
elasticityEquationsSet.DerivedCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

elasticityEquations = oc.Equations()
elasticityEquationsSet.EquationsCreateStart(elasticityEquations)
#elasticityEquations.SparsityTypeSet(oc.EquationsSparsityTypes.FULL)
elasticityEquations.SparsityTypeSet(oc.EquationsSparsityTypes.SPARSE)
#elasticityEquations.OutputTypeSet(oc.EquationsOutputTypes.NONE)
#elasticityEquations.OutputTypeSet(oc.EquationsOutputTypes.TIMING)
#elasticityEquations.OutputTypeSet(oc.EquationsOutputTypes.MATRIX)
elasticityEquations.OutputTypeSet(oc.EquationsOutputTypes.ELEMENT_MATRIX)
elasticityEquationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# PROBLEM
#-----------------------------------------------------------------------------------------------------------

elasticityProblem = oc.Problem()
elasticityProblemSpecification = [oc.ProblemClasses.ELASTICITY,
                                  oc.ProblemTypes.LINEAR_ELASTICITY,
                                  oc.ProblemSubtypes.NONE]
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
elasticitySolver = oc.Solver()
elasticityProblem.SolversCreateStart()
elasticityProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,elasticitySolver)
#elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.NONE)
#elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.MONITOR)
#elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.PROGRESS)
#elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.TIMING)
#elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.SOLVER)
elasticitySolver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
elasticitySolver.LinearTypeSet(oc.LinearSolverTypes.DIRECT)
elasticityProblem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create solver equations and add equations set to solver equations
elasticitySolver = oc.Solver()
elasticitySolverEquations = oc.SolverEquations()
elasticityProblem.SolverEquationsCreateStart()
elasticityProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,elasticitySolver)
elasticitySolver.SolverEquationsGet(elasticitySolverEquations)
#elasticitySolverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.FULL)
elasticitySolverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.SPARSE)
elasticityEquationsSetIndex = elasticitySolverEquations.EquationsSetAdd(elasticityEquationsSet)
elasticityProblem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

elasticityBoundaryConditions = oc.BoundaryConditions()
elasticitySolverEquations.BoundaryConditionsCreateStart(elasticityBoundaryConditions)

if (numberOfDimensions == 2):
    # Set left hand edge to be built in. 
    for yNodeIdx in range(0,numberOfYNodes):
        nodeNumber = yNodeIdx*numberOfXNodes+1
        nodeDomain = decomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                 oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                 oc.BoundaryConditionsTypes.FIXED,0.0)
            elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                 oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                 oc.BoundaryConditionsTypes.FIXED,0.0)
            if (haveHermite):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)                
        if (boundaryConditionType == DIRICHLET_BCS):
            #Set downward displacement on the right hand edge 
            nodeNumber = numberOfNodes
            nodeDomain = decomposition.NodeDomainGet(1,nodeNumber)
            if (nodeDomain == computationalNodeNumber):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                     oc.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,
                                                     oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                     oc.BoundaryConditionsTypes.FIXED,MAX_DISPLACEMENT)
        else:
            #Set downward force on the right hand edge
            if (numberOfDimensions == 2):
                nodeNumber = numberOfNodes
                nodeDomain = decomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.T,1,
                                                         oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                         oc.BoundaryConditionsTypes.FIXED,0.0)
                    elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.T,1,
                                                         oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                         oc.BoundaryConditionsTypes.FIXED,MAX_FORCE)
else:
    #3D - Use analytic
    elasticitySolverEquations.BoundaryConditionsAnalytic()

elasticitySolverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVE
#-----------------------------------------------------------------------------------------------------------

elasticityProblem.Solve()

# Calculate the derived fields
elasticityEquationsSet.DerivedVariableCalculate(oc.EquationsSetDerivedTensorTypes.SMALL_STRAIN)
elasticityEquationsSet.DerivedVariableCalculate(oc.EquationsSetDerivedTensorTypes.CAUCHY_STRESS)
elasticityEquationsSet.DerivedVariableCalculate(oc.EquationsSetDerivedTensorTypes.ELASTIC_WORK)

#-----------------------------------------------------------------------------------------------------------
# OUTPUT
#-----------------------------------------------------------------------------------------------------------

if(numberOfDimensions == 3):
    oc.AnalyticAnalysis_Output(elasticityDependentField,'CantileverEndLoad')

fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("LinearCantilever","FORTRAN")
fields.ElementsExport("LinearCantilever","FORTRAN")
fields.Finalise()

# Finalise OpenCMISS
oc.Finalise()
