!> Main program
PROGRAM MooneyRivlinInCellMLExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_AUTO_CREATE=1 ! 1=yes   0=no
  
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_SPATIAL_COORDINATES=3
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: QUADRATIC_BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: LINEAR_BASIS_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_XI_COORDINATES=3
  INTEGER(CMISSIntg), PARAMETER :: QUADRATIC_MESH_COMPONENT_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: LINEAR_MESH_COMPONENT_NUMBER=2

  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_COMPONENTS=3

  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_COMPONENTS=3

  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_COMPONENTS=2

  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_VARIABLES=4
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS=4
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN=6

  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_GAUSS_POINTS=3

  INTEGER(CMISSIntg), PARAMETER :: CELLML_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: CELLML_MODELS_FIELD_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: CELLML_INTERMEDIATE_FIELD_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: CELLML_PARAMETERS_FIELD_USER_NUMBER=15

  !Program types

  !Program variables
  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  CHARACTER(LEN=255) :: commandArgument
  LOGICAL  :: directoryExists = .FALSE.

  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: mpiIError
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex  
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber

  INTEGER(CMISSIntg) :: nodeNumber,nodeDomain,nodeIdx
  INTEGER(CMISSIntg), ALLOCATABLE :: bottomSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: leftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: rightSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: frontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: bottomNormalXi,leftNormalXi,rightNormalXi,frontNormalXi

  INTEGER(CMISSIntg) :: dependentVariableTypes(4)
  
  INTEGER(CMISSIntg) :: mooneyRivlinModelIndex
  INTEGER(CMISSIntg) :: cellMLIndex
  
  !CMISS variables

  TYPE(cmfe_BasisType) ::  linearBasis,quadraticBasis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_CellMLType) :: cellML
  TYPE(cmfe_CellMLEquationsType) :: cellMLEquations
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: dependentField,equationsSetField,fibreField,geometricField,materialsField
  TYPE(cmfe_FieldType) :: cellMLIntermediateField,cellMLModelsField,cellMLParametersField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: linearSolver,solver
  TYPE(cmfe_SolverType) :: cellMLSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: err

  !Get the command arguments
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments == 3) THEN
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
  ELSE
    numberOfGlobalXElements=1
    numberOfGlobalYElements=1
    numberOfGlobalZElements=1
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Problem_FiniteElementCalculate"],err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  CALL cmfe_Context_RandomSeedsSet(context,9999,err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(numberOfGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)
  CALL MPI_BCAST(numberOfGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)
  CALL MPI_BCAST(numberOfGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL cmfe_CoordinateSystem_TypeSet(coordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,err)
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,NUMBER_OF_SPATIAL_COORDINATES,err)
  CALL cmfe_CoordinateSystem_OriginSet(coordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],err)
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the CS to the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"Region",err)
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_CreateFinish(region,err)

  !Define basis functions - tri-Quadratic Lagrange and tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(quadraticBasis,err)
  CALL cmfe_Basis_CreateStart(QUADRATIC_BASIS_USER_NUMBER,context,quadraticBasis,err)
  CALL cmfe_Basis_TypeSet(quadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL cmfe_Basis_NumberOfXiSet(quadraticBasis,NUMBER_OF_XI_COORDINATES,err)
  CALL cmfe_Basis_InterpolationXiSet(quadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(quadraticBasis, &
    & [NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS],err)
  CALL cmfe_Basis_CreateFinish(quadraticBasis,err)

  CALL cmfe_Basis_Initialise(linearBasis,err)
  CALL cmfe_Basis_CreateStart(LINEAR_BASIS_USER_NUMBER,context,linearBasis,err)
  CALL cmfe_Basis_TypeSet(linearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL cmfe_Basis_NumberOfXiSet(linearBasis,NUMBER_OF_XI_COORDINATES,err)
  CALL cmfe_Basis_InterpolationXiSet(linearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(linearBasis, &
    & [NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS],err)
  CALL cmfe_Basis_CreateFinish(linearBasis,err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the basis 
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[quadraticBasis,linearBasis],err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH,WIDTH,HEIGHT],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements],err)
  !Finish the creation of the generated mesh in the region
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL cmfe_Decomposition_CalculateFacesSet(decomposition,.TRUE.,err)
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry) - quadratic interpolation
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  CALL cmfe_Field_TypeSet(geometricField,CMFE_FIELD_GEOMETRIC_TYPE,err)  
  CALL cmfe_Field_NumberOfVariablesSet(geometricField,FIELD_GEOMETRY_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL cmfe_Field_NumberOfComponentsSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_GEOMETRY_NUMBER_OF_COMPONENTS,err)  
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field - linear interpolation
  CALL cmfe_Field_Initialise(fibreField,err)
  CALL cmfe_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
  CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)        
  CALL cmfe_Field_geometricFieldSet(fibreField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(fibreField,FIELD_FIBRE_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL cmfe_Field_NumberOfComponentsSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_FIBRE_NUMBER_OF_COMPONENTS,err)  
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,LINEAR_MESH_COMPONENT_NUMBER,err) 
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,LINEAR_MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,LINEAR_MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_CreateFinish(fibreField,err)

  !Create a material field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_Field_CreateStart(FIELD_MATERIAL_USER_NUMBER,region,materialsField,err)
  CALL cmfe_Field_TypeSet(materialsField,CMFE_FIELD_MATERIAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(materialsField,decomposition,err)
  CALL cmfe_Field_geometricFieldSet(materialsField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(materialsField,FIELD_MATERIAL_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_NumberOfComponentsSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_MATERIAL_NUMBER_OF_COMPONENTS,err)
  !Default is CMFE_FIELD_NODE_BASED_INTERPOLATION
  CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
  CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,err)
  !CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
  !  & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
  !CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,2, &
  !  & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
  CALL cmfe_Field_CreateFinish(materialsField,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP,err)

  !Create the dependent field with 4 variables and the respective number of components
  !   1   U_Var_Type            4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   2   DELUDELN_Var_Type     4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   3   U1_Var_Type           6 components: 6 independent components of the strain tensor (quad interpol) [independent]
  !   4   U2_Var_Type           6 components: 6 independent components of the stress tensor (quad interpol) [dependent]
  dependentVariableTypes = [CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
    & CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_U2_VARIABLE_TYPE]
  CALL cmfe_Field_Initialise(dependentField,err)
  IF(DEPENDENT_FIELD_AUTO_CREATE/=1) THEN
    CALL cmfe_Field_CreateStart(FIELD_DEPENDENT_USER_NUMBER,region,dependentField,err)
    CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GENERAL_TYPE,err)
    CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
    CALL cmfe_Field_geometricFieldSet(dependentField,geometricField,err)
    CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
    CALL cmfe_Field_NumberOfVariablesSet(dependentField,FIELD_DEPENDENT_NUMBER_OF_VARIABLES,err)
    CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
    CALL cmfe_Field_VariableTypesSet(dependentField,dependentVariableTypes,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN,err)

    !Set interpolation for the components of the field variables. 
    !   Default is Node Based Interpolation - so there is nothing to be done for U_Variable_Type and DelUDelN_Variable_Type
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,1, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,2, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,3, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,4, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,5, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,6, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)

    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,1, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,2, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,3, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,4, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,5, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,6, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)

    !Set the corresponding mesh component
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,LINEAR_MESH_COMPONENT_NUMBER,err)

    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LINEAR_MESH_COMPONENT_NUMBER, &
      & err)

    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,4,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,5,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,6,QUADRATIC_MESH_COMPONENT_NUMBER,err)

    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,4,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,5,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,6,QUADRATIC_MESH_COMPONENT_NUMBER,err)

    !CALL cmfe_Field_ScalingTypeSet(dependentField,CMFE_FIELD_UNIT_SCALING,err)
    CALL cmfe_Field_CreateFinish(dependentField,err)
  ENDIF !DEPENDENT_FIELD_AUTO_CREATE

  !Create the equations_set
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err)
  IF(DEPENDENT_FIELD_AUTO_CREATE == 1) THEN
    CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  ENDIF
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(cellML,err)
  CALL cmfe_CellML_CreateStart(CELLML_USER_NUMBER,region,cellML,err)
  !Import a Mooney-Rivlin material law from a file
  CALL cmfe_CellML_ModelImport(cellML,"inputs/mooney_rivlin.xml",mooneyRivlinModelIndex,err)
  !Now we have imported the model we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E11",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E12",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E13",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E22",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E23",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E33",err)
  !CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/c1",err)
  !CALL cmfe_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/c2",err)
  !   - to get from the CellML side
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev11",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev12",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev13",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev22",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev23",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev33",err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(cellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(cellML,err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map the strain components
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E11",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E12",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E13",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,4,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E22",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,5,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E23",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,6,CMFE_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E33",CMFE_FIELD_VALUES_SET_TYPE,err)
  !Map the material parameters
  !CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,materialsField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
  !  & mooneyRivlinModelIndex,"equations/c1",CMFE_FIELD_VALUES_SET_TYPE,err)
  !CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,materialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE, &
  !  & mooneyRivlinModelIndex,"equations/c2",CMFE_FIELD_VALUES_SET_TYPE,err)
  !Map the stress components
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev11",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev12",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev13",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev22",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,4,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev23",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,5,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev33",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,6,CMFE_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(cellML,err)

  !Create the CellML models field
  CALL cmfe_Field_Initialise(cellMLModelsField,err)
  CALL cmfe_CellML_ModelsFieldCreateStart(cellML,CELLML_MODELS_FIELD_USER_NUMBER,cellMLModelsField,err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(cellML,err)

  !Create the CellML parameters field --- will be the strain field
  CALL cmfe_Field_Initialise(cellMLParametersField,err)
  CALL cmfe_CellML_ParametersFieldCreateStart(cellML,CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField,err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(cellML,err)

  !Create the CellML intermediate field --- will be the stress field
  CALL cmfe_Field_Initialise(cellMLIntermediateField,err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(cellML,CELLML_INTERMEDIATE_FIELD_USER_NUMBER,cellMLIntermediateField,err)
  CALL cmfe_CellML_IntermediateFieldCreateFinish(cellML,err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)   

  !Initialise dependent field from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSRP, &
    & err)

  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,0.0_CMISSRP, &
    & err)

  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
    & err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,0.0_CMISSRP, &
    & err)

  !Define the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE],problem,err)
   CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  CALL cmfe_ControlLoop_TypeSet(controlLoop,CMFE_CONTROL_SIMPLE_TYPE,err)
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)
  
  !Create the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver CellML equations
  CALL cmfe_Solver_Initialise(cellMLSolver,err)
  CALL cmfe_CellMLEquations_Initialise(cellMLEquations,err)
  CALL cmfe_Problem_CellMLEquationsCreateStart(problem,err)
  CALL cmfe_Solver_NewtonCellMLSolverGet(solver,cellMLSolver,err)
  CALL cmfe_Solver_CellMLEquationsGet(cellMLSolver,cellMLEquations,err)
  CALL cmfe_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)
  CALL cmfe_Problem_CellMLEquationsCreateFinish(problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)   
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  CALL cmfe_GeneratedMesh_SurfaceGet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,bottomSurfaceNodes,bottomNormalXi, &
    & err)
  CALL cmfe_GeneratedMesh_SurfaceGet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,leftSurfaceNodes,leftNormalXi, &
    & err)
  CALL cmfe_GeneratedMesh_SurfaceGet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,rightSurfaceNodes,rightNormalXi, &
    & err)
  CALL cmfe_GeneratedMesh_SurfaceGet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,frontSurfaceNodes,frontNormalXi, &
    & err)

  !Set x=0 nodes to no x displacment in x
  DO nodeIdx=1,SIZE(leftSurfaceNodes,1)
    nodeNumber=leftSurfaceNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    ENDIF
  ENDDO
  !Set x=WIDTH nodes to 10% x displacement
  DO nodeIdx=1,SIZE(rightSurfaceNodes,1)
    nodeNumber=rightSurfaceNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,1.1_CMISSRP*WIDTH,err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO nodeIdx=1,SIZE(frontSurfaceNodes,1)
    nodeNumber=frontSurfaceNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO nodeIdx=1,SIZE(bottomSurfaceNodes,1)
    nodeNumber=bottomSurfaceNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL cmfe_Problem_Solve(problem,err)

  INQUIRE(FILE="./results",EXIST=directoryExists)
  IF(.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution  
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"./results/MooneyRivlinInCellML","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"./results/MooneyRivlinInCellML","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)

    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP

  END SUBROUTINE HandleError

END PROGRAM MooneyRivlinInCellMLExample

