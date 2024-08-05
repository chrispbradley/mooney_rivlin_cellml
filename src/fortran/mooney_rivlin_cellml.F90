!> Main program
PROGRAM MooneyRivlinInCellMLExample

  USE OpenCMISS
#ifdef WITH_F08_MPI
  USE MPI_F08
#elif WITH_F90_MPI
  USE MPI
#endif

  IMPLICIT NONE

#ifdef WITH_F77_MPI
#include "mpif.h"
#endif

  !Test program parameters

  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_AUTO_CREATE=1 ! 1=yes   0=no
  
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP

  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_SPATIAL_COORDINATES=3
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: QUADRATIC_BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: LINEAR_BASIS_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1

  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_XI_COORDINATES=3
  INTEGER(OC_Intg), PARAMETER :: QUADRATIC_MESH_COMPONENT_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: LINEAR_MESH_COMPONENT_NUMBER=2

  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_COMPONENTS=3

  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_COMPONENTS=3

  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_COMPONENTS=2

  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_VARIABLES=4
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS=4
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN=6

  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_GAUSS_POINTS=3

  INTEGER(OC_Intg), PARAMETER :: CELLML_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: CELLML_MODELS_FIELD_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: CELLML_INTERMEDIATE_FIELD_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: CELLML_PARAMETERS_FIELD_USER_NUMBER=15

  !Program types

  !Program variables
  INTEGER(OC_Intg) :: numberOfArguments,argumentLength,status
  CHARACTER(LEN=255) :: commandArgument
  LOGICAL  :: directoryExists = .FALSE.

  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(OC_Intg) :: mpiIError
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex  
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber

  INTEGER(OC_Intg) :: nodeNumber,nodeDomain,nodeIdx
  INTEGER(OC_Intg), ALLOCATABLE :: bottomSurfaceNodes(:)
  INTEGER(OC_Intg), ALLOCATABLE :: leftSurfaceNodes(:)
  INTEGER(OC_Intg), ALLOCATABLE :: rightSurfaceNodes(:)
  INTEGER(OC_Intg), ALLOCATABLE :: frontSurfaceNodes(:)
  INTEGER(OC_Intg) :: bottomNormalXi,leftNormalXi,rightNormalXi,frontNormalXi

  INTEGER(OC_Intg) :: dependentVariableTypes(4)
  
  INTEGER(OC_Intg) :: mooneyRivlinModelIndex
  INTEGER(OC_Intg) :: cellMLIndex
  
  !CMISS variables

  TYPE(OC_BasisType) ::  linearBasis,quadraticBasis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_CellMLType) :: cellML
  TYPE(OC_CellMLEquationsType) :: cellMLEquations
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: dependentField,equationsSetField,fibreField,geometricField,materialsField
  TYPE(OC_FieldType) :: cellMLIntermediateField,cellMLModelsField,cellMLParametersField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: linearSolver,solver
  TYPE(OC_SolverType) :: cellMLSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(OC_Intg) :: err

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
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  !CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Problem_FiniteElementCalculate"],err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  CALL OC_Context_RandomSeedsSet(context,9999,err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_Bcast(numberOfGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)
  CALL MPI_Bcast(numberOfGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)
  CALL MPI_Bcast(numberOfGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiIError)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL OC_CoordinateSystem_TypeSet(coordinateSystem,OC_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,err)
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,NUMBER_OF_SPATIAL_COORDINATES,err)
  CALL OC_CoordinateSystem_OriginSet(coordinateSystem,[0.0_OC_RP,0.0_OC_RP,0.0_OC_RP],err)
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the CS to the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"Region",err)
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL OC_Region_CreateFinish(region,err)

  !Define basis functions - tri-Quadratic Lagrange and tri-Linear Lagrange
  CALL OC_Basis_Initialise(quadraticBasis,err)
  CALL OC_Basis_CreateStart(QUADRATIC_BASIS_USER_NUMBER,context,quadraticBasis,err)
  CALL OC_Basis_TypeSet(quadraticBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL OC_Basis_NumberOfXiSet(quadraticBasis,NUMBER_OF_XI_COORDINATES,err)
  CALL OC_Basis_InterpolationXiSet(quadraticBasis,[OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],err)
  CALL OC_Basis_QuadratureNumberOfGaussXiSet(quadraticBasis, &
    & [NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS],err)
  CALL OC_Basis_CreateFinish(quadraticBasis,err)

  CALL OC_Basis_Initialise(linearBasis,err)
  CALL OC_Basis_CreateStart(LINEAR_BASIS_USER_NUMBER,context,linearBasis,err)
  CALL OC_Basis_TypeSet(linearBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL OC_Basis_NumberOfXiSet(linearBasis,NUMBER_OF_XI_COORDINATES,err)
  CALL OC_Basis_InterpolationXiSet(linearBasis,[OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION,OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
  CALL OC_Basis_QuadratureNumberOfGaussXiSet(linearBasis, &
    & [NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS,NUMBER_OF_GAUSS_POINTS],err)
  CALL OC_Basis_CreateFinish(linearBasis,err)

  !Start the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the basis 
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,[quadraticBasis,linearBasis],err)
  !Define the mesh on the region
  CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH,WIDTH,HEIGHT],err)
  CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements],err)
  !Finish the creation of the generated mesh in the region
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL OC_Decomposition_CalculateFacesSet(decomposition,.TRUE.,err)
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry) - quadratic interpolation
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  CALL OC_Field_TypeSet(geometricField,OC_FIELD_GEOMETRIC_TYPE,err)  
  CALL OC_Field_NumberOfVariablesSet(geometricField,FIELD_GEOMETRY_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL OC_Field_NumberOfComponentsSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,FIELD_GEOMETRY_NUMBER_OF_COMPONENTS,err)  
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field - linear interpolation
  CALL OC_Field_Initialise(fibreField,err)
  CALL OC_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
  CALL OC_Field_DecompositionSet(fibreField,decomposition,err)        
  CALL OC_Field_geometricFieldSet(fibreField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(fibreField,FIELD_FIBRE_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL OC_Field_NumberOfComponentsSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,FIELD_FIBRE_NUMBER_OF_COMPONENTS,err)  
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,1,LINEAR_MESH_COMPONENT_NUMBER,err) 
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,2,LINEAR_MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,3,LINEAR_MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_CreateFinish(fibreField,err)

  !Create a material field and attach it to the geometric field - quadratic interpolation
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_Field_CreateStart(FIELD_MATERIAL_USER_NUMBER,region,materialsField,err)
  CALL OC_Field_TypeSet(materialsField,OC_FIELD_MATERIAL_TYPE,err)
  CALL OC_Field_DecompositionSet(materialsField,decomposition,err)
  CALL OC_Field_geometricFieldSet(materialsField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(materialsField,FIELD_MATERIAL_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_NumberOfComponentsSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,FIELD_MATERIAL_NUMBER_OF_COMPONENTS,err)
  !Default is OC_FIELD_NODE_BASED_INTERPOLATION
  CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_CONSTANT_INTERPOLATION,err)
  CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,2,OC_FIELD_CONSTANT_INTERPOLATION,err)
  !CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,1, &
  !  & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
  !CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,2, &
  !  & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
  CALL OC_Field_CreateFinish(materialsField,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2.0_OC_RP,err)
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,6.0_OC_RP,err)

  !Create the dependent field with 4 variables and the respective number of components
  !   1   U_Var_Type            4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   2   DELUDELN_Var_Type     4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   3   U1_Var_Type           6 components: 6 independent components of the strain tensor (quad interpol) [independent]
  !   4   U2_Var_Type           6 components: 6 independent components of the stress tensor (quad interpol) [dependent]
  dependentVariableTypes = [OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_DELUDELN_VARIABLE_TYPE, &
    & OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_U2_VARIABLE_TYPE]
  CALL OC_Field_Initialise(dependentField,err)
  IF(DEPENDENT_FIELD_AUTO_CREATE/=1) THEN
    CALL OC_Field_CreateStart(FIELD_DEPENDENT_USER_NUMBER,region,dependentField,err)
    CALL OC_Field_TypeSet(dependentField,OC_FIELD_GENERAL_TYPE,err)
    CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
    CALL OC_Field_geometricFieldSet(dependentField,geometricField,err)
    CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
    CALL OC_Field_NumberOfVariablesSet(dependentField,FIELD_DEPENDENT_NUMBER_OF_VARIABLES,err)
    CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
    CALL OC_Field_VariableTypesSet(dependentField,dependentVariableTypes,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_DISPL_PRESS,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE, &
      & FIELD_DEPENDENT_NUMBER_OF_COMPONENTS_STRESS_STRAIN,err)

    !Set interpolation for the components of the field variables. 
    !   Default is Node Based Interpolation - so there is nothing to be done for U_Variable_Type and DelUDelN_Variable_Type
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,1, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,2, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,3, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,4, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,5, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,6, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)

    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,1, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,2, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,3, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,4, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,5, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,6, &
      & OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,err)

    !Set the corresponding mesh component
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4,LINEAR_MESH_COMPONENT_NUMBER,err)

    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER, &
      & err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4,LINEAR_MESH_COMPONENT_NUMBER, &
      & err)

    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,4,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,5,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U1_VARIABLE_TYPE,6,QUADRATIC_MESH_COMPONENT_NUMBER,err)

    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,1,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,2,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,3,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,4,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,5,QUADRATIC_MESH_COMPONENT_NUMBER,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U2_VARIABLE_TYPE,6,QUADRATIC_MESH_COMPONENT_NUMBER,err)

    !CALL OC_Field_ScalingTypeSet(dependentField,OC_FIELD_UNIT_SCALING,err)
    CALL OC_Field_CreateFinish(dependentField,err)
  ENDIF !DEPENDENT_FIELD_AUTO_CREATE

  !Create the equations_set
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  CALL OC_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err)
  IF(DEPENDENT_FIELD_AUTO_CREATE == 1) THEN
    CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  ENDIF
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the CellML environment
  CALL OC_CellML_Initialise(cellML,err)
  CALL OC_CellML_CreateStart(CELLML_USER_NUMBER,region,cellML,err)
  !Import a Mooney-Rivlin material law from a file
  !CALL OC_CellML_ModelImport(cellML,"inputs/mooney_rivlin.xml",mooneyRivlinModelIndex,err)
  CALL OC_CellML_ModelImport(cellML,"inputs/mooney_rivlin_new.cellml",mooneyRivlinModelIndex,err)
  !Now we have imported the model we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E11",err)
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E12",err)
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E13",err)
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E22",err)
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E23",err)
  CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/E33",err)
  !CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/c1",err)
  !CALL OC_CellML_VariableSetAsKnown(cellML,mooneyRivlinModelIndex,"equations/c2",err)
  !   - to get from the CellML side
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev11",err)
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev12",err)
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev13",err)
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev22",err)
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev23",err)
  CALL OC_CellML_VariableSetAsWanted(cellML,mooneyRivlinModelIndex,"equations/Tdev33",err)
  !Finish the CellML environment
  CALL OC_CellML_CreateFinish(cellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL OC_CellML_FieldMapsCreateStart(cellML,err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map the strain components
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,1,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E11",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,2,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E12",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,3,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E13",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,4,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E22",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,5,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E23",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U1_VARIABLE_TYPE,6,OC_FIELD_VALUES_SET_TYPE, &
    & mooneyRivlinModelIndex,"equations/E33",OC_FIELD_VALUES_SET_TYPE,err)
  !Map the material parameters
  !CALL OC_CellML_CreateFieldToCellMLMap(cellML,materialsField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_VALUES_SET_TYPE, &
  !  & mooneyRivlinModelIndex,"equations/c1",OC_FIELD_VALUES_SET_TYPE,err)
  !CALL OC_CellML_CreateFieldToCellMLMap(cellML,materialsField,OC_FIELD_U_VARIABLE_TYPE,2,OC_FIELD_VALUES_SET_TYPE, &
  !  & mooneyRivlinModelIndex,"equations/c2",OC_FIELD_VALUES_SET_TYPE,err)
  !Map the stress components
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev11",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,1,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev12",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,2,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev13",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,3,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev22",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,4,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev23",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,5,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,mooneyRivlinModelIndex,"equations/Tdev33",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U2_VARIABLE_TYPE,6,OC_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL OC_CellML_FieldMapsCreateFinish(cellML,err)

  !Create the CellML models field
  CALL OC_Field_Initialise(cellMLModelsField,err)
  CALL OC_CellML_ModelsFieldCreateStart(cellML,CELLML_MODELS_FIELD_USER_NUMBER,cellMLModelsField,err)
  CALL OC_CellML_ModelsFieldCreateFinish(cellML,err)

  !Create the CellML parameters field --- will be the strain field
  CALL OC_Field_Initialise(cellMLParametersField,err)
  CALL OC_CellML_ParametersFieldCreateStart(cellML,CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField,err)
  CALL OC_CellML_ParametersFieldCreateFinish(cellML,err)

  !Create the CellML intermediate field --- will be the stress field
  CALL OC_Field_Initialise(cellMLIntermediateField,err)
  CALL OC_CellML_IntermediateFieldCreateStart(cellML,CELLML_INTERMEDIATE_FIELD_USER_NUMBER,cellMLIntermediateField,err)
  CALL OC_CellML_IntermediateFieldCreateFinish(cellML,err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)   

  !Initialise dependent field from undeformed geometry and set hydrostatic pressure
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,-8.0_OC_RP, &
    & err)

  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,5,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U1_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,6,0.0_OC_RP, &
    & err)

  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,5,0.0_OC_RP, &
    & err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U2_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,6,0.0_OC_RP, &
    & err)

  !Define the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE],problem,err)
   CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  CALL OC_ControlLoop_TypeSet(controlLoop,OC_CONTROL_SIMPLE_TYPE,err)
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)
  
  !Create the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_NewtonSolutionToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_NewtonRelativeToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver CellML equations
  CALL OC_Solver_Initialise(cellMLSolver,err)
  CALL OC_CellMLEquations_Initialise(cellMLEquations,err)
  CALL OC_Problem_CellMLEquationsCreateStart(problem,err)
  CALL OC_Solver_NewtonCellMLSolverGet(solver,cellMLSolver,err)
  CALL OC_Solver_CellMLEquationsGet(cellMLSolver,cellMLEquations,err)
  CALL OC_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)
  CALL OC_Problem_CellMLEquationsCreateFinish(problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)   
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  CALL OC_GeneratedMesh_SurfaceGet(generatedMesh,OC_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,bottomSurfaceNodes,bottomNormalXi, &
    & err)
  CALL OC_GeneratedMesh_SurfaceGet(generatedMesh,OC_GENERATED_MESH_REGULAR_LEFT_SURFACE,leftSurfaceNodes,leftNormalXi, &
    & err)
  CALL OC_GeneratedMesh_SurfaceGet(generatedMesh,OC_GENERATED_MESH_REGULAR_RIGHT_SURFACE,rightSurfaceNodes,rightNormalXi, &
    & err)
  CALL OC_GeneratedMesh_SurfaceGet(generatedMesh,OC_GENERATED_MESH_REGULAR_FRONT_SURFACE,frontSurfaceNodes,frontNormalXi, &
    & err)

  !Set x=0 nodes to no x displacment in x
  DO nodeIdx=1,SIZE(leftSurfaceNodes,1)
    nodeNumber=leftSurfaceNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
        & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    ENDIF
  ENDDO
  !Set x=WIDTH nodes to 10% x displacement
  DO nodeIdx=1,SIZE(rightSurfaceNodes,1)
    nodeNumber=rightSurfaceNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
        & OC_BOUNDARY_CONDITION_FIXED,1.1_OC_RP*WIDTH,err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO nodeIdx=1,SIZE(frontSurfaceNodes,1)
    nodeNumber=frontSurfaceNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
        & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO nodeIdx=1,SIZE(bottomSurfaceNodes,1)
    nodeNumber=bottomSurfaceNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
        & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    ENDIF
  ENDDO

  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL OC_Problem_Solve(problem,err)

  INQUIRE(FILE="./results",EXIST=directoryExists)
  IF(.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution  
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,"./results/MooneyRivlinInCellML","FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,"./results/MooneyRivlinInCellML","FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)

    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP

  END SUBROUTINE HandleError

END PROGRAM MooneyRivlinInCellMLExample

