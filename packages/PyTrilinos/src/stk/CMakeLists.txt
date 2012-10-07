# the stk project
PROJECT(PyPercept)

# we need python to build
FIND_PACKAGE(PythonLibs 2.2 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

FIND_PATH(MPI4PY_INCLUDE_DIR mpi4py HINTS ${PYTHON_INCLUDE_DIRS}/../../lib/python2.7/site-packages/mpi4py/include DOC "path to mpi4py include directory")

GET_FILENAME_COMPONENT(NETCDF_LIB_DIR ${TPL_Netcdf_LIBRARIES} PATH)
GET_FILENAME_COMPONENT(MPI_LIB_DIR ${MPI_EXEC} PATH)
STRING(REGEX REPLACE "bin$" "lib" MPI_LIB_DIR ${MPI_LIB_DIR})

#include STK related directories'
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}
          ${Shards_SOURCE_DIR}/src
          ${PyTrilinos_SOURCE_DIR}/doc/Doxygen
          ${Intrepid_SOURCE_DIR}/src/Shared
          ${Intrepid_SOURCE_DIR}/src
          ${Intrepid_SOURCE_DIR}/src/Cell
          ${Intrepid_SOURCE_DIR}/src/Discretization/Basis
          ${Intrepid_SOURCE_DIR}/src/Discretization/FunctionSpaceTools
          ${Intrepid_SOURCE_DIR}/src/Discretization/Integration
          ${Intrepid_BINARY_DIR}/src
          ${SEACAS_SOURCE_DIR}/libraries/ioss/src
          ${MPI4PY_INCLUDE_DIR}
          ${STK_SOURCE_DIR}/stk_adapt
          ${STK_SOURCE_DIR}/stk_adapt/sierra_element
          ${STK_SOURCE_DIR}/stk_search
          ${STK_SOURCE_DIR}/stk_expreval
          ${STK_SOURCE_DIR}/stk_io
          ${STK_SOURCE_DIR}/stk_io/stk_io/util
          ${STK_SOURCE_DIR}/stk_mesh
          ${STK_SOURCE_DIR}/stk_util
          ${STK_SOURCE_DIR}/stk_util/diag
          ${STK_SOURCE_DIR}/stk_util/stk_util/parallel
          ${STK_SOURCE_DIR}/stk_percept
          ${STK_SOURCE_DIR}/stk_percept/stk_percept
          ${STK_SOURCE_DIR}/stk_percept/stk_percept/norm
          ${STK_SOURCE_DIR}/stk_percept/stk_percept/fixtures
          ${STK_SOURCE_DIR}/stk_percept/stk_percept/function
          ${STK_SOURCE_DIR}/stk_percept/stk_percept/function/internal
          ${ThreadPool_SOURCE_DIR}/src
          ${STK_BINARY_DIR}/stk_percept/stk_percept
          ${STK_BINARY_DIR}/stk_util/stk_util
          ${Trios_BINARY_DIR}/ioss/src
          ${NumPy_INCLUDE_DIR}
          ${Boost_INCLUDE_DIRS}
          ${Boost_INCLUDE_DIRS}/..
          )

SET(PERCEPT_LIBS
   stk_percept
   stk_adapt
   stk_mesh_base
   stk_expreval
   stk_mesh_fem
   stk_io_util
   stk_util_util
   stk_util_diag
   stk_search
   intrepid
   shards
   teuchos
)

#
# Copy the pyfragments.swg file to the build directory
#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../pyfragments.swg.in
#  ${CMAKE_CURRENT_BINARY_DIR}/pyfragments.swg)

# create the swig custom command
SET(DOT_I_DEPENDS
    Norm.hpp
    Teuchos_RCPDecl.hpp
    Teuchos_RCP.hpp
    MDArray.hpp
    Name.hpp
    GenericFunction.hpp
    Function.hpp
    FieldFunction.hpp
    StringFunction.hpp
    Part.hpp
    FieldBase.hpp
    Types.hpp
    FieldRepository.hpp
    Gmesh_STKmesh_Fixture.hpp
    CellTopology.hpp
    FEMMetaData.hpp
    UniformRefiner.hpp
    UniformRefinerPattern.hpp
    Intrepid_FieldContainer.hpp
    PerceptMesh.hpp)

# Create the custom swig command.
# set PerceptMesh.i to generate c++ code
#SET(CMAKE_SWIG_FLAGS "-Wall")

SWIG_MODULE_GET_OUTDIR_AND_MODULE(PerceptMesh.i MODULE_OUTDIR MODULE_NAME)
SET_SOURCE_FILES_PROPERTIES(PerceptMesh.i PROPERTIES CPLUSPLUS ON OBJECT_DEPENDS "${DOT_I_DEPENDS}" )
SWIG_ADD_MODULE(PyPercept python PerceptMesh.i "${MODULE_OUTDIR}" "${MODULE_NAME}" )
SWIG_LINK_LIBRARIES(PyPercept ${PYTHON_LIBRARIES} ${PERCEPT_LIBS})

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/PyPercept/PerceptMesh.py DESTINATION site-packages)
INSTALL(TARGETS PyPercept DESTINATION site-packages)

set(dirs
  ${NETCDF_LIB_DIR}
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/nemesis
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src/exodusII
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src/generated
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src/init
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src/heartbeat
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/ioss/src/transform
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/exodus/cbind
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/chaco
  ${CMAKE_BINARY_DIR}/packages/seacas/libraries/aprepro_lib
  ${CMAKE_BINARY_DIR}/packages/intrepid/src
  ${CMAKE_BINARY_DIR}/packages/shards/src
  ${CMAKE_BINARY_DIR}/packages/sacado/src
  ${CMAKE_BINARY_DIR}/packages/ThreadPool/src
  ${CMAKE_BINARY_DIR}/packages/teuchos/src
  ${CMAKE_BINARY_DIR}/packages/stk/stk_util/stk_util/util
  ${CMAKE_BINARY_DIR}/packages/stk/stk_util/stk_util/use_cases
  ${CMAKE_BINARY_DIR}/packages/stk/stk_util/stk_util/diag
  ${CMAKE_BINARY_DIR}/packages/stk/stk_search/stk_search
  ${CMAKE_BINARY_DIR}/packages/stk/stk_percept/stk_percept
  ${CMAKE_BINARY_DIR}/packages/stk/stk_mesh/stk_mesh/fem
  ${CMAKE_BINARY_DIR}/packages/stk/stk_mesh/stk_mesh/base
  ${CMAKE_BINARY_DIR}/packages/stk/stk_io/stk_io
  ${CMAKE_BINARY_DIR}/packages/stk/stk_io/stk_io/util
  ${CMAKE_BINARY_DIR}/packages/stk/stk_expreval/stk_expreval
  ${CMAKE_BINARY_DIR}/packages/stk/stk_adapt/stk_adapt
  ${BoostLib_LIBRARY_DIRS}
  ${MPI_LIB_DIR}
)

INSTALL(CODE "
 include(${CMAKE_CURRENT_SOURCE_DIR}/FixupPyModule.cmake)
 set(percept_mesh_module \${CMAKE_INSTALL_PREFIX}/site-packages/${SWIG_MODULE_PyPercept_REAL_NAME}.so)
 fixup_py_module(\"\${CMAKE_INSTALL_PREFIX}/site-packages\" \"\${percept_mesh_module}\" \"${dirs}\")
")

