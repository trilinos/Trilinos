
cmake_minimum_required (VERSION 2.6)
project(ParaViewCatalystSierraAdaptor)

SET(PVCSA_VERSION "1_0_6")
SET(PVCSA_LIBRARY_NAME ParaViewCatalystSierraAdaptor_${PVCSA_VERSION})

SET(PVCSCA_VERSION "1_0_6")
SET(PVCSCA_LIBRARY_NAME ParaViewCatalystSierraCppAdaptor_${PVCSCA_VERSION})

SET(PVCSP_VERSION "1_0_6")
SET(PVCSP_LIBRARY_NAME ParaViewCatalystSierraParser_${PVCSP_VERSION})

find_package(ParaView REQUIRED COMPONENTS vtkPVCatalyst vtkPVPythonCatalyst)
include(${PARAVIEW_USE_FILE})

add_definitions( -DSTK_BUILT_IN_SIERRA )
add_definitions( -DFORTRAN_ONE_UNDERSCORE )

set(PATH_TO_SIERRA_CODE "NOT_SET_YET" CACHE STRING "Path To Sierra Code")
set(PATH_TO_PYTHON_LIBRARIES "NOT_SET_YET" CACHE STRING "Path To directory
containing Python shared library")

include_directories("${PATH_TO_SIERRA_CODE}/stk/stk_util")
include_directories("include")
link_directories("${PATH_TO_PYTHON_LIBRARIES}")

add_library(${PVCSA_LIBRARY_NAME} SHARED ParaViewCatalystSierraAdaptor.cxx vtkExodusIIMultiBlockDataSet.cxx)
target_link_libraries(${PVCSA_LIBRARY_NAME} vtkPVCatalyst)
target_link_libraries(${PVCSA_LIBRARY_NAME} vtkPVPythonCatalyst)
target_link_libraries(${PVCSA_LIBRARY_NAME} python2.7)

add_subdirectory(CppPipe)
include_directories(CppPipe)
add_library(${PVCSCA_LIBRARY_NAME} SHARED ParaViewCatalystSierraAdaptor.cxx vtkExodusIIMultiBlockDataSet.cxx)
set_target_properties(${PVCSCA_LIBRARY_NAME} PROPERTIES COMPILE_DEFINITIONS "USE_CPP_PIPE")
target_link_libraries(${PVCSCA_LIBRARY_NAME} vtkPVCatalyst)
target_link_libraries(${PVCSCA_LIBRARY_NAME} CppPipe)

add_subdirectory(parser)
add_library(${PVCSP_LIBRARY_NAME} SHARED ParaViewCatalystSierraParser.cxx)
target_link_libraries(${PVCSP_LIBRARY_NAME} catalystparser)
