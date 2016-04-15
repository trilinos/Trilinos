INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_Serial)
  ADD_SUBDIRECTORY(Serial)
ENDIF()

IF(Kokkos_ENABLE_OpenMP)
  ADD_SUBDIRECTORY(OpenMP)
ENDIF()

#IF(Kokkos_ENABLE_Cuda)
#  ADD_SUBDIRECTORY(Cuda)
#ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(FunctionSpaceToolsTestDataCopy
  SOURCE_FILES
    mass_HCURL_HEX_I1_FEM_01.dat
    mass_HCURL_HEX_I1_FEM_02.dat
    mass_HCURL_HEX_I1_FEM_03.dat
    mass_HCURL_HEX_I1_FEM_04.dat
    mass_HDIV_HEX_I1_FEM_01.dat
    mass_HDIV_HEX_I1_FEM_02.dat
    mass_HDIV_HEX_I1_FEM_03.dat
    mass_HDIV_HEX_I1_FEM_04.dat
    mass_TET_FEM_P1_01.dat
    mass_TET_FEM_P1_02.dat
    mass_TET_FEM_P1_03.dat
    mass_TET_FEM_P1_04.dat
    mass_fp_HCURL_HEX_I1_FEM_01.dat
    mass_fp_HCURL_HEX_I1_FEM_02.dat
    mass_fp_HCURL_HEX_I1_FEM_03.dat
    mass_fp_HCURL_HEX_I1_FEM_04.dat
    mass_fp_HDIV_HEX_I1_FEM_01.dat
    mass_fp_HDIV_HEX_I1_FEM_02.dat
    mass_fp_HDIV_HEX_I1_FEM_03.dat
    mass_fp_HDIV_HEX_I1_FEM_04.dat
    stiff_HCURL_HEX_I1_FEM_01.dat
    stiff_HCURL_HEX_I1_FEM_02.dat
    stiff_HCURL_HEX_I1_FEM_03.dat
    stiff_HCURL_HEX_I1_FEM_04.dat
    stiff_HDIV_HEX_I1_FEM_01.dat
    stiff_HDIV_HEX_I1_FEM_02.dat
    stiff_HDIV_HEX_I1_FEM_03.dat
    stiff_HDIV_HEX_I1_FEM_04.dat
    stiff_TET_FEM_P1_01.dat
    stiff_TET_FEM_P1_02.dat
    stiff_TET_FEM_P1_03.dat
    stiff_TET_FEM_P1_04.dat
    stiff_fp_HCURL_HEX_I1_FEM_01.dat
    stiff_fp_HCURL_HEX_I1_FEM_02.dat
    stiff_fp_HCURL_HEX_I1_FEM_03.dat
    stiff_fp_HCURL_HEX_I1_FEM_04.dat
    stiff_fp_HDIV_HEX_I1_FEM_01.dat
    stiff_fp_HDIV_HEX_I1_FEM_02.dat
    stiff_fp_HDIV_HEX_I1_FEM_03.dat
    stiff_fp_HDIV_HEX_I1_FEM_04.dat
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/testdata"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/testdata"
#  EXEDEPS Test_01 #Test_02
)

