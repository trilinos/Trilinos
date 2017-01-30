ADD_SUBDIRECTORY(StridedSolve)
ADD_SUBDIRECTORY(AddMultiplyPrecs)
ADD_SUBDIRECTORY(BuildPreconditioner)
ADD_SUBDIRECTORY(BasicExamples)

IF (ML_ENABLE_TekoSmoothers)
   ADD_SUBDIRECTORY(ML-Teko-Coupling)
ENDIF()

# fetch data
SET(input_files
  nslhs_test.mm
  nsrhs_test.mm
  )
IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  LIST(APPEND input_files
    nsjac.mm
    )
ENDIF()
include(ExternalData)
set(ExternalData_URL_TEMPLATES
  "https://sourceforge.net/projects/trilinos-data/files/%(algo)/%(hash)"
  )
foreach(file ${input_files})
  ExternalData_Expand_Arguments(
    tekoExampleFetchData
    OUT_DATA DATA{${CMAKE_CURRENT_LIST_DIR}/data/${file}}
    )
endforeach()
ExternalData_Add_Target(tekoExampleFetchData)
