
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   AddMultiply_Example
   SOURCES
      Driver.cpp
   COMM mpi
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(addmultiply_example_copyfiles
   SOURCE_FILES
      solverparams.xml
   )

# fetch data
SET(input_files
  lhs_test.mm
  rhs_test.mm
  modified.mm
  )
include(ExternalData)
set(ExternalData_URL_TEMPLATES
  "https://sourceforge.net/projects/trilinos-data/files/%(algo)/%(hash)"
  )
foreach(file ${input_files})
  ExternalData_Expand_Arguments(
    tekoAddMultiplyPrecsFetchData
    OUT_DATA DATA{${CMAKE_CURRENT_LIST_DIR}/${file}}
    )
endforeach()
ExternalData_Add_Target(tekoAddMultiplyPrecsFetchData)
