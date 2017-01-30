
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   ml_ns
   SOURCES
      ml_teko.cpp
   COMM mpi
   )

TRIBITS_ADD_EXECUTABLE(
   ml_single
   SOURCES
      ml_teko_single.cpp
   COMM mpi
   )

TRIBITS_ADD_EXECUTABLE(
   teko_ml_prec
   SOURCES
      teko_ml_prec.cpp
   COMM mpi
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ml_teko_copy
  SOURCE_FILES
    ml_teko.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS ml_ns
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ml_teko_copy_bigdataset
  SOURCE_FILES
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/big_dataset"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/big_dataset"
  EXEDEPS ml_ns ml_single
  )

# fetch data
SET(input_files
    small_dataset/B.mm
    small_dataset/Bt.mm
    small_dataset/C.mm
    small_dataset/F.mm
    small_dataset/mass_0.mm
    small_dataset/mass_1.mm
    small_dataset/xcoords.mm
    small_dataset/ycoords.mm
    #
    big_dataset/CompP1.serial.md5
    big_dataset/B.mm
    big_dataset/Bt.mm
    big_dataset/C.mm
    big_dataset/F.mm
    big_dataset/xcoords.mm
    big_dataset/ycoords.mm
  )
include(ExternalData)
set(ExternalData_URL_TEMPLATES
  "https://sourceforge.net/projects/trilinos-data/files/%(algo)/%(hash)"
  )
foreach(file ${input_files})
  ExternalData_Expand_Arguments(
    tekoMLFetchData
    OUT_DATA DATA{${CMAKE_CURRENT_LIST_DIR}/${file}}
    )
endforeach()
ExternalData_Add_Target(tekoMLFetchData)
