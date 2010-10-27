INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
   ml_ns
   SOURCES 
      ml_teko.cpp
   COMM mpi
   )

PACKAGE_ADD_EXECUTABLE(
   ml_single
   SOURCES 
      ml_teko_single.cpp
   COMM mpi
   )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ml_teko_copy
  SOURCE_FILES 
    ml_teko.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS ml_ns
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ml_teko_copy_smalldataset
  SOURCE_FILES 
    B.mm
    Bt.mm
    C.mm
    F.mm
    mass_0.mm
    mass_1.mm
    xcoords.mm
    ycoords.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/small_dataset"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/small_dataset"
  EXEDEPS ml_ns ml_single
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ml_teko_copy_bigdataset
  SOURCE_FILES 
    B.mm
    Bt.mm
    C.mm
    F.mm
    xcoords.mm
    ycoords.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/big_dataset"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/big_dataset"
  EXEDEPS ml_ns ml_single
  )
