ADD_SUBDIRECTORY(singleFrequency)
ADD_SUBDIRECTORY(multipleFrequency)

ROL_COPY_FILES_TO_BINARY_DIR(
  OEDHelmholtzOCTRealDataCopy
  SOURCE_FILES
    meshfiles/octagon.jou
    samples/hole_points_10.txt samples/hole_weights_10.txt
    samples/hole_points_20.txt samples/hole_weights_20.txt
    samples/hole_points_5000.txt samples/hole_weights_5000.txt
    samples/points_10.txt samples/weights_10.txt
    samples/points_20.txt samples/weights_20.txt
    samples/points_5000.txt samples/weights_5000.txt
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)

