ADD_SUBDIRECTORY(models)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  ProjectAxisymmetricDataCopy
  SOURCE_FILES
    mesh/axi-diffuser-tri-clip.jou
    mesh/README.md
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)
