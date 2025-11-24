ADD_SUBDIRECTORY(models)

ROL_COPY_FILES_TO_BINARY_DIR(
  ProjectAxisymmetricDataCopy
  SOURCE_FILES
    mesh/axi-diffuser-tri-clip.jou
    mesh/axi-nonconvex.jou
    mesh/axi-spherical.jou
    mesh/README.md
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)
