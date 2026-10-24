# Python files to install
FILE(GLOB PyKrinoPyFiles *.py)

file(COPY ${PyKrinoPyFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pykrino/.)

INSTALL(FILES
  ${PyKrinoPyFiles}
  DESTINATION ${PyKrino_INSTALL_DIR})
  
ADD_SUBDIRECTORY(bindings)
