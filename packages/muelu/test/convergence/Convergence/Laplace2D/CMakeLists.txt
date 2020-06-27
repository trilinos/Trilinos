FILE(GLOB convFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} problem.xml smoother1.xml smoother2.xml)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Convergence_Laplace2D_cp
  SOURCE_FILES ${convFiles}
  )
#ENDFILE
