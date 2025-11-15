FILE(
  GLOB convFiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  problem.xml
  jacobi.xml
  chebyshev.xml
  vCycle.xml
  wCycle.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Convergence_Laplace2D_cp SOURCE_FILES ${convFiles})
#ENDFILE
