INCLUDE(AddSubdirectories)
INCLUDE(PackageCopyFilesToBinaryDir)
  
ADD_SUBDIRECTORIES(
  solvers
  adapters
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(Amesos2_Superlu_Solve_Tests
  SOURCE_FILES ./matrices/arc130.mtx
	       ./matrices/gemat12.mtx
	       ./matrices/sherman3.mtx
	       ./matrices/young1c.mtx
	       ./matrices/amesos2_test_mat1.mtx
	       ./matrices/amesos2_test_mat2.mtx
	       ./matrices/amesos2_test_mat3.mtx
	       ./matrices/amesos2_test_mat4.mtx
	       ./solvers/solvers_test.xml
  )

