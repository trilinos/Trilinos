INCLUDE(AddSubdirectories)
INCLUDE(PackageCopyFilesToBinaryDir)
  
ADD_SUBDIRECTORIES(
  solvers
  adapters
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(Amesos2_Superlu_Solve_Tests
  SOURCE_FILES ./matrices/arc130.mtx
               ./matrices/ash219.mtx
	       ./matrices/bcsstk01.mtx
	       ./matrices/bcsstk18.mtx
	       ./matrices/bcsstm01.mtx
	       ./matrices/beacxc.mtx
	       ./matrices/gemat12.mtx
	       ./matrices/gr_30_30.mtx
	       ./matrices/sherman3.mtx
	       ./matrices/young1c.mtx
	       ./matrices/amesos2_test_mat1.mtx
  )
