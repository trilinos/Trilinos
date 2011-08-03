INCLUDE(AddSubdirectories)
INCLUDE(PackageCopyFilesToBinaryDir)

ADD_SUBDIRECTORIES(
  solvers
  adapters
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyMatrixTestFiles
  SOURCE_FILES ./matrices/arc130.mtx
               ./matrices/orsirr_2.mtx
               ./matrices/gemat12.mtx
               ./matrices/sherman3.mtx
               ./matrices/young1c.mtx
	       ./matrices/amesos2_test_mat0.mtx
	       ./matrices/amesos2_test_mat0_complex.mtx
               ./matrices/amesos2_test_mat1.mtx
               ./matrices/amesos2_test_mat2.mtx
               ./matrices/amesos2_test_mat3.mtx
               ./matrices/amesos2_test_mat4.mtx
  )

