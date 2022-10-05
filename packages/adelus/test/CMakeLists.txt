IF(Adelus_ENABLE_ZCPLX OR Adelus_ENABLE_DREAL)
  ADD_SUBDIRECTORY(vector_random)
  ADD_SUBDIRECTORY(vector_random_fs)
  ADD_SUBDIRECTORY(vector_random_mc)
ENDIF()

IF(Adelus_ENABLE_Teuchos AND (Adelus_ENABLE_ZCPLX OR Adelus_ENABLE_DREAL))
  ADD_SUBDIRECTORY(perf_test)
ENDIF()

