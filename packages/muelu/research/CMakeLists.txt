# based on experimental code
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  IF(${PACKAGE_NAME}_ENABLE_BROKEN_TESTS)
    ADD_SUBDIRECTORY(q2q1)
  ENDIF()
  IF(${PACKAGE_NAME}_ENABLE_ADDITIVE_VARIANT)
    MESSAGE(STATUS "MG Additive Variant enabled")
    ADD_SUBDIRECTORY(max)
  ENDIF()
  ADD_SUBDIRECTORY(semicoarsening)
  ADD_SUBDIRECTORY(tawiesn/crada)
ENDIF()


