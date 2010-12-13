
ADD_SUBDIRECTORY(simple_example)
ADD_SUBDIRECTORY(square_mesh)

IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
   ADD_SUBDIRECTORY(assembly_engine)
ENDIF()
