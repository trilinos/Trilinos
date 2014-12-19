IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
   TpetraClassic_ENABLE_Thrust AND TpetraClassic_ENABLE_ThreadPool)
  ADD_SUBDIRECTORY(clone)
ENDIF()

ADD_SUBDIRECTORY(memory)
ADD_SUBDIRECTORY(multiplesolve)
ADD_SUBDIRECTORY(separatesetups)
ADD_SUBDIRECTORY(stratimikos)
ADD_SUBDIRECTORY(levelwrap)

