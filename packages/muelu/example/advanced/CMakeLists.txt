IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
   KokkosClassic_ENABLE_Thrust AND KokkosClassic_ENABLE_ThreadPool)
  ADD_SUBDIRECTORY(clone)
ENDIF()

ADD_SUBDIRECTORY(memory)
ADD_SUBDIRECTORY(multiplesolve)
ADD_SUBDIRECTORY(separatesetups)
ADD_SUBDIRECTORY(stratimikos)

