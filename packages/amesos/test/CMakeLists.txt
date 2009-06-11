#need to add scripts dirs

ADD_SUBDIRECTORY(Test_Basic)
ADD_SUBDIRECTORY(Test_Epetra_RowMatrix)
ADD_SUBDIRECTORY(Test_Epetra_CrsMatrix)
ADD_SUBDIRECTORY(Test_MultipleSolves)
ADD_SUBDIRECTORY(Test_Detailed)
ADD_SUBDIRECTORY(Test_LAPACK)
ADD_SUBDIRECTORY(Test_KLU)
ADD_SUBDIRECTORY(Test_UMFPACK)
ADD_SUBDIRECTORY(Test_Singular)
ADD_SUBDIRECTORY(Test_SuperLU)
ADD_SUBDIRECTORY(Test_SuperLU_DIST)
#ADD_SUBDIRECTORY(Test_DSCPACK)
ADD_SUBDIRECTORY(TestOptions)

#There appears to be a problem with the version of mumps
#we have to test against. Commenting this out for now.
#once we get a newer version of mumps installed and working
#we can revisit getting this test to work.
#ADD_SUBDIRECTORY(Test_MUMPS)

IF(${PACKAGE_NAME}_ENABLE_Performance)
  ADD_SUBDIRECTORY(Test_Performance)
ENDIF()
