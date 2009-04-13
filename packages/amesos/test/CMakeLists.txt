ADD_SUBDIRECTORY(Test_KLU)
ADD_SUBDIRECTORY(Test_UMFPACK)
ADD_SUBDIRECTORY(Test_SuperLU)

#There appears to be a problem with the version of mumps
#we have to test against. Commenting this out for now.
#once we get a newer version of mumps installed and working
#we can revisit getting this test to work.
#ADD_SUBDIRECTORY(Test_MUMPS)
