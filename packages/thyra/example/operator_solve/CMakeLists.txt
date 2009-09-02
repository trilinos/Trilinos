INCLUDE(PackageAddExecutableAndTest)


#ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
#IF (Teuchos_ENABLE_COMPLEX)
#  PACKAGE_ADD_EXECUTABLE_AND_TEST(
#    silly1DFFT_serial
#    SOURCES silly1DFFT_serial.cpp
#    ARGS --verbose
#    COMM serial
#    )
#ENDIF()
