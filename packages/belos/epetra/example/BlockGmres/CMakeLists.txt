
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)


IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockGmres_Galeri_Ex
    SOURCES BlockGmresEpetraGaleriEx.cpp
    COMM serial mpi
    ARGS "--verbose --nx=10 --ortho=DGKS"
         "--verbose --nx=10 --ortho=ICGS"
         "--verbose --nx=10 --ortho=IMGS"
    STANDARD_PASS_OUTPUT
    )


ENDIF (${PACKAGE_NAME}_ENABLE_Galeri)

