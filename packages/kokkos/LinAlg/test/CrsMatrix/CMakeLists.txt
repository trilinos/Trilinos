INCLUDE(PackageAddExecutableAndTest)   
PACKAGE_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultMultiply  
  SOURCES 
    CrsMatrix_DefaultMultiply.cpp  
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial  
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
PACKAGE_ADD_EXECUTABLE_AND_TEST( CrsMatrix_UnitTests 
  SOURCES   
    CrsMatrix_UnitTests.cpp 
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp 
  COMM serial 
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
