INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos)

ADD_SUBDIRECTORY(remEssBCs)

ASSERT_DEFINED(FEI_ENABLE_Teuchos)
IF(FEI_ENABLE_Teuchos)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../teuchos/src)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../teuchos/src)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    fei_test
    SOURCES fei_test.cpp
    ARGS "-d . -test unit_tests"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(FEITestRemEssBCsCopyFiles
    DEST_FILES
    pll_RemEssBCs.Aztec.parameters pll_RemEssBCs.i pll_RemEssBCs.in.1.0
    pll_RemEssBCs.in.3.0 pll_RemEssBCs.in.3.1 pll_RemEssBCs.in.3.2
    pll_RemEssBCs.input pll_RemEssBCs.soln.node.correct.1.1.0
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/remEssBCs
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/remEssBCs
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestRemEssBCsCopyFiles
    NAME fei_test_remEssBCs
    ARGS "-d ./remEssBCs -i pll_RemEssBCs.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  #PACKAGE_ADD_TEST(
  #  fei_test
  #  NAME fei_test_elemDOF_Aztec
  #  ARGS "-d ./elemDOF_Aztec -i elemDOF.i"
  #  COMM serial mpi
  #  NUM_MPI_PROCS 1-2
  #  )
  # 2009/01/09: rabartl: Above test timed out on my machine

#  PACKAGE_ADD_TEST(
#    fei_test
#    NAME fei_test_lagrange_20quad
#    ARGS "-d ./cr_lagrange_20quad -i lagrange.i"
#    COMM serial mpi
#    NUM_MPI_PROCS 1,2,4
#    )
#
#  PACKAGE_ADD_TEST(
#    fei_test
#    NAME fei_test_cr_slave_20quad
#    ARGS "-d ./cr_slave_20quad -i slave.i"
#    COMM serial mpi
#    NUM_MPI_PROCS 1,2,4
#    )
#
#  PACKAGE_ADD_TEST(
#    fei_test
#    NAME fei_test_cr_lagrange_20quad_AZLSC
#    ARGS "-d ./cr_lagrange_20quad_AZLSC -i lagrange.i"
#    COMM serial mpi
#    NUM_MPI_PROCS 1,2,4
#    )
#
#  PACKAGE_ADD_TEST(
#    fei_test
#    NAME fei_test_cr_lagrange_20quad_old
#    ARGS "-d ./cr_lagrange_20quad_old -i lagrange.i"
#    COMM serial mpi
#    NUM_MPI_PROCS 1,2,4
#    )

  #PACKAGE_ADD_TEST(
  #  fei_test
  #  NAME fei_test_cr_multifield_slave
  #  ARGS "-d ./cr_multifield_slave -i mf.i"
  #  COMM serial mpi
  #  NUM_MPI_PROCS 1-3
  #  )
  # 2009/01/09: rabartl: Above test timed out on my machine

ENDIF()
