INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos)

ADD_SUBDIRECTORY(remEssBCs)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../teuchos/src)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  fei_test
  SOURCES fei_test.cpp
  ARGS "-d . -test unit_tests"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

IF(FEI_ENABLE_AztecOO)

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
    NAME remEssBCs
    ARGS "-d ./remEssBCs -i pll_RemEssBCs.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestRemEssBCsCopyFiles
    NAME remEssBCs
    ARGS "-d ./remEssBCs -i pll_RemEssBCs.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(FEITestRemEssBCsTrilinosCopyFiles
    DEST_FILES
    pll_RemEssBCs.Aztec.parameters pll_RemEssBCs.i pll_RemEssBCs.in.1.0
    pll_RemEssBCs.in.3.0 pll_RemEssBCs.in.3.1 pll_RemEssBCs.in.3.2
    pll_RemEssBCs.input pll_RemEssBCs.soln.node.correct.1.1.0
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/remEssBCs_Trilinos
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/remEssBCs_Trilinos
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestRemEssBCsTrilinosCopyFiles
    NAME remEssBCs_Trilinos
    ARGS "-d ./remEssBCs_Trilinos -i pll_RemEssBCs.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestRemEssBCsTrilinosCopyFiles
    NAME remEssBCs_Trilinos
    ARGS "-d ./remEssBCs_Trilinos -i pll_RemEssBCs.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(FEITestElemDOFAztecCopyFiles
    DEST_FILES
    elemDOF.i elemDOF.in.1.0 elemDOF.in.2.0 elemDOF.in.2.1
    elemDOF.input elemDOF.soln.elem.correct.1.1.0
    elemDOF.soln.node.correct.1.1.0 parameters.Aztec
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/elemDOF_Aztec
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/elemDOF_Aztec
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestElemDOFAztecCopyFiles
    NAME elemDOF_Aztec
    ARGS "-d ./elemDOF_Aztec -i elemDOF.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestElemDOFAztecCopyFiles
    NAME elemDOF_Aztec
    ARGS "-d ./elemDOF_Aztec -i elemDOF.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(FEITestCRLagrange20QuadCopyFiles
    DEST_FILES
    lagrange.i lagrange.in.1.0 lagrange.in.2.0 lagrange.in.2.1
    lagrange.in.4.0 lagrange.in.4.1 lagrange.in.4.2 lagrange.in.4.3
    lagrange.input lagrange.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_lagrange_20quad
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_lagrange_20quad
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadCopyFiles
    NAME lagrange_20quad
    ARGS "-d ./cr_lagrange_20quad -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
 
  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadCopyFiles
    NAME lagrange_20quad
    ARGS "-d ./cr_lagrange_20quad -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )
 
  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadCopyFiles
    NAME lagrange_20quad
    ARGS "-d ./cr_lagrange_20quad -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )
 
  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestCRLagrange20QuadAZLSCCopyFiles
    DEST_FILES
    lagrange.i lagrange.in.1.0 lagrange.in.2.0 lagrange.in.2.1
    lagrange.in.4.0 lagrange.in.4.1 lagrange.in.4.2 lagrange.in.4.3
    lagrange.input lagrange.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_lagrange_20quad_AZLSC
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_lagrange_20quad_AZLSC
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadAZLSCCopyFiles
    NAME lagrange_20quad_AZ
    ARGS "-d ./cr_lagrange_20quad_AZLSC -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadAZLSCCopyFiles
    NAME lagrange_20quad_AZ
    ARGS "-d ./cr_lagrange_20quad_AZLSC -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrange20QuadAZLSCCopyFiles
    NAME lagrange_20quad_AZ
    ARGS "-d ./cr_lagrange_20quad_AZLSC -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestCRLagrage20QuadOldCopyFiles
    DEST_FILES
    lagrange.i lagrange.in.1.0 lagrange.in.2.0 lagrange.in.2.1
    lagrange.in.4.0 lagrange.in.4.1 lagrange.in.4.2 lagrange.in.4.3
    lagrange.input lagrange.soln.node.correct.1.1.0 parameters.Aztec
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_lagrange_20quad_old
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_lagrange_20quad_old
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrage20QuadOldCopyFiles
    NAME lagrange_20quad_old
    ARGS "-d ./cr_lagrange_20quad_old -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrage20QuadOldCopyFiles
    NAME lagrange_20quad_old
    ARGS "-d ./cr_lagrange_20quad_old -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRLagrage20QuadOldCopyFiles
    NAME lagrange_20quad_old
    ARGS "-d ./cr_lagrange_20quad_old -i lagrange.i"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestMultifieldCopyFiles
    DEST_FILES
    mf.i mf.in.1.0 mf.in.2.0 mf.in.2.1 mf.in.3.0 mf.in.3.1 mf.in.3.2
    mf.input mf.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/multifield
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/multifield
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldCopyFiles
    NAME multifield
    ARGS "-d ./multifield -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldCopyFiles
    NAME multifield
    ARGS "-d ./multifield -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldCopyFiles
    NAME multifield
    ARGS "-d ./multifield -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestMultifieldVBRCopyFiles
    DEST_FILES
    mf.i mf.in.1.0 mf.in.2.0 mf.in.2.1 mf.in.3.0 mf.in.3.1 mf.in.3.2
    mf.input mf.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/multifield_vbr
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/multifield_vbr
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRCopyFiles
    NAME multifield_vbr
    ARGS "-d ./multifield_vbr -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRCopyFiles
    NAME multifield_vbr
    ARGS "-d ./multifield_vbr -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRCopyFiles
    NAME multifield_vbr
    ARGS "-d ./multifield_vbr -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestMultifieldVBRAZCopyFiles
    DEST_FILES
    mf.i mf.in.1.0 mf.in.2.0 mf.in.2.1 mf.in.3.0 mf.in.3.1 mf.in.3.2
    mf.input mf.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/multifield_vbr_az
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/multifield_vbr_az
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRAZCopyFiles
    NAME multifield_vbr_az
    ARGS "-d ./multifield_vbr_az -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRAZCopyFiles
    NAME multifield_vbr_az
    ARGS "-d ./multifield_vbr_az -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestMultifieldVBRAZCopyFiles
    NAME multifield_vbr_az
    ARGS "-d ./multifield_vbr_az -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestCRMultifieldSlaveCopyFiles
    DEST_FILES
    mf.i mf.in.1.0 mf.in.2.0 mf.in.2.1 mf.in.3.0 mf.in.3.1 mf.in.3.2
    mf.input mf.soln.node.correct.1.1.0 parameters.Trilinos
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_multifield_slave
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_multifield_slave
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRMultifieldSlaveCopyFiles
    NAME multifield_slave
    ARGS "-d ./cr_multifield_slave -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRMultifieldSlaveCopyFiles
    NAME multifield_slave
    ARGS "-d ./cr_multifield_slave -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRMultifieldSlaveCopyFiles
    NAME multifield_slave
    ARGS "-d ./cr_multifield_slave -i mf.i"
    COMM serial mpi
    NUM_MPI_PROCS 3
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestCRSlave20QuadCopyFiles
    DEST_FILES
    parameters.Trilinos slave.i slave.in.1.0 slave.in.2.0 slave.in.2.1
    slave.in.4.0 slave.in.4.1 slave.in.4.2 slave.in.4.3 slave.input
    slave.soln.node.correct.1.1.0
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_slave_20quad
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_slave_20quad
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadCopyFiles
    NAME slave_20quad
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadCopyFiles
    NAME slave_20quad
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadCopyFiles
    NAME slave_20quad
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR( FEITestCRSlave20QuadOldCopyFiles
    DEST_FILES
    parameters.Trilinos slave.i slave.in.1.0 slave.in.2.0 slave.in.2.1
    slave.in.4.0 slave.in.4.1 slave.in.4.2 slave.in.4.3 slave.input
    slave.soln.node.correct.1.1.0
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cr_slave_20quad_old
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cr_slave_20quad_old
    EXEDEPS fei_test
  )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadOldCopyFiles
    NAME slave_20quad_old
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadOldCopyFiles
    NAME slave_20quad_old
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  PACKAGE_ADD_TEST(
    fei_test FEITestCRSlave20QuadOldCopyFiles
    NAME slave_20quad_old
    ARGS "-d ./cr_slave_20quad -i slave.i"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

ENDIF()

