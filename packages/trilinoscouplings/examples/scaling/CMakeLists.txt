
IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
   ${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_AztecOO AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND ${PACKAGE_NAME}_ENABLE_ML AND
   ${PACKAGE_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson_NoFE
    SOURCES example_Poisson_NoFE.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson
    SOURCES example_Poisson.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_BlockMaterials
    SOURCES example_Poisson_BlockMaterials.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
   Example_Poisson_BlockMaterials
   NAME "Poisson2D_BlockMaterials"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

  TRIBITS_ADD_TEST(
   Example_Poisson_BlockMaterials
   NAME "Poisson_BlockMaterials_RTC"
   ARGS "--meshfile=Poisson_RTC.xml"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson2D
    SOURCES example_Poisson2D.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson2D_p2
    SOURCES example_Poisson2D_p2.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson2D_p2_tpetra
    SOURCES example_Poisson2D_p2_tpetra.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    COMM mpi # Fails to build in serial (see Trilinos #465)
    )

  TRIBITS_ADD_EXECUTABLE(
    Example_Maxwell
    SOURCES example_Maxwell.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
      Example_Maxwell
      NAME Example_Maxwell_ML
      ARGS "--solverName=ML"
      NUM_MPI_PROCS 1
      COMM mpi serial
      )

  IF(${PACKAGE_NAME}_ENABLE_Stratimikos)
    TRIBITS_ADD_TEST(
      Example_Maxwell
      NAME Example_Maxwell_ML_Stratimikos
      ARGS "--solverName=ML-Stratimikos"
      NUM_MPI_PROCS 1
      COMM mpi serial
      )
  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_Belos AND
      ${PACKAGE_NAME}_ENABLE_MueLu)

    IF(MueLu_ENABLE_Epetra)
      TRIBITS_ADD_TEST(
        Example_Maxwell
        NAME Example_Maxwell_MueLu
        ARGS "--solverName=MueLu"
        NUM_MPI_PROCS 1
        COMM mpi serial
        )

      TRIBITS_ADD_TEST(
        Example_Maxwell
        NAME Example_Maxwell_MueLu
        ARGS "--solverName=MueLu"
        NUM_MPI_PROCS 4
        COMM mpi
        )

      IF(${PACKAGE_NAME}_ENABLE_Stratimikos)
        TRIBITS_ADD_TEST(
          Example_Maxwell
          NAME Example_Maxwell_MueLu_Stratimikos
          ARGS "--solverName=MueLu-Stratimikos"
          NUM_MPI_PROCS 1
          COMM mpi serial
          )

        TRIBITS_ADD_TEST(
          Example_Maxwell
          NAME Example_Maxwell_MueLu_Stratimikos
          ARGS "--solverName=MueLu-Stratimikos"
          NUM_MPI_PROCS 4
          COMM mpi
          )
      ENDIF()
    ENDIF() # MueLu_ENABLE_Epetra
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyMaxwellFiles
    SOURCE_FILES Maxwell.xml Ninja.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_GradDiv
    SOURCES example_GradDiv.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_CVFEM
    SOURCES example_CVFEM.cpp
    NUM_MPI_PROCS 1
    )

  ##### Higher order family of tests #####
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson2D_pn
    SOURCES example_Poisson2D_pn.cpp
    TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson2D_pn_tpetra
    SOURCES example_Poisson2D_pn_tpetra.cpp
    TrilinosCouplings_Pamgen_Utils.cpp
    COMM mpi
    )

 TRIBITS_ADD_TEST(
   Example_Poisson2D_pn
   NAME "Poisson2D_pn=2"
   ARGS "Poisson2D.xml"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   Example_Poisson2D_pn
   NAME "Poisson2D_pn=3"
   ARGS "Poisson2D_p3.xml"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   Example_Poisson2D_pn
   NAME "Poisson2D_pn=4"
   ARGS "Poisson2D_p4.xml"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   Example_Poisson2D_pn_tpetra
   NAME "Poisson2D_pn=4_tpetra"
   ARGS "--mesh=Poisson2D_p4.xml --solver=muelu_ho.xml"
   COMM mpi
   NUM_MPI_PROCS 1
   )



  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyPoissonFiles
  SOURCE_FILES Poisson.xml Poisson2D.xml Poisson2D_p3.xml Poisson2D_p4.xml Maxwell.xml CVFEM.xml ML_ho.xml muelu_ho.xml  muelu_ho_lin.xml Poisson_RTC.xml GradDiv.xml Poisson_avatar.xml avatar.names avatar.trees bounds.data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()



IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
   ${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_AztecOO AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND ${PACKAGE_NAME}_ENABLE_ML AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND ${PACKAGE_NAME}_ENABLE_Teko)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Darcy
    SOURCES example_Darcy.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_DivLSFEM
    SOURCES example_DivLSFEM.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_CurlLSFEM
    SOURCES example_CurlLSFEM.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_StabilizatedADR
    SOURCES example_StabilizedADR.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    NUM_MPI_PROCS 1
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES ADR.xml
  	       Darcy.xml
               CurlLSFEMin.xml
               CurlLSFEMblock_in_block.xml
               DivLSFEMin.xml

               ML_nonsym.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
   ${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_AztecOO AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND ${PACKAGE_NAME}_ENABLE_ML AND
   ${PACKAGE_NAME}_ENABLE_STKIO AND ${PACKAGE_NAME}_ENABLE_STKMesh)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_STK
    SOURCES example_Poisson_stk.cpp
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyMeshFilesSTK
  SOURCE_FILES unit_cube_10int_hex.exo
               unit_cube_5int_tet.exo
	       unit_cube_int_tri.exo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()

# ik 20 Mar 2012, mfh 17 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
# Belos and Tpetra also both have required dependencies on Teuchos.
IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_NoFE_Tpetra
    SOURCES example_Poisson_NoFE_Tpetra.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )
ENDIF()

# ik 20 Mar 2012, mfh 17 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
IF(${PACKAGE_NAME}_ENABLE_AztecOO AND
   ${PACKAGE_NAME}_ENABLE_Epetra AND
   ${PACKAGE_NAME}_ENABLE_EpetraExt AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_ML AND
   ${PACKAGE_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_NoFE_Epetra
    SOURCES example_Poisson_NoFE_Epetra.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )
ENDIF()

# mfh 18 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
# Belos and Tpetra also both have required dependencies on Teuchos.
# The executable may use MueLu only if the MueLu package is enabled.
IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    IntrepidPoisson_Pamgen_Tpetra
    SOURCES IntrepidPoisson_Pamgen_Tpetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )
ENDIF()

# mfh 24 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
# Belos both also have required dependencies on Teuchos.
IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Epetra)

  # The example has an option to use MueLu, if it is enabled.  _If_
  # MueLu is enabled, the example needs MueLu to support Epetra.  It
  # is possible to turn off MueLu's support for Epetra, so we have to
  # check here whether we are allowed to use MueLu with Epetra.
  IF(${PACKAGE_NAME}_ENABLE_MueLu)
    # This variable only exists if MueLu is enabled.
    ASSERT_DEFINED(MueLu_ENABLE_Epetra)
    SET(TrilinosCouplings_ENABLE_EXAMPLE_IntrepidPoisson_Pamgen_Epetra ${MueLu_ENABLE_Epetra})
  ELSE()
    # The example can still build without MueLu.
    SET(TrilinosCouplings_ENABLE_EXAMPLE_IntrepidPoisson_Pamgen_Epetra ON)
  ENDIF()

  IF(TrilinosCouplings_ENABLE_EXAMPLE_IntrepidPoisson_Pamgen_Epetra)
    TRIBITS_ADD_EXECUTABLE(
      IntrepidPoisson_Pamgen_Epetra
      SOURCES IntrepidPoisson_Pamgen_Epetra_main.cpp
        TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
        TrilinosCouplings_EpetraIntrepidPoissonExample.cpp
        TrilinosCouplings_EpetraIntrepidPoissonExample_SolveWithBelos.cpp
        TrilinosCouplings_Pamgen_Utils.cpp
      COMM serial mpi
      )
  ENDIF()
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_AztecOO AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Epetra AND
   ${PACKAGE_NAME}_ENABLE_ML AND
   ${PACKAGE_NAME}_ENABLE_MueLu)

  # The example requires MueLu and needs MueLu to support Tpetra.  It
  # is possible to enable MueLu, yet turn off MueLu's support for
  # Epetra, so we have to check here whether we are allowed to use
  # MueLu with Epetra.  We can't check above because if MueLu is not
  # enabled, the MueLu_ENABLE_Epetra variable does not exist.

  IF(MueLu_ENABLE_Epetra)
    TRIBITS_ADD_EXECUTABLE(
      IntrepidPoisson_Pamgen_EpetraAztecOO
      SOURCES IntrepidPoisson_Pamgen_EpetraAztecOO_main.cpp
        TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
        TrilinosCouplings_EpetraIntrepidPoissonExample.cpp
        TrilinosCouplings_EpetraIntrepidPoissonExample_SolveWithAztecOO.cpp
        TrilinosCouplings_Pamgen_Utils.cpp
      COMM serial mpi
      )
  ENDIF() # MueLu_ENABLE_Epetra
ENDIF()

# Copy over example xml file for [E|T]petraIntrepidPoisson example
TRIBITS_COPY_FILES_TO_BINARY_DIR(copyIntrepidPoissonFiles
  DEST_FILES   intrepid_poisson.xml intrepid_poisson_muelu_cheby.xml intrepid_poisson_pamgen.xml
  )

# mfh 14 Feb 2013
#
# Just build this executable; don't run it as a test.
#
# Tpetra has a required dependency on Teuchos and Kokkos,
# so we don't have to name those packages explicitly here.
# This test uses the default Kokkos Node type.
IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Tpetra AND TPETRA_INST_INT_INT)
  TRIBITS_ADD_EXECUTABLE(
    epetraTpetraImportBenchmark
    SOURCES epetraTpetraImportBenchmark.cpp
    COMM serial mpi
    )
ENDIF()

# lbv 22 March 2018
#
# Just build this executable; don't run it as a test (yet).
#
IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    StructuredIntrepidPoisson_Pamgen_Tpetra
    SOURCES StructuredIntrepidPoisson_Pamgen_Tpetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_TpetraIntrepidStructuredPoissonExample.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyMueLuStructuredFiles
    SOURCE_FILES muelu_structured_aggregation.xml
                 muelu_geometric_interpolation.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDIF()

# pbo 19 June 2018
#
# Just build this executable; don't run it as a test (yet).
#
IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    HybridIntrepidPoisson3D_Pamgen_Tpetra
    SOURCES HybridIntrepidPoisson3D_Pamgen_Tpetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_TpetraIntrepidHybridPoisson3DExample.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )
  TRIBITS_ADD_EXECUTABLE(
    HybridIntrepidPoisson2D_Pamgen_Tpetra
    SOURCES HybridIntrepidPoisson2D_Pamgen_Tpetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_TpetraIntrepidHybridPoisson2DExample.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
    COMM serial mpi
    )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyMueLuHybridFiles
   SOURCE_FILES muelu_hybrid_aggregation_2D.xml muelu_hybrid_aggregation_3D.xml
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
   DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
   )
ENDIF()


# Just build this executable; don't run it as a test.
IF (${PACKAGE_NAME}_ENABLE_Pamgen AND 
    ${PACKAGE_NAME}_ENABLE_SEACASExodus AND
    ${PACKAGE_NAME}_ENABLE_SEACASNemesis)
  message("-- TrilinosCouplings: Building rtpamgen")
   TRIBITS_ADD_EXECUTABLE(
     rtpamgen
     SOURCES rtpamgen.cpp
     COMM serial mpi
     )
   TRIBITS_COPY_FILES_TO_BINARY_DIR(TestRtpamgen
     SOURCE_FILES test.rtp
     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
     DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
     )

   
  TRIBITS_ADD_TEST(
   rtpamgen
   NAME "RTpamgen_test"
   ARGS "--file=test.rtp"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )


ENDIF()