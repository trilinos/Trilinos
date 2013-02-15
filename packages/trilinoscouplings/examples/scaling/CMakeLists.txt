INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Intrepid AND Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson_NoFE
    SOURCES example_Poisson_NoFE.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson
    SOURCES example_Poisson.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson_BlockMaterials
    SOURCES example_Poisson_BlockMaterials.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson2D
    SOURCES example_Poisson2D.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Maxwell
    SOURCES example_Maxwell.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )



  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyPoissonFiles
  SOURCE_FILES Poisson.xml Poisson2D.xml Maxwell.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  ) 
ENDIF()

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Intrepid AND Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_Pamgen AND Trilinos_ENABLE_Teko)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Darcy
    SOURCES example_Darcy.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml ifpack teko
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_DivLSFEM
    SOURCES example_DivLSFEM.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_CurlLSFEM
    SOURCES example_CurlLSFEM.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_StabilizatedADR
    SOURCES example_StabilizedADR.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
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

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Intrepid AND Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_STK)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_STK
    SOURCES example_Poisson_stk.cpp
    DEPLIBS epetra epetraext aztecoo amesos ml
    ADD_DIR_TO_NAME
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyMeshFiles
  SOURCE_FILES unit_cube_10int_hex.exo
               unit_cube_5int_tet.exo
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
IF(Trilinos_ENABLE_Belos AND
   Trilinos_ENABLE_Intrepid AND
   Trilinos_ENABLE_Pamgen AND
   Trilinos_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_NoFE_Tpetra
    SOURCES example_Poisson_NoFE_Tpetra.cpp
    COMM serial mpi
    DEPLIBS belos intrepid pamgen tpetra
    ADD_DIR_TO_NAME
    )
ENDIF()

# ik 20 Mar 2012, mfh 17 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
IF(Trilinos_ENABLE_AztecOO AND
   Trilinos_ENABLE_Epetra AND
   Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Intrepid AND
   Trilinos_ENABLE_ML AND
   Trilinos_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson_NoFE_Epetra
    SOURCES example_Poisson_NoFE_Epetra.cpp
    COMM serial mpi
    DEPLIBS aztecoo epetra epetraext intrepid ml pamgen
    ADD_DIR_TO_NAME
    )
ENDIF()

# mfh 18 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
# Belos and Tpetra also both have required dependencies on Teuchos.
IF(Trilinos_ENABLE_Belos AND
   Trilinos_ENABLE_Intrepid AND
   Trilinos_ENABLE_Pamgen AND
   Trilinos_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    IntrepidPoisson_Pamgen_Tpetra
    SOURCES IntrepidPoisson_Pamgen_Tpetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample.cpp
      TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
    COMM serial mpi
    DEPLIBS belos intrepid pamgen tpetra
    ADD_DIR_TO_NAME
    )
ENDIF()

# mfh 24 May 2012
#
# Just build this executable; don't run it as a test (yet).
#
# Intrepid has required dependencies on Teuchos, Shards, and Sacado,
# so we don't have to name those packages explicitly here.
# Belos both also have required dependencies on Teuchos.
IF(Trilinos_ENABLE_Belos AND
   Trilinos_ENABLE_Intrepid AND
   Trilinos_ENABLE_Pamgen AND
   Trilinos_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE(
    IntrepidPoisson_Pamgen_Epetra
    SOURCES IntrepidPoisson_Pamgen_Epetra_main.cpp
      TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
      TrilinosCouplings_EpetraIntrepidPoissonExample.cpp
      TrilinosCouplings_EpetraIntrepidPoissonExample_SolveWithBelos.cpp
    COMM serial mpi
    DEPLIBS belos intrepid pamgen epetra
    ADD_DIR_TO_NAME
    )
ENDIF()

# mfh 14 Feb 2013
#
# Just build this executable; don't run it as a test.
#
# Tpetra has a required dependency on Teuchos and Kokkos,
# so we don't have to name those packages explicitly here.
IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    epetraTpetraImportBenchmark
    SOURCES epetraTpetraImportBenchmark.cpp
    COMM serial mpi
    DEPLIBS tpetra epetra
    ADD_DIR_TO_NAME
    )
ENDIF()

