TRIBITS_ADD_LIBRARY(
	Iocatalyst_TestDriverMain
HEADERS
        Iocatalyst_BlockMesh.h
        Iocatalyst_BlockMeshSet.h
        Iocatalyst_DatabaseIOTest.h
SOURCES
        Iocatalyst_BlockMesh.C
        Iocatalyst_BlockMeshSet.C
        Iocatalyst_DatabaseIOTest.C
        Iocatalyst_TestDriverMain.C
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_BlockMeshTest
 SOURCES
    Iocatalyst_BlockMeshTest.C
)

TRIBITS_ADD_TEST(
	Iocatalyst_BlockMeshTest
	NAME Iocatalyst_BlockMeshTest
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_BlockMeshSetTest
 SOURCES
    Iocatalyst_BlockMeshSetTest.C
)

TRIBITS_ADD_TEST(
	Iocatalyst_BlockMeshSetTest
	NAME Iocatalyst_BlockMeshSetTest
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_LoggingTest
 SOURCES
    Iocatalyst_LoggingTest.C
)

TRIBITS_ADD_TEST(
	Iocatalyst_LoggingTest
	NAME Iocatalyst_LoggingTest
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_ManagerTest
 SOURCES
    Iocatalyst_ManagerTest.C
)

TRIBITS_ADD_TEST(
	Iocatalyst_ManagerTest
	NAME Iocatalyst_ManagerTest
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_StructuredBlockTest
 SOURCES
    Iocatalyst_StructuredBlockTest.C
)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_ElementBlockTest
 SOURCES
    Iocatalyst_ElementBlockTest.C
)

SET(CATALYST_IOSS_ENTITY_TESTS
    Iocatalyst_StructuredBlockTest
    Iocatalyst_ElementBlockTest
)

FOREACH(CATALYST_TEST ${CATALYST_IOSS_ENTITY_TESTS})
    FOREACH(MPI_RANKS 1 4)
        TRIBITS_ADD_TEST(
	        ${CATALYST_TEST}
	        NAME ${CATALYST_TEST}
	        NUM_MPI_PROCS ${MPI_RANKS}
        )
    ENDFOREACH()
ENDFOREACH()

function(catalyst_test_ioshell_generated test_name
                                         ioshell_gen_args
                                         ioshell_output_file_name
                                         test_time
                                         test_time_step
                                         num_procs)

set(CATALYST_FNAME catalyst_time_${test_time_step}_np_${num_procs}_${ioshell_output_file_name})
set(IOSHELL_FNAME ioshell_time_${test_time_step}_np_${num_procs}_${ioshell_output_file_name})

if(num_procs GREATER 1)
  set(CATALYST_FNAME_EPU ${CATALYST_FNAME}.${num_procs}.0)
else()
  set(CATALYST_FNAME_EPU ${CATALYST_FNAME})
endif()

TRIBITS_ADD_ADVANCED_TEST(${test_name}
   TEST_0 EXEC io_shell ARGS --in_type generated ${ioshell_gen_args}
          ${ioshell_output_file_name}
     DIRECTORY ../main
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_1 EXEC decomp ARGS -p ${num_procs} ${ioshell_output_file_name}
     DIRECTORY ../../../../scripts
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_2 EXEC io_shell ARGS -out_type catalyst ${ioshell_output_file_name}
          catalyst.bin
     DIRECTORY ../main
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${num_procs}
   TEST_3 EXEC io_shell ARGS -in_type catalyst catalyst.bin
          ${CATALYST_FNAME}
     DIRECTORY ../main
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${num_procs}
   TEST_4 EXEC epu ARGS -auto ${CATALYST_FNAME_EPU}
     DIRECTORY ../../../../applications/epu
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_5 EXEC io_shell ARGS -select_times ${test_time}
          ${ioshell_output_file_name} ${IOSHELL_FNAME}
     DIRECTORY ../main
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_6 EXEC exodiff ARGS -pedantic ${IOSHELL_FNAME}
          ${CATALYST_FNAME}
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   ENVIRONMENT CATALYST_DATA_DUMP_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/${test_name}
     IOSS_PROPERTIES=CATALYST_READER_TIME_STEP=${test_time_step}
  )
endfunction()

foreach(NUM_PROCS 1 4)

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_hex_MPI_${NUM_PROCS}"
    "10x10x10+times:4+variables:element,2,nodal,3"
    "ioshell_10x10x10.g" "3" 3 ${NUM_PROCS})

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_tets_MPI_${NUM_PROCS}"
    "10x10x10+tets:+times:2+variables:element,2,nodal,3"
    "ioshell_10x10x10_tets.g" "1" 1 ${NUM_PROCS})

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_pyramids_MPI_${NUM_PROCS}"
    "10x10x10+pyramids:+times:2+variables:element,2,nodal,3"
    "ioshell_10x10x10_pyramids.g" "1" 1 ${NUM_PROCS})

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_shell_MPI_${NUM_PROCS}"
    "10x10x10+shell:xX:+times:2+variables:element,2,nodal,3"
    "ioshell_10x10x10_shell.g" "1" 1 ${NUM_PROCS})

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_nodeset_MPI_${NUM_PROCS}"
    "10x10x10+nodeset:xX:+times:2+variables:element,2,nodal,3,nodeset,4"
    "ioshell_10x10x10_nodeset.g" "1" 1 ${NUM_PROCS})

  catalyst_test_ioshell_generated(
    "Iocatalyst_10x10x10_sideset_MPI_${NUM_PROCS}"
    "10x10x10+sideset:xX:+times:2+variables:element,2,nodal,3,sideset,4"
    "ioshell_10x10x10_sideset.g" "1" 1 ${NUM_PROCS})

endforeach()

function(catalyst_test_ioshell_exodus_file test_name
                                           input_file
                                           test_time
                                           test_time_step
                                           num_procs)

set(CATALYST_FNAME catalyst_time_${test_time_step}_np_${num_procs}_${input_file})
set(IOSHELL_FNAME ioshell_time_${test_time_step}_np_${num_procs}_${input_file})

if(num_procs GREATER 1)
  set(CATALYST_FNAME_EPU ${CATALYST_FNAME}.${num_procs}.0)
else()
  set(CATALYST_FNAME_EPU ${CATALYST_FNAME})
endif()

TRIBITS_ADD_ADVANCED_TEST(${test_name}
    TEST_0 EXEC decomp ARGS -p
           ${num_procs} ${input_file}
      DIRECTORY ../../../../scripts
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_1 EXEC io_shell ARGS -out_type catalyst
           ${input_file} catalyst.bin
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS ${num_procs}
    TEST_2 EXEC io_shell ARGS -in_type catalyst
           catalyst.bin ${CATALYST_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS ${num_procs}
    TEST_3 EXEC epu ARGS -auto ${CATALYST_FNAME_EPU}
      DIRECTORY ../../../../applications/epu
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_4 EXEC io_shell ARGS -select_times
           ${test_time} ${input_file} ${IOSHELL_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_5 EXEC exodiff ARGS -pedantic
           ${IOSHELL_FNAME} ${CATALYST_FNAME}
      DIRECTORY ../../../../applications/exodiff
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    ENVIRONMENT CATALYST_DATA_DUMP_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/${test_name}
      IOSS_PROPERTIES=CATALYST_READER_TIME_STEP=${test_time_step}
    ADDED_TEST_NAME_OUT FileTestExodus
)

set_tests_properties(${FileTestExodus}
   PROPERTIES FIXTURES_SETUP WriteTestConduit)

endfunction()

set(EXO_TEST_FILES
    cube.g
    two-block.g
    8-block.g
    can.ex2)

foreach(fname ${EXO_TEST_FILES})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../main/test/${fname}
    ${fname} COPYONLY)
endforeach()

foreach(NUM_PROCS 1 4)
  catalyst_test_ioshell_exodus_file(
    "Iocatalyst_cube_g_MPI_${NUM_PROCS}"
    "cube.g" "1" 1 ${NUM_PROCS})
  catalyst_test_ioshell_exodus_file(
    "Iocatalyst_two_block_g_MPI_${NUM_PROCS}"
    "two-block.g" "1" 1 ${NUM_PROCS})
  catalyst_test_ioshell_exodus_file(
    "Iocatalyst_eight_block_g_MPI_${NUM_PROCS}"
    "8-block.g" "0.05" 5 ${NUM_PROCS})
  catalyst_test_ioshell_exodus_file(
    "Iocatalyst_can_ex2_MPI_${NUM_PROCS}"
    "can.ex2" "0.00219992" 22 ${NUM_PROCS})
endforeach()

if ( NETCDF_NCDUMP_BINARY )

function(catalyst_test_ioshell_generated_exodus_file test_name
                                                     gen_exo_file_command
                                                     gen_exo_file_args
                                                     gen_exo_file_name
                                                     test_time
                                                     test_time_step)

set(CATALYST_FNAME catalyst_time_${test_time_step}_${gen_exo_file_name})
set(IOSHELL_FNAME ioshell_time_${test_time_step}_${gen_exo_file_name})
set(INPUT_FILE_PATH ${gen_exo_file_name})

TRIBITS_ADD_ADVANCED_TEST(${test_name}
    TEST_0 EXEC ${gen_exo_file_command} ARGS ${gen_exo_file_args}
      DIRECTORY ../../../exodus/test
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_1 EXEC io_shell ARGS -out_type catalyst ${INPUT_FILE_PATH} catalyst.bin
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_2 EXEC io_shell ARGS -in_type catalyst catalyst.bin ${CATALYST_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_3 EXEC io_shell ARGS -select_times ${test_time} ${INPUT_FILE_PATH} ${IOSHELL_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_4 EXEC exodiff ARGS -pedantic ${IOSHELL_FNAME} ${CATALYST_FNAME}
      DIRECTORY ../../../../applications/exodiff
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    ENVIRONMENT CATALYST_DATA_DUMP_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/${test_name}
      IOSS_PROPERTIES=CATALYST_READER_TIME_STEP=${test_time_step}
)

endfunction()

catalyst_test_ioshell_generated_exodus_file("Iocatalyst_test_blob_exo_MPI_1"
  "testwt-blob" "" "test-blob.exo" "0.02" 1)

catalyst_test_ioshell_generated_exodus_file("Iocatalyst_test_seacas_exodus_exoIIC_MPI_1"
  "SEACASExodus_ExoIICTests.exe" "CreateEdgeFace" "edgeFace.exo" "1" 0)

endif()

IF (TPL_ENABLE_CGNS)

if ( CGNS_CGNSDIFF_BINARY )
   # Do nothing. Set earlier in config process...
else()
   find_program (CGNS_CGNSDIFF_BINARY
	     NAME cgnsdiff
	     PATHS
	       ${CGNS_BINARY_DIR}
               ${CGNS_LIBRARY_DIRS}/../bin
               $ENV{ACCESS}/bin
               $ENV{CGNS_DIR}/bin
               $ENV{CGNS_DIR}/cgnsdiff
               ENV PATH
	        ${CGNS_DIR}/bin
	     )
endif()

if ( CGNS_CGNSDIFF_BINARY )

function(catalyst_test_ioshell_cgns_file test_name
                                         input_file
                                         test_time
                                         test_time_step)

set(CATALYST_FNAME catalyst_time_${test_time_step}_${input_file})
set(IOSHELL_FNAME ioshell_time_${test_time_step}_${input_file})
set(INPUT_FILE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../main/test/${input_file})

TRIBITS_ADD_ADVANCED_TEST(${test_name}
    TEST_0 EXEC io_shell ARGS -debug -out_type catalyst ${INPUT_FILE_PATH} catalyst.bin
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_1 EXEC io_shell ARGS -in_type catalyst catalyst.bin ${CATALYST_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_2 EXEC io_shell ARGS -select_times ${test_time} ${INPUT_FILE_PATH} ${IOSHELL_FNAME}
      DIRECTORY ../main
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    TEST_3 EXEC ${CGNS_CGNSDIFF_BINARY} ARGS -d ${IOSHELL_FNAME} ${CATALYST_FNAME}
      NOEXEPREFIX NOEXESUFFIX
      NUM_MPI_PROCS 1
    ENVIRONMENT CATALYST_DATA_DUMP_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/${test_name}
      IOSS_PROPERTIES=CATALYST_READER_TIME_STEP=${test_time_step}
    ADDED_TEST_NAME_OUT FileTestCGNS
)

set_tests_properties(${FileTestCGNS}
   PROPERTIES FIXTURES_SETUP WriteTestConduit)

endfunction()

catalyst_test_ioshell_cgns_file(
  "Iocatalyst_sparc1_cgns_MPI_1" "sparc1.cgns" "15.992" 1)

TRIBITS_ADD_EXECUTABLE(
    Iocatalyst_ConduitReadTest
 SOURCES
    Iocatalyst_ConduitReadTest.C
)

TRIBITS_ADD_TEST(
	Iocatalyst_ConduitReadTest
  NAME Iocatalyst_ConduitReadTest
	NUM_MPI_PROCS 1
  ADDED_TESTS_NAMES_OUT ConduitReadTest
)

set_tests_properties(${ConduitReadTest}
  PROPERTIES FIXTURES_REQUIRED WriteTestConduit)

endif()
endif()
