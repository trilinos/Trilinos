SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.C)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)

list(APPEND HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/io_info_lib_export.h")

TRIBITS_ADD_LIBRARY(
  io_info_lib
  HEADERS io_info.h info_interface.h
  SOURCES io_info.C volume.C info_interface.C
  DEPLIBS Ioss
  )

generate_export_header(io_info_lib)
target_include_directories(io_info_lib
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>")

ASSERT_DEFINED(TPL_ENABLE_CGNS)
IF (TPL_ENABLE_CGNS)
TRIBITS_ADD_EXECUTABLE(
  struc_to_unstruc
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES struc_to_unstruc.C
  INSTALLABLE
  )

TRIBITS_ADD_EXECUTABLE(
  cgns_decomp
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES cgns_decomp.C
  INSTALLABLE
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  cth_pressure_map
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES cth_pressure_map.C vector3d.C
  INSTALLABLE
  )
TRIBITS_ADD_EXECUTABLE(
  io_shell
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES io_shell.C shell_interface.C
  INSTALLABLE
  )
TRIBITS_ADD_EXECUTABLE(
  shell_to_hex
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES shell_to_hex.C vector3d.C
  INSTALLABLE
  )
TRIBITS_ADD_EXECUTABLE(
  io_info
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES io_info_main.C
  INSTALLABLE
  )

IF (${PACKAGE_NAME}_ENABLE_SEACASExodus)
TRIBITS_ADD_EXECUTABLE(
  io_modify
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES io_modify.C modify_interface.C
  INSTALLABLE
  )
endif()

TRIBITS_ADD_EXECUTABLE(
  sphgen
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES sphgen.C
  INSTALLABLE
  )

TRIBITS_ADD_EXECUTABLE(
  skinner
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES skinner.C skinner_interface.C
  INSTALLABLE
  )

if (TPL_ENABLE_MPI)
  IF (TPL_Netcdf_PARALLEL)
    set(DECOMP_ARG "--rcb")
    set(JOIN_ARG   "--compose=mpiio")
    set(PJOIN_ARG   "--compose=pnetcdf")
    set(NPROCS "1")
    set(SERIALIZE_ARG "--serialize_io_size=2")
  else()
    set(NPROCS "1")
  endif()
endif()

if (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodiff)
TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus32
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.g
       DIRECTORY ../../../../applications/exodiff
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

if (TPL_ENABLE_MPI)
  IF (TPL_Netcdf_PARALLEL)
    TRIBITS_ADD_ADVANCED_TEST(exodus_fpp_serialize
       TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${SERIALIZE_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32-ser.g
          NOEXEPREFIX NOEXESUFFIX
          NUM_MPI_PROCS ${NPROCS}
       TEST_1 EXEC io_shell ARGS ${JOIN_ARG} ${SERIALIZE_ARG} 8-block32-ser.g 8-block-32-serj.g
          NOEXEPREFIX NOEXESUFFIX
          NUM_MPI_PROCS ${NPROCS}
       TEST_2 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block-32-serj.g
          DIRECTORY ../../../../applications/exodiff
          NOEXEPREFIX NOEXESUFFIX
          NUM_MPI_PROCS 1
    COMM mpi
    XHOSTTYPE Windows
    )
  endif()

  TRIBITS_ADD_ADVANCED_TEST(create_path_fpp
     TEST_0 CMND rm ARGS -rf create_path_test_directory
     TEST_1 EXEC io_shell ARGS --in_type generated 100x100x100 create_path_test_directory/1/2/3/4/5/tmp.g
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS ${NPROCS}
     TEST_2 CMND ls ARGS -1 create_path_test_directory/1/2/3/4/5/
     )

  TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus64
     TEST_0 EXEC io_shell ARGS -64 ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.g
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS ${NPROCS}
     TEST_1 EXEC exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.g
       DIRECTORY ../../../../applications/exodiff
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS 1

    COMM mpi serial
    )

  TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus32_pnetcdf
     TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${PJOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64p.g
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS ${NPROCS}
     TEST_1 EXEC exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64p.g
       DIRECTORY ../../../../applications/exodiff
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS 1
    COMM mpi serial
    )
ENDIF()

IF (SEACASIoss_ENABLE_THREADSAFE)
  TRIBITS_ADD_EXECUTABLE(
    io_shell_ts
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES io_shell_ts.C shell_interface.C
    INSTALLABLE
    )

  if (TPL_ENABLE_MPI)
     SET(IOSHELL_ARG --compose --in_type generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4)
  else()
     SET(IOSHELL_ARG --in_type generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4)
  endif()
  TRIBITS_ADD_ADVANCED_TEST(
    io_shell_ts
    TEST_0 NOEXEPREFIX NOEXESUFFIX EXEC io_shell_ts ARGS ${IOSHELL_ARG} gen-shell-ts.g
    TEST_1 NOEXEPREFIX NOEXESUFFIX EXEC io_shell    ARGS ${IOSHELL_ARG} gen-shell.g
    TEST_2 EXEC exodiff ARGS -pedantic gen-shell-ts.g gen-shell.g
       DIRECTORY ../../../../applications/exodiff
       NOEXEPREFIX NOEXESUFFIX
       NUM_MPI_PROCS 1
    COMM mpi serial
    FINAL_PASS_REGULAR_EXPRESSION
  )
ENDIF()
ENDIF()

IF (TPL_ENABLE_Pamgen OR Trilinos_ENABLE_Pamgen)
# This test makes sure that the adjacency information (what element blocks are adjacent to what other
# element blocks and what sidesets touch which element blocks) is the same for reading an exodus mesh
# and for reading a internally generated pamgen mesh (the exodus mesh is also from the same pamgen input
# file, so the meshes should be identical). Should be able to just compare the io_info output from both
# runs, but on some platforms, we are getting some KOKKOS/CUDA/? debug information which messes up the
# comparison.  So...we just grep the output files for the word 'element' which gets the adjacency information
# that we need to ensure the same information is generated on both (and some other information also)
TRIBITS_ADD_ADVANCED_TEST(pamgen_exodus_io_info
   TEST_0 EXEC io_shell ARGS --in_type pamgen ${CMAKE_CURRENT_SOURCE_DIR}/test/pamgen.i pamgen.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_1 EXEC io_info ARGS --db_type exodus --adjacencies pamgen.e
     OUTPUT_FILE info_exodus.tmp
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_2 EXEC io_info ARGS --db_type pamgen --adjacencies ${CMAKE_CURRENT_SOURCE_DIR}/test/pamgen.i
     OUTPUT_FILE info_pamgen.tmp
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
   TEST_3 CMND grep ARGS element info_pamgen.tmp OUTPUT_FILE info_pamgen.out
   TEST_4 CMND grep ARGS element info_exodus.tmp OUTPUT_FILE info_exodus.out
   TEST_5 CMND diff ARGS info_exodus.out info_pamgen.out
  COMM mpi serial
  )
ENDIF()

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
TRIBITS_ADD_ADVANCED_TEST(exodus64_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} --64 ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d 8-block64.cgns ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.cgns.gold
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(generated32_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type generated 2x2x8+times:2+variables:nodal,2,element,3 generated32.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d generated32.cgns ${CMAKE_CURRENT_SOURCE_DIR}/test/generated.cgns.gold
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(generated64_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type generated --64 2x2x8+times:2+variables:nodal,2,element,3 generated64.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d generated64.cgns ${CMAKE_CURRENT_SOURCE_DIR}/test/generated.cgns.gold
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus32_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d 8-block32.cgns ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.cgns.gold
  COMM mpi serial
  )
endif()

IF (TPL_Netcdf_Enables_PNetcdf OR NOT TPL_ENABLE_MPI)
TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} --compare ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g"
  NAME io_shell_compare_same_database
  NOEXEPREFIX NOEXESUFFIX
  NUM_MPI_PROCS ${NPROCS}
  COMM mpi serial
 )
ENDIF()

TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} --compare ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns"
  NAME io_shell_compare_same_database_cgns
  NOEXEPREFIX NOEXESUFFIX
  COMM mpi serial
 )

TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} ${JOIN_ARG} --out_type=null --in_type=generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4 null_file"
  NAME generated_to_null
  NOEXEPREFIX NOEXESUFFIX
  COMM mpi serial
 )

IF (ENABLE_ExoNull)
TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} ${JOIN_ARG} --out_type=exonull --in_type=generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4 null_file"
  NAME generated_to_exonull
  NOEXEPREFIX NOEXESUFFIX
  COMM mpi serial
 )
ENDIF()

if (TPL_ENABLE_MPI)
  set(EXTERNAL "--compose=external")
endif()

TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} ${JOIN_ARG} ${EXTERNAL} --out_type=null --in_type=generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4 null_file"
  NAME generated_to_null_fpp
  NOEXEPREFIX NOEXESUFFIX
  COMM mpi serial
 )

IF (ENABLE_ExoNull)
TRIBITS_ADD_TEST(
  io_shell ARGS "${DECOMP_ARG} ${JOIN_ARG} ${EXTERNAL} --out_type=exonull --in_type=generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4 null_file"
  NAME generated_to_exonull_fpp
  NOEXEPREFIX NOEXESUFFIX
  COMM mpi serial
 )
ENDIF()

if (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodiff)
TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} 8-block.cgns 8-block.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 EXEC exodiff ARGS -stat -pedantic -map -dumpmap ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block.g
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_file_per_state_to_exodus
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} --file_per_state ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block-link.cgns
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} 8-block-link.cgns 8-block-rt.g
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block-rt.g
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus_single_block
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/cube.g cube.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} cube.cgns cube.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 EXEC exodiff ARGS -stat -pedantic -map ${CMAKE_CURRENT_SOURCE_DIR}/test/cube.g cube.g
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus_double_block
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/two-block.g two-block.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} two-block.cgns two-block.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 EXEC exodiff ARGS -stat -pedantic -map ${CMAKE_CURRENT_SOURCE_DIR}/test/two-block.g two-block.g
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

IF (NOT SEACASIoss_ENABLE_THREADSAFE)
TRIBITS_ADD_ADVANCED_TEST(cgns_bc_no_parents_to_exodus
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/BC_without_parents.cgns exodus_bc_parents.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_bc_parents.gold exodus_bc_parents.e
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(cgns_Q2_bc_no_parents_to_exodus
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/BC_Q2_without_parents.cgns exodus_Q2_bc_parents.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_Q2_bc_parents.gold exodus_Q2_bc_parents.e
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

IF (disabled)
# This requires a patched version of NetCDF.  Once the PR is accepted, the test will be re-enabled.
TRIBITS_ADD_ADVANCED_TEST(exodus_append
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} --Maximum_Time=0.055 ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g base.g
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --Minimum_Time=0.055 ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g append.g
     NUM_MPI_PROCS 1
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} --append_after_time 0.055 append.g base.g
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_3 EXEC exodiff ARGS -map ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g base.g
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  )
endif()

if ( CGNS_CGNSDIFF_BINARY )
TRIBITS_ADD_ADVANCED_TEST(cgns_Q2_parents_no_parent_compare_random_decomp
   TEST_0 EXEC io_shell ARGS --random ${CMAKE_CURRENT_SOURCE_DIR}/test/BC_Q2_without_parents.cgns BC_Q2_p1.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS --random ${CMAKE_CURRENT_SOURCE_DIR}/test/BC_Q2_parents.cgns BC_Q2_p2.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d BC_Q2_p1.cgns BC_Q2_p2.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  XHOSTTYPE Windows
  )
ENDIF()
ENDIF()

TRIBITS_ADD_ADVANCED_TEST(unstructured_cgns_to_exodus
   TEST_0 EXEC io_shell ARGS ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/generated.cgns exodus_generated.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_generated.gold exodus_generated.e
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi serial
  XHOSTTYPE Windows
  )

IF (TPL_ENABLE_MPI)
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_to_unstructured
   TEST_0 EXEC struc_to_unstruc ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns exodus_struc.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_struc.gold exodus_struc.e
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM mpi
  XHOSTTYPE Windows
  )
ELSE()
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_to_unstructured
   TEST_0 EXEC struc_to_unstruc ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns exodus_struc.e
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_struc.gold exodus_struc.e
     DIRECTORY ../../../../applications/exodiff
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1
  COMM serial
  XHOSTTYPE Windows
  )
ENDIF()
ENDIF()

if ( CGNS_CGNSDIFF_BINARY )
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_with_flush
   TEST_0 EXEC io_shell ARGS --flush_interval 1 ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns sparc1.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --flush_interval 1 sparc1.cgns sparc2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns sparc2.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_self_looping
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/self-looping-grid.cgns self-looping-grid.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/self-looping-grid.cgns self-looping-grid.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_zgc
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/5blocks-struc.cgns 5blocks.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS 5blocks.cgns 5blocks2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/5blocks-struc.cgns 5blocks2.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_periodic_roundtrip_fpp
   TEST_0 EXEC io_shell ARGS --compose external ${CMAKE_CURRENT_SOURCE_DIR}/test/self-looping-grid.cgns self-looping-fpp.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --external --compose external self-looping-fpp.cgns self-looping2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 EXEC io_shell ARGS --external self-looping2.cgns self-looping3.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_3 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/self-looping-grid.cgns self-looping3.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_fpp
   TEST_0 EXEC io_shell ARGS --compose external ${CMAKE_CURRENT_SOURCE_DIR}/test/5blocks-struc.cgns 5blocks.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --external --compose external 5blocks.cgns 5blocks2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 EXEC io_shell ARGS --external 5blocks2.cgns 5blocks2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_3 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/5blocks-struc.cgns 5blocks2.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_fpp2
   TEST_0 EXEC io_shell ARGS --compose external ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns sparc2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --external sparc2.cgns sparc3.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns sparc3.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_fields
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/transient-fields.cgns transient-fields-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/transient-fields.cgns transient-fields-out.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_embedded_nodeblock_fields
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/structured_with_nodeblock_fields.cgns structured_with_nodeblock_fields_out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/structured_with_nodeblock_fields.cgns structured_with_nodeblock_fields_out.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(unstructured_2d_cgns
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/circle-square-2d.cgns 2d-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/circle-square-2d.cgns 2d-out.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_fpp_vertex_cell_fields
   TEST_0 EXEC io_shell ARGS --compose external ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc_fields.cgns sparc_fields.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --external --compose external sparc_fields.cgns sparc_fields2.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 EXEC io_shell ARGS --external sparc_fields2.cgns sparc_fields_compose.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_3 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc_fields.cgns sparc_fields_compose.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_append
   TEST_0 EXEC io_shell ARGS --Maximum_Time=0.25 ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns append_fields-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --Minimum_Time=0.26 --append_after_time=0.25 ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns append_fields-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns append_fields-out.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_reordered_zones_fields
   TEST_0 EXEC io_shell ARGS --reverse --flush_interval 1 ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns multiple_zones_fields-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/multiple_zones_fields.cgns multiple_zones_fields-out.cgns
     NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_single_fpp_single_fields
   TEST_0 EXEC io_shell ARGS --compose external ${CMAKE_CURRENT_SOURCE_DIR}/test/blottner.cgns blottner-fpp.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS --external blottner-fpp.cgns blottner-single.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/blottner.cgns blottner-single.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

if (TPL_ENABLE_MPI)
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_pre_split_join
   TEST_0 EXEC io_shell ARGS --external ${CMAKE_CURRENT_SOURCE_DIR}/test/pre-split-old.cgns pre-split-old.cgns
     NUM_MPI_PROCS 4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/pre-split_cgns.gold pre-split-old.cgns
     NUM_MPI_PROCS 1
  COMM mpi
)
endif()

IF (TPL_ENABLE_Pamgen OR Trilinos_ENABLE_Pamgen)
TRIBITS_ADD_ADVANCED_TEST(pamgen_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type pamgen ${CMAKE_CURRENT_SOURCE_DIR}/test/pamgen.i pamgen.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/pamgen.cgns.gold pamgen.cgns
  COMM mpi serial
  )
ENDIF()

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_assembly_copy
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/assemblies.cgns assemblies-out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/assemblies.cgns assemblies-out.cgns
     NUM_MPI_PROCS 1
  COMM mpi
  )

if (TPL_ENABLE_MPI)
   TRIBITS_ADD_ADVANCED_TEST(structured_cgns_assembly_copy_fpp
      TEST_0 EXEC io_shell ARGS --compose=external ${CMAKE_CURRENT_SOURCE_DIR}/test/assemblies.cgns assemblies-fpp-out.cgns
        NUM_MPI_PROCS 1-4
        NOEXEPREFIX NOEXESUFFIX
      TEST_1 EXEC io_shell ARGS --external assemblies-fpp-out.cgns assemblies-fpp-join.cgns
        NUM_MPI_PROCS 1-4
        NOEXEPREFIX NOEXESUFFIX
      TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/assemblies.cgns assemblies-fpp-join.cgns
        NUM_MPI_PROCS 1
     COMM mpi
     )
endif()
endif()

ENDIF()

IF (TPL_ENABLE_MPI AND ${PROJECT_NAME}_ENABLE_Zoltan)
TRIBITS_ADD_TEST(
   io_info ARGS --config
   NAME io_info_config_has_zoltan
   NOEXEPREFIX NOEXESUFFIX
   PASS_REGULAR_EXPRESSION "Zoltan Library is Available"
   NUM_MPI_PROCS 1
   COMM mpi serial
  )
ENDIF()

IF (TPL_ENABLE_MPI AND TPL_ENABLE_ParMETIS)
TRIBITS_ADD_TEST(
   io_info ARGS --config
   NAME io_info_config_has_parmetis
   NOEXEPREFIX NOEXESUFFIX
   PASS_REGULAR_EXPRESSION "ParMetis Library Version"
   NUM_MPI_PROCS 1
   COMM mpi serial
  )
ENDIF()
