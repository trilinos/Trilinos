
SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.C)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)

TRIBITS_ADD_LIBRARY(
  io_info_lib
  HEADERS io_info.h info_interface.h
  SOURCES io_info.C volume.C info_interface.C
  DEPLIBS Ioss
  )

ASSERT_DEFINED(TPL_ENABLE_CGNS)
IF (TPL_ENABLE_CGNS)
TRIBITS_ADD_EXECUTABLE(
  struc_to_unstruc
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES struc_to_unstruc.C
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  cth_pressure_map
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES cth_pressure_map.C vector3d.C
  )
TRIBITS_ADD_EXECUTABLE(
  io_shell
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES io_shell.C shell_interface.C
  )
TRIBITS_ADD_EXECUTABLE(
  shell_to_hex
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES shell_to_hex.C vector3d.C
  )
TRIBITS_ADD_EXECUTABLE(
  io_info
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES io_info_main.C
  )
TRIBITS_ADD_EXECUTABLE(
  sphgen
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES sphgen.C
  )

TRIBITS_ADD_EXECUTABLE(
  skinner
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES skinner.C skinner_interface.C
  )

if (TPL_ENABLE_MPI)
  IF (TPL_Netcdf_PARALLEL)
    set(DECOMP_ARG "--rcb")
    set(JOIN_ARG   "--compose=mpiio")
    set(PJOIN_ARG   "--compose=pnetcdf")
    set(NPROCS "1")
  else()
    set(NPROCS "1")
  endif()
endif()

TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus32
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.g
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus64
   TEST_0 EXEC io_shell ARGS -64 ${DECOMP_ARG} ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.g
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus32_to_exodus32_pnetcdf
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${PJOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64p.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64p.g
  COMM mpi serial
  )

IF (SEACASIoss_ENABLE_THREADSAFE)
  TRIBITS_ADD_EXECUTABLE(
    io_shell_ts
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES io_shell_ts.C shell_interface.C
    )

  SET(IOSHELL_ARG --in_type generated 10x10x10+shell:xXyYzZ+sideset:xXyY+times:2+variables:element,2,nodal,3,sideset,4)
IF ($(TPL_ENABLE_MPI))
  TRIBITS_ADD_ADVANCED_TEST(
    io_shell_ts
    TEST_0 NOEXEPREFIX NOEXESUFFIX EXEC io_shell_ts ARGS ${IOSHELL_ARG} gen-shell-ts.g
    TEST_1 NOEXEPREFIX NOEXESUFFIX EXEC io_shell    ARGS ${IOSHELL_ARG} gen-shell.g
    TEST_2 CMND ../../../../applications/epu/epu ARGS -auto gen-shell-ts.g.4.0
    TEST_3 CMND ../../../../applications/epu/epu ARGS -auto gen-shell.g.4.0
    TEST_4 CMND ../../../../applications/exodiff/exodiff ARGS -pedantic gen-shell-ts.g gen-shell.g
    COMM mpi serial
    FINAL_PASS_REGULAR_EXPRESSION
  )
ELSE()
  TRIBITS_ADD_ADVANCED_TEST(
    io_shell_ts
    TEST_0 NOEXEPREFIX NOEXESUFFIX EXEC io_shell_ts ARGS ${IOSHELL_ARG} gen-shell-ts.g
    TEST_1 NOEXEPREFIX NOEXESUFFIX EXEC io_shell    ARGS ${IOSHELL_ARG} gen-shell.g
    TEST_4 CMND ../../../../applications/exodiff/exodiff ARGS -pedantic gen-shell-ts.g gen-shell.g
    COMM mpi serial
    FINAL_PASS_REGULAR_EXPRESSION
  )
ENDIF()
ENDIF()


install_executable(cth_pressure_map)
IF (TPL_ENABLE_CGNS)
install_executable(struc_to_unstruc)
ENDIF()
install_executable(io_shell)
install_executable(shell_to_hex)
install_executable(io_info)
install_executable(sphgen)
install_executable(skinner)

IF (TPL_ENABLE_CGNS)

if ( CGNS_CGNSLIST_BINARY )
   # Do nothing. Set earlier in config process...
else()
   find_program (CGNS_CGNSLIST_BINARY
	     NAME cgnslist
	     PATHS
	       ${CGNS_BINARY_DIR}
               ${CGNS_LIBRARY_DIRS}/../bin
               $ENV{ACCESS}/bin
               $ENV{CGNS_DIR}/bin
               $ENV{CGNS_DIR}/cgnslist
               ENV PATH
	        ${CGNS_DIR}/bin
	     )
endif()

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


if ( CGNS_CGNSLIST_BINARY )
TRIBITS_ADD_ADVANCED_TEST(exodus64_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} --64 ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block64.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ${CGNS_CGNSLIST_BINARY} ARGS -ld 8-block64.cgns
     OUTPUT_FILE exodus64_to_unstructured_cgns.tmp
   TEST_2 CMND grep ARGS -v Information exodus64_to_unstructured_cgns.tmp OUTPUT_FILE exodus64_to_unstructured_cgns.out
   TEST_3 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_to_unstructured_cgns.gold exodus64_to_unstructured_cgns.out
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(generated32_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type generated 2x2x8+times:2+variables:nodal,2,element,3 generated32.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSLIST_BINARY} ARGS -ld generated32.cgns
     OUTPUT_FILE generated32_to_unstructured_cgns.tmp
   TEST_2 CMND grep ARGS -v Information generated32_to_unstructured_cgns.tmp OUTPUT_FILE generated32_to_unstructured_cgns.out
   TEST_3 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/generated_to_unstructured_cgns.gold generated32_to_unstructured_cgns.out
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(generated64_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type generated --64 2x2x8+times:2+variables:nodal,2,element,3 generated64.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSLIST_BINARY} ARGS -ld generated64.cgns
     OUTPUT_FILE generated64_to_unstructured_cgns.tmp
   TEST_2 CMND grep ARGS -v Information generated64_to_unstructured_cgns.tmp OUTPUT_FILE generated64_to_unstructured_cgns.out
   TEST_3 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/generated_to_unstructured_cgns.gold generated64_to_unstructured_cgns.out
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus32_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block32.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ${CGNS_CGNSLIST_BINARY} ARGS -ld 8-block32.cgns
     OUTPUT_FILE exodus32_to_unstructured_cgns.tmp
   TEST_2 CMND grep ARGS -v Information exodus32_to_unstructured_cgns.tmp OUTPUT_FILE exodus32_to_unstructured_cgns.out
   TEST_3 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_to_unstructured_cgns.gold exodus32_to_unstructured_cgns.out
  COMM mpi serial
  )
endif()

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} 8-block.cgns 8-block.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic -map -dumpmap ${CMAKE_CURRENT_SOURCE_DIR}/test/8-block.g 8-block.g
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus_single_block
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/cube.g cube.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} cube.cgns cube.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic -map ${CMAKE_CURRENT_SOURCE_DIR}/test/cube.g cube.g
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(exodus_to_unstructured_cgns_to_exodus_double_block
   TEST_0 EXEC io_shell ARGS ${DECOMP_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/two-block.g two-block.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 EXEC io_shell ARGS ${JOIN_ARG} two-block.cgns two-block.g
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_2 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic -map ${CMAKE_CURRENT_SOURCE_DIR}/test/two-block.g two-block.g
  COMM mpi serial
  )

TRIBITS_ADD_ADVANCED_TEST(unstructured_cgns_to_exodus
   TEST_0 EXEC io_shell ARGS ${JOIN_ARG} ${CMAKE_CURRENT_SOURCE_DIR}/test/generated.cgns exodus_generated.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_generated.gold exodus_generated.e
  COMM mpi serial
  )

IF ($(TPL_ENABLE_MPI))
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_to_unstructured
   TEST_0 EXEC struc_to_unstruc ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns exodus_struc.e
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS ${NPROCS}
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_struc.gold exodus_struc.e
  COMM mpi
  )
ELSE()
TRIBITS_ADD_ADVANCED_TEST(structured_cgns_to_unstructured
   TEST_0 EXEC struc_to_unstruc ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/sparc1.cgns exodus_struc.e
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 CMND ../../../../applications/exodiff/exodiff ARGS -stat -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/exodus_struc.gold exodus_struc.e
  COMM serial
  )
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

TRIBITS_ADD_ADVANCED_TEST(structured_cgns_roundtrip_bump
   TEST_0 EXEC io_shell ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/bump_str.cgns tmp.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC io_shell ARGS tmp.cgns bump_str_out.cgns
     NUM_MPI_PROCS 1-4
     NOEXEPREFIX NOEXESUFFIX
   TEST_2 CMND ${CGNS_CGNSDIFF_BINARY} ARGS -d ${CMAKE_CURRENT_SOURCE_DIR}/test/bump_str.cgns bump_str_out.cgns
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

IF (TPL_ENABLE_Pamgen)
TRIBITS_ADD_ADVANCED_TEST(pamgen_to_unstructured_cgns
   TEST_0 EXEC io_shell ARGS --in_type pamgen ${CMAKE_CURRENT_SOURCE_DIR}/testpamgen.i pamgen.cgns
     NOEXEPREFIX NOEXESUFFIX
     NUM_MPI_PROCS 1-4
   TEST_1 CMND ${CGNS_CGNSLIST_BINARY} ARGS -ld pamgen.cgns
     OUTPUT_FILE pamgen_to_unstructured_cgns.out
   TEST_2 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/pamgen_to_unstructured_cgns.gold pamgen_to_unstructured_cgns.out
  COMM mpi serial
  )
ENDIF(TPL_ENABLE_Pamgen)
endif()

ENDIF(TPL_ENABLE_CGNS)
