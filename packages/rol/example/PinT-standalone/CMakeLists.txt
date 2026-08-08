# example/PinT-standalone/CMakeLists.txt
# Standalone (no-Trilinos) build of the PinT FD validation test under real MPI.
# Reached only when STANDALONE_ROL AND ENABLE_MPI (see example/CMakeLists.txt).

set(PINT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../PinT/parabolic-control)

add_executable(parabolic-control_PinTConstraint_FD_Test
               ${PINT_SRC_DIR}/PinTConstraint_FD_Test.cpp)

target_include_directories(parabolic-control_PinTConstraint_FD_Test PRIVATE
                           ${PINT_SRC_DIR}              # dynamicConstraint.hpp
                           ${ROL_MPI_ADAPTER_INCLUDES}) # adapters/mpi/src/{vector,function}

# rol carries all core/compatibility includes transitively (PUBLIC).
# MPI::MPI_CXX supplies mpi.h + the MPI libs. Link MPI to the EXE only, never to 'rol'.
target_link_libraries(parabolic-control_PinTConstraint_FD_Test PRIVATE
                      rol MPI::MPI_CXX)

# Stage the XML input next to the executable (FD test reads input_ex01.xml).
# NOTE: COPYONLY runs at configure time; editing the XML needs a reconfigure to re-stage.
configure_file(${PINT_SRC_DIR}/input_ex01.xml
               ${CMAKE_CURRENT_BINARY_DIR}/input_ex01.xml COPYONLY)

# ctest: run on 2 ranks via mpirun. PinTCommunicators needs a power-of-two rank count.
# Test name is hand-chosen (not build-dir-derived, since the source lives elsewhere).
if(ENABLE_TESTS)
  add_test(NAME parabolic-control_PinTConstraint_FD_Test
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:parabolic-control_PinTConstraint_FD_Test>
                   PrintItAll ${MPIEXEC_POSTFLAGS})
  set_tests_properties(parabolic-control_PinTConstraint_FD_Test PROPERTIES
                       PASS_REGULAR_EXPRESSION "TEST PASSED"
                       PROCESSORS 2)
endif()

# Manual run (power-of-two ranks only):
#   mpirun -np 2 ./example/PinT-standalone/parabolic-control_PinTConstraint_FD_Test PrintItAll
# Expect rank-0 line: 'End Result: TEST PASSED'

# --- Unit tests for the new session + shims (independent of PinTConstraint) ---
# Linking 'rol' brings the session + shim headers transitively (PUBLIC includes);
# MPI::MPI_CXX brings mpi.h. No adapter include dirs needed (no PinT headers used).
add_executable(mpi_session_shim_unittest mpi_session_shim_unittest.cpp)
target_link_libraries(mpi_session_shim_unittest PRIVATE rol MPI::MPI_CXX)

if(ENABLE_TESTS)
  # np=1 covers the single-rank session path (getNProc()==1, sum is identity);
  # np=2 covers genuine multi-rank reductions.
  foreach(np 1 2)
    add_test(NAME mpi_session_shim_unittest_np${np}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                     ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:mpi_session_shim_unittest> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(mpi_session_shim_unittest_np${np} PROPERTIES
                         PASS_REGULAR_EXPRESSION "TEST PASSED"
                         PROCESSORS ${np})
  endforeach()
endif()
