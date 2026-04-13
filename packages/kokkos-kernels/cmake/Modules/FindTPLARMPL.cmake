# Both the armpl_mp and armpl libraries define the same public symbol names.
# In order to link against the openmp  armpl symbols, instruct cmake to link against armpl_mp.
# In order to link against the default armpl symbols, instruct cmake to link against armpl.
if(KOKKOSKERNELS_INST_EXECSPACE_OPENMP)
  set(ARMPL_LIB armpl_mp)
else()
  set(ARMPL_LIB armpl)
endif()

if(ARMPL_LIBRARY_DIRS AND ARMPL_LIBRARIES)
  kokkoskernels_find_imported(ARMPL INTERFACE LIBRARIES ${ARMPL_LIBRARIES} LIBRARY_PATHS ${ARMPL_LIBRARY_DIRS})
elseif(ARMPL_LIBRARIES)
  kokkoskernels_find_imported(ARMPL INTERFACE LIBRARIES ${ARMPL_LIBRARIES})
elseif(ARMPL_LIBRARY_DIRS)
  kokkoskernels_find_imported(ARMPL INTERFACE LIBRARIES amath ${ARMPL_LIB} LIBRARY_PATHS ${ARMPL_LIBRARY_DIRS})
elseif(DEFINED ENV{ARMPL_DIR})
  set(ARMPL_BUILD $ENV{ARMPL_BUILD})
  set(ARMPL_ROOT $ENV{ARMPL_DIR})
  kokkoskernels_find_imported(ARMPL INTERFACE
    LIBRARIES     amath ${ARMPL_LIB}
    LIBRARY_PATHS ${ARMPL_ROOT}/lib
    HEADERS       armpl.h
    HEADER_PATHS  ${ARMPL_ROOT}/include)
else()
  find_package(ARMPL REQUIRED)
  kokkoskernels_create_imported_tpl(ARMPL INTERFACE LINK_LIBRARIES ${ARMPL_LIBRARIES})
endif()

try_compile(
  KOKKOSKERNELS_TRY_COMPILE_ARMPL ${KOKKOSKERNELS_TOP_BUILD_DIR}/tpl_tests
  ${KOKKOSKERNELS_TOP_SOURCE_DIR}/cmake/compile_tests/armpl.cpp
  LINK_LIBRARIES -l${ARMPL_LIB} -lgfortran -lamath -lm
  OUTPUT_VARIABLE KOKKOSKERNELS_TRY_COMPILE_ARMPL_OUT)
if(NOT KOKKOSKERNELS_TRY_COMPILE_ARMPL)
  message(FATAL_ERROR "KOKKOSKERNELS_TRY_COMPILE_ARMPL_OUT=${KOKKOSKERNELS_TRY_COMPILE_ARMPL_OUT}")
else()
  # KokkosKernels::ARMPL is an alias to the ARMPL target.
  # Let's add in the libgfortran and libm dependencies for users here.
  get_target_property(ARMPL_INTERFACE_LINK_LIBRARIES KokkosKernels::ARMPL INTERFACE_LINK_LIBRARIES)
  set(ARMPL_INTERFACE_LINK_LIBRARIES "${ARMPL_INTERFACE_LINK_LIBRARIES};-lgfortran;-lm")
  set_target_properties(ARMPL PROPERTIES INTERFACE_LINK_LIBRARIES "${ARMPL_INTERFACE_LINK_LIBRARIES}")
endif()
