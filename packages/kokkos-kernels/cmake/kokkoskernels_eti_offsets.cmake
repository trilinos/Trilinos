if(KOKKOSKERNELS_HAS_TRILINOS)
  # In a Trilinos build, size_t is the default offset because this is what Tpetra uses
  # TODO: update this when Tpetra can use different offsets
  set(KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
  set(KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT OFF)
else()
  # But in a standalone KokkosKernels build, int is the default offset type
  # This provides the maximum TPL compatibility
  set(KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT OFF)
  set(KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
endif()

set(OFFSETS OFFSET_INT OFFSET_SIZE_T)
set(OFFSET_INT_CPP_TYPE int)
set(OFFSET_SIZE_T_CPP_TYPE size_t)
#GLOBAL_SET(KokkosKernels_INST_OFFSET_INT64_T_DEFAULT  OFF)

kokkoskernels_add_option("INST_OFFSET_INT" ${KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the offset type int.  This option is KokkosKernels_INST_OFFSET_INT=OFF by default. Default: ${KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT}")

kokkoskernels_add_option("INST_OFFSET_SIZE_T" ${KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the offset type size_t.  This option is KokkosKernels_INST_OFFSET_SIZE_T=ON by default. Default: ${KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT}")

if(KOKKOSKERNELS_INST_OFFSET_INT)
  list(APPEND OFFSET_LIST "int")
endif()

if(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  list(APPEND OFFSET_LIST "size_t")
endif()
