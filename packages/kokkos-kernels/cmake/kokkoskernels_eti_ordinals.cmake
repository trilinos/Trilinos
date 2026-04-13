set(ORDINALS ORDINAL_INT ORDINAL_INT64_T)
set(ORDINAL_INT_CPP_TYPE int)
set(ORDINAL_INT64_T_CPP_TYPE int64_t)

set(KOKKOSKERNELS_INST_ORDINAL_INT_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
set(KOKKOSKERNELS_INST_ORDINAL_INT64_T_DEFAULT OFF)

kokkoskernels_add_option("INST_ORDINAL_INT" ${KOKKOSKERNELS_INST_ORDINAL_INT_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the ordinal type int.  This option is KokkosKernels_INST_ORDINAL_INT=ON by default. Default: ON")

kokkoskernels_add_option("INST_ORDINAL_INT64_T" ${KOKKOSKERNELS_INST_ORDINAL_INT64_T_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the ordinal type int64_t.  This option is KokkosKernels_INST_ORDINAL_INT64_T=OFF by default. Default: OFF")

if(KOKKOSKERNELS_INST_ORDINAL_INT)
  list(APPEND ORDINAL_LIST "int")
endif()

if(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
  list(APPEND ORDINAL_LIST "int64_t")
endif()
