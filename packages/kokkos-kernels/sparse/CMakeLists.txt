list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse/impl)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/sparse/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse/tpls)

#Add tpls source file
list(APPEND SOURCES sparse/tpls/KokkosKernels_tpl_handles.cpp)

##################
#                #
# ETI generation #
#                #
##################

macro(gen_sparse_eti FUNC_NAME SUBDIR)
  kokkoskernels_generate_eti(${FUNC_NAME} ${SUBDIR}
    COMPONENTS  sparse
    HEADER_LIST ETI_HEADERS
    SOURCE_LIST SOURCES
    TYPE_LISTS  FLOATS ORDINALS OFFSETS LAYOUTS DEVICES)
endmacro()

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros
gen_sparse_eti(Sparse_spmv_struct        spmv)
gen_sparse_eti(Sparse_spmv_mv_struct     spmv)
gen_sparse_eti(Sparse_spmv_bsrmatrix     spmv)
gen_sparse_eti(Sparse_spmv_mv_bsrmatrix  spmv)
gen_sparse_eti(Sparse_spmv               spmv)
gen_sparse_eti(Sparse_spmv_mv            spmv)
gen_sparse_eti(Sparse_spgemm_symbolic    spgemm_symbolic)
gen_sparse_eti(Sparse_spgemm_numeric     spgemm_numeric)
kokkoskernels_generate_eti(Sparse_spgemm_noreuse spgemm_noreuse
  COMPONENTS  sparse
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS ORDINALS OFFSETS DEVICES) # do not do LAYOUTS for spgemm_noreuse
gen_sparse_eti(Sparse_bspgemm_numeric    bspgemm_numeric)
gen_sparse_eti(Sparse_spgemm_jacobi      spgemm_jacobi)
# NOTE: SpAdd symbolic doesn't use scalars directly,
# but it needs the type to use handles.
gen_sparse_eti(Sparse_spadd_symbolic         spadd_symbolic)
gen_sparse_eti(Sparse_spadd_numeric          spadd_numeric)
gen_sparse_eti(Sparse_spiluk_symbolic        spiluk_symbolic)
gen_sparse_eti(Sparse_spiluk_numeric         spiluk_numeric)
gen_sparse_eti(Sparse_par_ilut_symbolic      par_ilut_symbolic)
gen_sparse_eti(Sparse_par_ilut_numeric       par_ilut_numeric)
gen_sparse_eti(Sparse_gmres                  gmres)
gen_sparse_eti(Sparse_sptrsv_symbolic        sptrsv_symbolic)
gen_sparse_eti(Sparse_sptrsv_solve           sptrsv_solve)
gen_sparse_eti(Sparse_trsv                   trsv)
gen_sparse_eti(Sparse_gauss_seidel_symbolic  gauss_seidel_symbolic)
gen_sparse_eti(Sparse_gauss_seidel_numeric   gauss_seidel_numeric)
gen_sparse_eti(Sparse_gauss_seidel_apply     gauss_seidel_apply)
