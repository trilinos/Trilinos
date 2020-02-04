#!/bin/bash
KOKKOSKERNELS_PATH=$1
cd ${KOKKOSKERNELS_PATH}/src/impl
mkdir generated_specializations_hpp
mkdir generated_specializations_cpp

#sgs
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse_ml.bash gauss_seidel_symbolic KokkosSparse_gauss_seidel_symbolic KokkosSparse_gauss_seidel_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse_ml.bash gauss_seidel_numeric KokkosSparse_gauss_seidel_numeric KokkosSparse_gauss_seidel_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse_ml.bash gauss_seidel_apply KokkosSparse_gauss_seidel_apply KokkosSparse_gauss_seidel_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}


#spgemm_symbolic
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse_ml.bash spgemm_symbolic KokkosSparse_spgemm_symbolic KokkosSparse_spgemm_symbolic_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
#spgemm_numeric
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse_ml.bash spgemm_numeric  KokkosSparse_spgemm_numeric  KokkosSparse_spgemm_numeric_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}

#trsv
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash trsv KokkosSparse_trsv KokkosSparse_trsv_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}

#sptrsv
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash sptrsv_solve KokkosSparse_sptrsv_solve KokkosSparse_sptrsv_solve_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash sptrsv_symbolic KokkosSparse_sptrsv_symbolic KokkosSparse_sptrsv_symbolic_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}

#spmv
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spmv KokkosSparse_spmv KokkosSparse_spmv_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spmv KokkosSparse_spmv_mv KokkosSparse_spmv_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}


#spmv_struct
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spmv KokkosSparse_spmv_struct KokkosSparse_spmv_struct_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spmv KokkosSparse_spmv_mv_struct KokkosSparse_spmv_struct_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}

#spiluk
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spiluk_numeric KokkosSparse_spiluk_numeric KokkosSparse_spiluk_numeric_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function_sparse.bash spiluk_symbolic KokkosSparse_spiluk_symbolic KokkosSparse_spiluk_symbolic_spec.hpp KokkosSparse ${KOKKOSKERNELS_PATH}

#abs
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash abs KokkosBlas1_abs KokkosBlas1_abs_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash abs KokkosBlas1_abs_mv KokkosBlas1_abs_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#axpby
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash axpby KokkosBlas1_axpby KokkosBlas1_axpby_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash axpby KokkosBlas1_axpby_mv KokkosBlas1_axpby_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#dot
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash dot KokkosBlas1_dot KokkosBlas1_dot_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash dot KokkosBlas1_dot_mv KokkosBlas1_dot_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#iamax
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash iamax KokkosBlas1_iamax KokkosBlas1_iamax_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash iamax KokkosBlas1_iamax_mv KokkosBlas1_iamax_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#mult
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash mult KokkosBlas1_mult KokkosBlas1_mult_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash mult KokkosBlas1_mult_mv KokkosBlas1_mult_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#nrm1
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm1 KokkosBlas1_nrm1 KokkosBlas1_nrm1_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm1 KokkosBlas1_nrm1_mv KokkosBlas1_nrm1_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#nrm2
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm2 KokkosBlas1_nrm2 KokkosBlas1_nrm2_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm2 KokkosBlas1_nrm2_mv KokkosBlas1_nrm2_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#nrm2w
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm2w KokkosBlas1_nrm2w KokkosBlas1_nrm2w_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrm2w KokkosBlas1_nrm2w_mv KokkosBlas1_nrm2w_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#nrminf
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrminf KokkosBlas1_nrminf KokkosBlas1_nrminf_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash nrminf KokkosBlas1_nrminf_mv KokkosBlas1_nrminf_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#reciprocal
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash reciprocal KokkosBlas1_reciprocal KokkosBlas1_reciprocal_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash reciprocal KokkosBlas1_reciprocal_mv KokkosBlas1_reciprocal_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#scal
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash scal KokkosBlas1_scal KokkosBlas1_scal_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash scal KokkosBlas1_scal_mv KokkosBlas1_scal_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#sum
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash sum KokkosBlas1_sum KokkosBlas1_sum_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash sum KokkosBlas1_sum_mv KokkosBlas1_sum_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#update
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash update KokkosBlas1_update KokkosBlas1_update_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash update KokkosBlas1_update_mv KokkosBlas1_update_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#gemv
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash gemv KokkosBlas2_gemv KokkosBlas2_gemv_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#gemm
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash gemm KokkosBlas3_gemm KokkosBlas3_gemm_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}

#gesv
${KOKKOSKERNELS_PATH}/scripts/generate_specialization_function.bash gesv KokkosBlas_gesv KokkosBlas_gesv_spec.hpp KokkosBlas ${KOKKOSKERNELS_PATH}
