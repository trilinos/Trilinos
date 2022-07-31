#ifndef TEST_BLAS_HPP
#define TEST_BLAS_HPP

#include "Test_Blas_gesv.hpp"
#include "Test_Blas_trtri.hpp"

// Blas 1
#include "Test_Blas1_abs.hpp"
#include "Test_Blas1_asum.hpp"
#include "Test_Blas1_axpby.hpp"
#include "Test_Blas1_axpy.hpp"
#include "Test_Blas1_dot.hpp"
#include "Test_Blas1_iamax.hpp"
#include "Test_Blas1_mult.hpp"
#include "Test_Blas1_nrm1.hpp"
#include "Test_Blas1_nrm2_squared.hpp"
#include "Test_Blas1_nrm2.hpp"
#include "Test_Blas1_nrm2w_squared.hpp"
#include "Test_Blas1_nrm2w.hpp"
#include "Test_Blas1_nrminf.hpp"
#include "Test_Blas1_reciprocal.hpp"
#include "Test_Blas1_scal.hpp"
#include "Test_Blas1_sum.hpp"
#include "Test_Blas1_update.hpp"

// Team Blas 1
#include "Test_Blas1_team_abs.hpp"
#include "Test_Blas1_team_axpby.hpp"
#include "Test_Blas1_team_axpy.hpp"
#include "Test_Blas1_team_dot.hpp"
#include "Test_Blas1_team_mult.hpp"
#include "Test_Blas1_team_nrm2.hpp"
#include "Test_Blas1_team_scal.hpp"
#include "Test_Blas1_team_update.hpp"

// Blas 2
#include "Test_Blas2_gemv.hpp"

// Team Blas 2
#include "Test_Blas2_team_gemv.hpp"

// Blas 3
#include "Test_Blas3_gemm.hpp"
#include "Test_Blas3_trmm.hpp"
#include "Test_Blas3_trsm.hpp"

// TPLs
#include "Test_Blas_rocblas.hpp"

#endif  // TEST_BLAS_HPP
