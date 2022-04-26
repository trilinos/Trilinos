#ifndef TEST_SPARSE_HPP
#define TEST_SPARSE_HPP

#include "Test_Sparse_block_gauss_seidel.hpp"
#include "Test_Sparse_CrsMatrix.hpp"
#include "Test_Sparse_BlockCrsMatrix.hpp"
#include "Test_Sparse_BsrMatrix.hpp"
#include "Test_Sparse_findRelOffset.hpp"
#include "Test_Sparse_gauss_seidel.hpp"
#include "Test_Sparse_replaceSumInto.hpp"
#include "Test_Sparse_replaceSumIntoLonger.hpp"
#include "Test_Sparse_spadd.hpp"
#include "Test_Sparse_spgemm_jacobi.hpp"
#include "Test_Sparse_spgemm.hpp"
#include "Test_Sparse_spiluk.hpp"
#include "Test_Sparse_spmv.hpp"
//#include "Test_Sparse_spmv_blockcrs.hpp"
//#include "Test_Sparse_spmv_bsr.hpp"
#include "Test_Sparse_sptrsv.hpp"
#include "Test_Sparse_trsv.hpp"

// TPL specific tests, these require
// particular pairs of backend and TPL
// to actually define tests.

#include "Test_Sparse_Utils_cusparse.hpp"

#include "Test_Sparse_rocsparse.hpp"

#endif  // TEST_SPARSE_HPP
