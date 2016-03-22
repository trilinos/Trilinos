#ifndef __TACHO_CHOL_UPPER_HPP__
#define __TACHO_CHOL_UPPER_HPP__


/// \file Tacho_Chol_Upper.hpp
/// \brief Upper Cholesky factorization variations
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// testing task-data parallelism
//#include "chol_u_unblocked_dummy.hpp"

// triple for loop
//#include "Tacho_Chol_Upper_Unblocked_var1.hpp"
//#include "Tacho_Chol_Upper_Unblocked_var2.hpp"

// tools for supernodal algorithms
#include "Tacho_Chol_Upper_ExternalLapack.hpp"
#include "Tacho_Chol_Upper_DenseByBlocks.hpp"

// partitioned block algorithms: see control.hpp
//#include "chol_u_right_look_by_blocks.hpp"
//#include "chol_u_nested_dense_block.hpp"
//#include "chol_u_nested_dense_by_blocks.hpp"


#endif
