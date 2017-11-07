#ifndef __TACHO_TRSM_LEFT_UPPER_CONJTRANS_HPP__
#define __TACHO_TRSM_LEFT_UPPER_CONJTRANS_HPP__

/// \file Tacho_Trsm_Left_Upper_ConjTrans.hpp
/// \brief Triangular solve 
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

//#include "trsm_l_u_ct_for_tri_solve_blocked.hpp"

// - Dense TRSM
#include "Tacho_Trsm_Left_Upper_ConjTrans_ExternalBlas.hpp"
#include "Tacho_Trsm_Left_Upper_ConjTrans_DenseByBlocks.hpp"

// - Sparse TRSM
#include "Tacho_Trsm_Left_Upper_ConjTrans_SparseSparseUnblocked.hpp"
//#include "Tacho_Trsm_Left_Upper_ConjTrans_SparseDenseUnblocked.hpp"

// - TRSM for supernodal algorithms
#include "Tacho_Trsm_Left_Upper_ConjTrans_SparseSparseSuperNodes.hpp"
#include "Tacho_Trsm_Left_Upper_ConjTrans_SparseSparseSuperNodesByBlocks.hpp"

#include "Tacho_Trsm_Left_Upper_ConjTrans_SparseDenseSuperNodes.hpp"

#endif
