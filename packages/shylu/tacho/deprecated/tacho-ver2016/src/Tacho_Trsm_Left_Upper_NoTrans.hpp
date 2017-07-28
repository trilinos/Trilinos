#ifndef __TACHO_TRSM_LEFT_UPPER_NOTRANS_HPP__
#define __TACHO_TRSM_LEFT_UPPER_NOTRANS_HPP__

/// \file Tacho_Trsm_Left_Upper_NoTrans.hpp
/// \brief Triangular solve 
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///


// - Dense TRSM
#include "Tacho_Trsm_Left_Upper_NoTrans_ExternalBlas.hpp"
//#include "Tacho_Trsm_Left_Upper_NoTrans_DenseByBlocks.hpp"

// - Sparse TRSM
//#include "Tacho_Trsm_Left_Upper_NoTrans_SparseSparseUnblocked.hpp"
//#include "Tacho_Trsm_Left_Upper_NoTrans_SparseDenseUnblocked.hpp"

// - TRSM for supernodal algorithms
//#include "Tacho_Trsm_Left_Upper_NoTrans_SparseSparseSuperNodes.hpp"
//#include "Tacho_Trsm_Left_Upper_NoTrans_SparseSparseSuperNodesByBlocks.hpp"

#include "Tacho_Trsm_Left_Upper_NoTrans_SparseDenseSuperNodes.hpp"

#endif
