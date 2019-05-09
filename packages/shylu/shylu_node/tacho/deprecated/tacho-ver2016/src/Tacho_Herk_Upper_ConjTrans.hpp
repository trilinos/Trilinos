#ifndef __TACHO_HERK_UPPER_CONJTRANS_HPP__
#define __TACHO_HERK_UPPER_CONJTRANS_HPP__

/// \file Tacho_Herk_Upper_ConjTrans.hpp
/// \brief Hermitian rank-k update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// - Dense HERK
#include "Tacho_Herk_Upper_ConjTrans_ExternalBlas.hpp"
#include "Tacho_Herk_Upper_ConjTrans_DenseByBlocks.hpp"

// - Sparse HERK
#include "Tacho_Herk_Upper_ConjTrans_SparseSparseUnblocked.hpp"

// - HERK for supernodal algorithms
#include "Tacho_Herk_Upper_ConjTrans_SparseSparseSuperNodes.hpp"
#include "Tacho_Herk_Upper_ConjTrans_SparseSparseSuperNodesByBlocks.hpp"

#endif
