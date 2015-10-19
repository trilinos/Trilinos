#pragma once
#ifndef __CONTROL_HPP__
#define __CONTROL_HPP__

#include "util.hpp"

/// \file control.hpp
/// \brief A collection of control trees composing high-level variants of algorithms.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
using namespace std;

namespace Tacho {

  // forward declaration for control tree
  template<int ArgAlgo, int ArgVariant>
  struct Control {
    static constexpr int Self[2] = { ArgAlgo, ArgVariant };
  };

  // - Level based incomplete CholByBlocks Variant 1
  // * partitioned block matrix : sparse
  // * a block : sparse
  // * nested block structure : not applicable
  // * dependence is made using the right look algorithm
  template<> struct Control<AlgoChol::ByBlocks,Variant::One> {
    static constexpr int Chol[2] = { AlgoChol::UnblockedOpt, Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::ForFactorBlocked, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::ForFactorBlocked, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::ForFactorBlocked, Variant::One };
  };

  // - Complete CholByBlocks Variant 1
  // * partitioned block matrix : sparse
  // * a block : sparse
  // * dependence is made using the right look algorithm
  // * no nested block structure

  // * sparse block partitioned matrix
  // * a block has one-level nested sparse blocks
  // * a block is sparse
  // * nested block structure : 
  //   - a diagonal block may have nested sparse blocks 
  //   - nested blocks are sparse
  // * dependence is made using the right look algorithm
  template<> struct Control<AlgoChol::ByBlocks,Variant::Two> {
    static constexpr int Chol[2] = { AlgoChol::HierByBlocks, Variant::One }; // respawn algorithm

    // for testing here we use blocked version
    static constexpr int Trsm[2] = { AlgoTrsm::ForFactorBlocked, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::ForFactorBlocked, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::ForFactorBlocked, Variant::One };
  };

  // - Hierarchical CholByBlocks Variant 1
  // * partitioned block matrix : sparse
  // * a block : sparse
  // * dependence is made using the right look algorithm
  // * no nested block structure

  // * sparse block partitioned matrix
  // * a block has one-level nested sparse blocks
  // * a block is sparse
  // * nested block structure : 
  //   - a diagonal block may have nested sparse blocks 
  //   - nested blocks are sparse
  // * dependence is made using the right look algorithm
  template<> struct Control<AlgoChol::HierByBlocks,Variant::One> {
    static constexpr int CholScalar[2]   = { AlgoChol::UnblockedOpt, Variant::One };
    static constexpr int CholByBlocks[2] = { AlgoChol::ByBlocks, Variant::One }; 
  };

}

#endif
