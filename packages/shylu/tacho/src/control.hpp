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
  template<int ArgAlgo>
  struct Control {
    static const int Self = ArgAlgo;
  };

  // -- level based incomplete CholByBlocks
  // * sparse block partitioned matrix
  // * a block is sparse
  // * dependence is made using the right look algorithm
  template<> struct Control<AlgoChol::ByBlocksVar1> {
    static const int Chol = AlgoChol::UnblockedOpt1;
    static const int Trsm = AlgoTrsm::ForFactorBlocked;
    static const int Herk = AlgoHerk::ForFactorBlocked;
    static const int Gemm = AlgoGemm::ForFactorBlocked;
  };

  // // -- complte CholByBlocks
  // // * sparse block partitioned matrix
  // // * a block has one-level nested sparse blocks
  // // * dependence is made using the right look algorithm
  // template<> struct Control<AlgoChol::ByBlocksVar2> {
    
  // };

}

#endif
