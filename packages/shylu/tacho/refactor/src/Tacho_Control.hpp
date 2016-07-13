#ifndef __TACHO_CONTROL_HPP__
#define __TACHO_CONTROL_HPP__

/// \file Tacho_Control.hpp
/// \brief A collection of control trees composing high-level variants of algorithms.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  // forward declaration for control tree
  template<int ArgAlgo, int ArgVariant>
  struct Control {
    static constexpr int Self[2] = { ArgAlgo, ArgVariant };
  };

  // ----------------------------------------------------------------------------------

  // - DenseByBlocks
  template<> struct Control<AlgoGemm::DenseByBlocks,Variant::One> {
    static constexpr int Gemm[2] = { AlgoGemm::ExternalBlas, Variant::One };
  };

  template<> struct Control<AlgoGemm::DenseByBlocks,Variant::Two> {
    static constexpr int Gemm[2] = { AlgoGemm::InternalBlas, Variant::One };
  };

  template<> struct Control<AlgoHerk::DenseByBlocks,Variant::One> {
    static constexpr int Herk[2] = { AlgoHerk::ExternalBlas, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::ExternalBlas, Variant::One };
  };

  template<> struct Control<AlgoHerk::DenseByBlocks,Variant::Two> {
    static constexpr int Herk[2] = { AlgoHerk::InternalBlas, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::InternalBlas, Variant::One };
  };

  template<> struct Control<AlgoTrsm::DenseByBlocks,Variant::One> {
    static constexpr int Gemm[2] = { AlgoGemm::ExternalBlas, Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::ExternalBlas, Variant::One };
  };

  template<> struct Control<AlgoTrsm::DenseByBlocks,Variant::Two> {
    static constexpr int Gemm[2] = { AlgoGemm::InternalBlas, Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::InternalBlas, Variant::One };
  };

  template<> struct Control<AlgoChol::DenseByBlocks,Variant::One> {
    static constexpr int Chol[2] = { AlgoChol::ExternalLapack, Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::ExternalBlas,   Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::ExternalBlas,   Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::ExternalBlas,   Variant::One };
  };

  template<> struct Control<AlgoChol::DenseByBlocks,Variant::Two> {
    static constexpr int Chol[2] = { AlgoChol::ExternalLapack, Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::InternalBlas,   Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::InternalBlas,   Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::InternalBlas,   Variant::One };
  };

  // - SparseByBlocks
  template<> struct Control<AlgoChol::ByBlocks,Variant::One> {
    static constexpr int Chol[2] = { AlgoChol::Unblocked,             Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::SparseSparseUnblocked, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::SparseSparseUnblocked, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseSparseUnblocked, Variant::One };
  };

  // - SuperNodalByblocks
  template<> struct Control<AlgoChol::ByBlocks,Variant::Two> {
    static constexpr int Chol[2] = { AlgoChol::SuperNodes,    Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::SparseSparseSuperNodes, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::SparseSparseSuperNodes, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseSparseSuperNodes, Variant::One };
  };

  // - Fine grained SuperNodalByblocks
  template<> struct Control<AlgoChol::ByBlocks,Variant::Three> {
    static constexpr int Chol[2] = { AlgoChol::SuperNodesByBlocks,             Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::SparseSparseSuperNodesByBlocks, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::SparseSparseSuperNodesByBlocks, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseSparseSuperNodesByBlocks, Variant::One };
  };

  // - SparseByBlocks
  template<> struct Control<AlgoTriSolve::ByBlocks,Variant::One> {
    static constexpr int Trsm[2] = { AlgoTrsm::SparseDenseUnblocked, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseDenseUnblocked, Variant::One };
  };

  // - SuperNodalByblocks
  template<> struct Control<AlgoTriSolve::ByBlocks,Variant::Two> {
    static constexpr int Trsm[2] = { AlgoTrsm::SparseDenseSuperNodes, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseDenseSuperNodes, Variant::One };
  };

  // - Fine grained SuperNodalByblocks
  template<> struct Control<AlgoTriSolve::ByBlocks,Variant::Three> {
    static constexpr int Trsm[2] = { AlgoTrsm::SparseDenseSuperNodesByBlocks, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseDenseSuperNodesByBlocks, Variant::One };
  };


  // - ByBlocksSerial (test only)
  template<> struct Control<AlgoTriSolve::ByBlocksSerial,Variant::One> {
    static constexpr int Trsm[2] = { AlgoTrsm::SparseDenseSuperNodes, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::SparseDenseSuperNodes, Variant::One };
  };

  // // ----------------------------------------------------------------------------------

}

#endif
