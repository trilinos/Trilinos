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
    static constexpr int Chol[2] = { AlgoChol::Unblocked,        Variant::One };
    static constexpr int Trsm[2] = { AlgoTrsm::ForFactorization, Variant::One };
    static constexpr int Herk[2] = { AlgoHerk::ForFactorization, Variant::One };
    static constexpr int Gemm[2] = { AlgoGemm::ForFactorization, Variant::One };
  };

  // // - CholByBlocks Variant 2
  // // * diagonal blocks have nested dense blocks
  // template<> struct Control<AlgoChol::ByBlocks,Variant::Two> {
  //   static constexpr int Chol[2] = { AlgoChol::NestedDenseBlock, Variant::One }; 
  //   static constexpr int Trsm[2] = { AlgoTrsm::ForFactorBlocked, Variant::One };
  //   static constexpr int Herk[2] = { AlgoHerk::ForFactorBlocked, Variant::One };
  //   static constexpr int Gemm[2] = { AlgoGemm::ForFactorBlocked, Variant::One };
  // };

  // // - CholByBlocks Variant 3
  // // * all blocks have nested dense blocks (full supernodal algorithm)
  // // template<> struct Control<AlgoChol::ByBlocks,Variant::Three> {
  // //   static constexpr int Chol[2] = { AlgoChol::NestedDenseBlock, Variant::One }; 
  // //   static constexpr int Trsm[2] = { AlgoTrsm::NestedDenseBlock, Variant::One };
  // //   static constexpr int Herk[2] = { AlgoHerk::NestedDenseBlock, Variant::One };
  // //   static constexpr int Gemm[2] = { AlgoGemm::NestedDenseBlock, Variant::One };
  // // };

  // // - CholByBlocks Variant 4
  // // * diagonal blocks have nested hier dense blocks (hierarchical task scheduling)
  // template<> struct Control<AlgoChol::ByBlocks,Variant::Four> {
  //   static constexpr int Chol[2] = { AlgoChol::NestedDenseByBlocks, Variant::One }; 
  //   static constexpr int Trsm[2] = { AlgoTrsm::ForFactorBlocked,    Variant::One };
  //   static constexpr int Herk[2] = { AlgoHerk::ForFactorBlocked,    Variant::One };
  //   static constexpr int Gemm[2] = { AlgoGemm::ForFactorBlocked,    Variant::One };
  // };

  // // - CholByBlocks Variant 5
  // // * diagonal blocks have nested hier dense blocks (hierarchical task scheduling)
  // // template<> struct Control<AlgoChol::ByBlocks,Variant::Four> {
  // //   static constexpr int Chol[2] = { AlgoChol::NestedDenseByBlocks, Variant::One }; 
  // //   static constexpr int Trsm[2] = { AlgoTrsm::NestedDenseByBlocks, Variant::One };
  // //   static constexpr int Herk[2] = { AlgoHerk::NestedDenseByBlocks, Variant::One };
  // //   static constexpr int Gemm[2] = { AlgoGemm::NestedDenseByBlocks, Variant::One };
  // // };

  // // ----------------------------------------------------------------------------------

  // // - CholNestedDenseBlock
  // // * branch control between sparse and dense operations
  // template<> struct Control<AlgoChol::NestedDenseBlock,Variant::One> {
  //   static constexpr int CholSparse[2] = { AlgoChol::UnblockedOpt,   Variant::One };
  //   static constexpr int CholDense[2]  = { AlgoChol::ExternalLapack, Variant::One }; 
  // };

  // // - CholNestedDenseBlock
  // // * branch control between sparse and dense operations
  // template<> struct Control<AlgoChol::NestedDenseByBlocks,Variant::One> {
  //   static constexpr int CholSparse[2]        = { AlgoChol::UnblockedOpt,  Variant::One };
  //   static constexpr int CholDenseByBlocks[2] = { AlgoChol::DenseByBlocks, Variant::One }; 
  // };

  // // ----------------------------------------------------------------------------------

}

#endif
