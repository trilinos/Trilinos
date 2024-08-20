// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_MULTIVECTORLOCALGATHERSCATTER_HPP
#define IFPACK2_DETAILS_MULTIVECTORLOCALGATHERSCATTER_HPP

/// \file Ifpack2_Details_MultiVectorLocalGatherScatter.hpp
/// \brief Declaration and definition of the
///   Ifpack2::Details::MultiVectorLocalGatherScatter class.

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"

namespace Ifpack2 {
namespace Details {

/// \class MultiVectorLocalGatherScatter
/// \brief Implementation detail of Ifpack2::Container subclasses.
///
/// \warning This is an implementation detail of subclasses of
///   Container.  This class may cease to exist or change its
///   interface at any time.
///
/// DenseContainer and SparseContainer use this class to copy data
/// between the input ordering of the global operator's domain and
/// range (of apply() and weightedApply()), and the ordering of the
/// local operator.  The latter ordering is a permuted subset of the
/// former.  Hence, we've chosen the terms "gather" and "scatter" to
/// refer to the copy operation to resp. from the local operator's
/// ordering.
///
/// <tt>MV_in</tt> and <tt>MV_out</tt> are possibly different
/// specializations of Tpetra::MultiVector.  <tt>MV_in</tt>
/// corresponds to the input and output arguments of apply() and
/// weightedApply() in Container, and <tt>MV_out</tt> to the input and
/// output arguments of the local operator.  The two specializations
/// of Tpetra::MultiVector may have entirely different template
/// parameters, even different <tt>Scalar</tt>, <tt>LocalOrdinal</tt>,
/// or <tt>GlobalOrdinal</tt> types.  This is a good way to experiment
/// with mixed-precision computation, for example.  Since
/// <tt>MV_in</tt> and <tt>MV_out</tt> may be different types, it
/// makes sense to implement "local gather / scatter" as a separate
/// class that uses the public interface of Tpetra::MultiVector,
/// rather than an instance method (which would have to be templated).
template<class MV_in, class MV_out>
class MultiVectorLocalGatherScatter {
public:
  typedef typename MV_in::scalar_type InScalar;
  typedef typename MV_out::scalar_type OutScalar;
  typedef typename MV_in::local_ordinal_type LO;
  typedef typename MV_in::global_ordinal_type GO;
  typedef typename MV_in::node_type NO;

  /**************/
  /* MV <==> MV */
  /**************/
  void
  gather (MV_out& X_out,
          const MV_in& X_in,
          const Teuchos::ArrayView<const LO> perm) const
  {
    using Teuchos::ArrayRCP;
    const size_t numRows = X_out.getLocalLength ();
    const size_t numVecs = X_in.getNumVectors ();
    Kokkos::fence();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const InScalar> X_in_j = X_in.getData(j);
      ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for (size_t i = 0; i < numRows; ++i) {
        const size_t i_perm = perm[i];
        X_out_j[i] = X_in_j[i_perm];
      }
    }
    X_out.modify_host();
  }

  //Gather blocks (contiguous groups of blockSize rows)
  //X_out and X_in are point indexed, but perm uses block indices.
  //So X_out.getLocalLength() / blockSize gives the number of blocks.
  void
  gatherBlock (
        MV_out& X_out,
        const MV_in& X_in,
        const Teuchos::ArrayView<const LO> perm,
        LO blockSize) const
  {
    using Teuchos::ArrayRCP;
    const size_t numBlocks = X_out.getLocalLength() / blockSize;
    const size_t numVecs = X_in.getNumVectors ();
    Kokkos::fence();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const InScalar> X_in_j = X_in.getData(j);
      ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for (size_t i = 0; i < numBlocks; ++i) {
        const size_t i_perm = perm[i];
        for (LO k = 0; k < blockSize; k++) {
          X_out_j[i * blockSize + k] = X_in_j[i_perm * blockSize + k];
        }
      }
    }
    X_out.modify_host();
  }

  void
  scatter (MV_in& X_in,
           const MV_out& X_out,
           const Teuchos::ArrayView<const LO> perm) const
  {
    using Teuchos::ArrayRCP;
    const size_t numRows = X_out.getLocalLength();
    const size_t numVecs = X_in.getNumVectors();
    Kokkos::fence();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<InScalar> X_in_j = X_in.getDataNonConst(j);
      ArrayRCP<const OutScalar> X_out_j = X_out.getData(j);
      for (size_t i = 0; i < numRows; ++i) {
        const size_t i_perm = perm[i];
        X_in_j[i_perm] = X_out_j[i];
      }
    }
    X_out.modify_host();
  }

  void
  scatterBlock (
        MV_in& X_in,
        const MV_out& X_out,
        const Teuchos::ArrayView<const LO> perm,
        LO blockSize) const
  {
    using Teuchos::ArrayRCP;
    const size_t numBlocks = X_out.getLocalLength() / blockSize;
    const size_t numVecs = X_in.getNumVectors ();
    Kokkos::fence();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const InScalar> X_in_j = X_in.getData(j);
      ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for (size_t i = 0; i < numBlocks; ++i) {
        const size_t i_perm = perm[i];
        for (LO k = 0; k < blockSize; k++) {
           X_in_j[i_perm * blockSize + k] = X_out_j[i * blockSize + k];
        }
      }
    }
    X_out.modify_host();
  }

  /******************/
  /* View <==> View */
  /******************/
  template<typename InView, typename OutView>
  void gatherViewToView(OutView X_out,
                        const InView X_in,
                        const Teuchos::ArrayView<const LO> perm) const
  {
    //note: j is col, i is row
    Kokkos::fence(); // demonstrated via unit test failure
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < X_out.extent(0); ++i) {
        const LO i_perm = perm[i];
        X_out(i, j) = X_in(i_perm, j);
      }
    }
  }

  template<typename InView, typename OutView>
  void scatterViewToView(InView X_in,
                         const OutView X_out,
                         const Teuchos::ArrayView<const LO> perm) const
  {
    Kokkos::fence();
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < X_out.extent(0); ++i) {
        const LO i_perm = perm[i];
        X_in(i_perm, j) = X_out(i, j);
      }
    }
  }

  template<typename InView, typename OutView>
  void gatherViewToViewBlock(OutView X_out,
                             const InView X_in,
                             const Teuchos::ArrayView<const LO> perm,
                             LO blockSize) const
  {
    //note: j is col, i is row
    Kokkos::fence();
    size_t numBlocks = X_out.extent(0) / blockSize;
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < numBlocks; ++i) {
        const LO i_perm = perm[i];
        for(LO k = 0; k < blockSize; k++) {
          X_out(i * blockSize + k, j) = X_in(i_perm * blockSize + k, j);
        }
      }
    }
  }

  template<typename InView, typename OutView>
  void scatterViewToViewBlock(InView X_in,
                              const OutView X_out,
                              const Teuchos::ArrayView<const LO> perm,
                              LO blockSize) const
  {
    //note: j is col, i is row
    Kokkos::fence();
    size_t numBlocks = X_out.extent(0) / blockSize;
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < numBlocks; ++i) {
        const LO i_perm = perm[i];
        for(LO k = 0; k < blockSize; k++) {
          X_in(i_perm * blockSize + k, j) = X_out(i * blockSize + k, j);
        }
      }
    }
  }

  /*******************************/
  /* MV <==> View specialization */
  /*******************************/
  template<typename InView>
  void gatherMVtoView(MV_out X_out,
                      InView X_in,
                      const Teuchos::ArrayView<const LO> perm) const
  {
    //note: j is col, i is row
    Kokkos::fence();
    size_t numRows = X_out.getLocalLength();
    for(size_t j = 0; j < X_out.getNumVectors(); ++j) {
      Teuchos::ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for(size_t i = 0; i < numRows; ++i) {
        const LO i_perm = perm[i];
        X_out_j[i] = X_in(i_perm, j);
      }
    }
  }

  template<typename InView>
  void scatterMVtoView(InView X_in,
                       MV_out X_out,
                       const Teuchos::ArrayView<const LO> perm) const
  {
    size_t numRows = X_out.getLocalLength(); 
    Kokkos::fence();
    for(size_t j = 0; j < X_in.extent(1); ++j) {
      Teuchos::ArrayRCP<const OutScalar> X_out_j = X_out.getData(j);
      for(size_t i = 0; i < numRows; ++i) {
        const LO i_perm = perm[i];
        X_in(i_perm, j) = X_out_j[i];
      }
    }
  }

  template<typename InView>
  void gatherMVtoViewBlock(MV_out X_out,
                           InView X_in,
                           const Teuchos::ArrayView<const LO> perm,
                           LO blockSize) const
  {
    //note: j is col, i is row
    size_t numBlocks = X_out.getLocalLength() / blockSize;
    Kokkos::fence();
    for(size_t j = 0; j < X_out.getNumVectors(); ++j) {
      Teuchos::ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for(size_t i = 0; i < numBlocks; ++i) {
        const LO i_perm = perm[i];
        for(LO k = 0; k < blockSize; k++) {
          X_out_j[i * blockSize + k] = X_in(i_perm * blockSize + k, j);
        }
      }
    }
  }

  template<typename InView>
  void scatterMVtoViewBlock(InView X_in,
                            MV_out X_out,
                            const Teuchos::ArrayView<const LO> perm,
                            LO blockSize) const
  {
    size_t numBlocks = X_out.getLocalLength() / blockSize;
    Kokkos::fence();
    for(size_t j = 0; j < X_in.extent(1); ++j) {
      Teuchos::ArrayRCP<const OutScalar> X_out_j = X_out.getData(j);
      for(size_t i = 0; i < numBlocks; ++i) {
        const LO i_perm = perm[i];
        for(LO k = 0; k < blockSize; k++) {
          X_in(i_perm * blockSize + k, j) = X_out_j[i * blockSize + k];
        }
      }
    }
  }
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_MULTIVECTORLOCALGATHERSCATTER_HPP
