/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
          const Teuchos::ArrayView<const LO>& perm) const
  {
    using Teuchos::ArrayRCP;
    const size_t numRows = X_out.getLocalLength ();
    const size_t numVecs = X_in.getNumVectors ();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const InScalar> X_in_j = X_in.getData(j);
      ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for (size_t i = 0; i < numRows; ++i) {
        const size_t i_perm = perm[i];
        X_out_j[i] = X_in_j[i_perm];
      }
    }
  }

  void
  scatter (MV_in& X_in,
           const MV_out& X_out,
           const Teuchos::ArrayView<const LO>& perm) const
  {
    using Teuchos::ArrayRCP;
    const size_t numRows = X_out.getLocalLength();
    const size_t numVecs = X_in.getNumVectors();
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<InScalar> X_in_j = X_in.getDataNonConst(j);
      ArrayRCP<const OutScalar> X_out_j = X_out.getData(j);
      for (size_t i = 0; i < numRows; ++i) {
        const size_t i_perm = perm[i];
        X_in_j[i_perm] = X_out_j[i];
      }
    }
  }

  /******************/
  /* View <==> View */
  /******************/
  template<typename InView, typename OutView>
  void gatherViewToView(OutView& X_out,
                        InView& X_in,
                        const Teuchos::ArrayView<const LO>& perm) const
  {
    //note: j is col, i is row
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < X_out.extent(0); ++i) {
        const LO i_perm = perm[i];
        X_out(i, j) = X_in(i_perm, j);
      }
    }
  }

  template<typename InView, typename OutView>
  void scatterViewToView(InView& X_in,
                         OutView& X_out,
                         const Teuchos::ArrayView<const LO>& perm) const
  {
    for(size_t j = 0; j < X_out.extent(1); ++j) {
      for(size_t i = 0; i < X_out.extent(0); ++i) {
        const LO i_perm = perm[i];
        X_in(i_perm, j) = X_out(i, j);
      }
    }
  }

  /*******************************/
  /* MV <==> View specialization */
  /*******************************/
  template<typename InView>
  void gatherMVtoView(MV_out& X_out,
                      InView& X_in,
                      const Teuchos::ArrayView<const LO>& perm) const
  {
    //note: j is col, i is row
    for(size_t j = 0; j < X_out.getNumVectors(); ++j) {
      Teuchos::ArrayRCP<OutScalar> X_out_j = X_out.getDataNonConst(j);
      for(size_t i = 0; i < X_out.getLocalLength(); ++i) {
        const LO i_perm = perm[i];
        X_out_j[i] = X_in(i_perm, j);
      }
    }
  }

  template<typename InView>
  void scatterMVtoView(InView& X_in,
                       MV_out& X_out,
                       const Teuchos::ArrayView<const LO>& perm) const
  {
    for(size_t j = 0; j < X_in.extent(1); ++j) {
      Teuchos::ArrayRCP<const OutScalar> X_out_j = X_out.getData(j);
      for(size_t i = 0; i < X_out.getLocalLength(); ++i) {
        const LO i_perm = perm[i];
        X_in(i_perm, j) = X_out_j[i];
      }
    }
  }
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_MULTIVECTORLOCALGATHERSCATTER_HPP
