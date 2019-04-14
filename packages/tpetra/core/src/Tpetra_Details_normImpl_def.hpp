// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_NORMIMPL_DEF_HPP
#define TPETRA_DETAILS_NORMIMPL_DEF_HPP

/// \file Tpetra_Details_normImpl_def.hpp
/// \brief Definition of the Tpetra::Details::normImpl function

#include "Tpetra_Details_normImpl_decl.hpp" // for the enum
// #include "Tpetra_Details_Behavior.hpp"
// #include "Tpetra_Details_isInterComm.hpp"
// #include "Tpetra_Details_Profiling.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "KokkosBlas.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

template<class RV, class XMV>
void
lclNormImpl (const RV& normsOut,
             const XMV& X,
             const size_t numVecs,
             const Teuchos::ArrayView<const size_t>& whichVecs,
             const bool constantStride,
             const EWhichNorm whichNorm)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  using mag_type = typename RV::non_const_value_type;

  static_assert (static_cast<int> (RV::Rank) == 1,
                 "Tpetra::MultiVector::lclNormImpl: "
                 "The first argument normsOut must have rank 1.");
  static_assert (Kokkos::Impl::is_view<XMV>::value,
                 "Tpetra::MultiVector::lclNormImpl: "
                 "The second argument X is not a Kokkos::View.");
  static_assert (static_cast<int> (XMV::Rank) == 2,
                 "Tpetra::MultiVector::lclNormImpl: "
                 "The second argument X must have rank 2.");

  const size_t lclNumRows = static_cast<size_t> (X.extent (0));
  TEUCHOS_TEST_FOR_EXCEPTION
    (lclNumRows != 0 && constantStride &&
     static_cast<size_t> (X.extent (1)) != numVecs,
     std::logic_error, "Constant Stride X's dimensions are " << X.extent (0)
     << " x " << X.extent (1) << ", which differ from the local dimensions "
     << lclNumRows << " x " << numVecs << ".  Please report this bug to "
     "the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (lclNumRows != 0 && ! constantStride &&
     static_cast<size_t> (X.extent (1)) < numVecs,
     std::logic_error, "Strided X's dimensions are " << X.extent (0) << " x "
     << X.extent (1) << ", which are incompatible with the local dimensions "
     << lclNumRows << " x " << numVecs << ".  Please report this bug to "
     "the Tpetra developers.");

  if (lclNumRows == 0) {
    const mag_type zeroMag = Kokkos::ArithTraits<mag_type>::zero ();
    Kokkos::deep_copy (normsOut, zeroMag);
  }
  else { // lclNumRows != 0
    if (constantStride) {
      if (whichNorm == NORM_INF) {
        KokkosBlas::nrminf (normsOut, X);
      }
      else if (whichNorm == NORM_ONE) {
        KokkosBlas::nrm1 (normsOut, X);
      }
      else if (whichNorm == NORM_TWO) {
        KokkosBlas::nrm2_squared (normsOut, X);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "Should never get here!");
      }
    }
    else { // not constant stride
      // NOTE (mfh 15 Jul 2014, 11 Apr 2019) This does a kernel launch
      // for every column.  It might be better to have a kernel that
      // does the work all at once.  On the other hand, we don't
      // prioritize performance of MultiVector views of noncontiguous
      // columns.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t X_col = constantStride ? k : whichVecs[k];
        if (whichNorm == NORM_INF) {
          KokkosBlas::nrminf (subview (normsOut, k),
                              subview (X, ALL (), X_col));
        }
        else if (whichNorm == NORM_ONE) {
          KokkosBlas::nrm1 (subview (normsOut, k),
                            subview (X, ALL (), X_col));
        }
        else if (whichNorm == NORM_TWO) {
          KokkosBlas::nrm2_squared (subview (normsOut, k),
                                    subview (X, ALL (), X_col));
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, "Should never get here!");
        }
      } // for each column
    } // constantStride
  } // lclNumRows != 0
}

// Kokkos::parallel_for functor for applying square root to each
// entry of a 1-D Kokkos::View.
template<class ViewType>
class SquareRootFunctor {
public:
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;

  SquareRootFunctor (const ViewType& theView) :
    theView_ (theView)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type& i) const
  {
    typedef typename ViewType::non_const_value_type value_type;
    typedef Kokkos::Details::ArithTraits<value_type> KAT;
    theView_(i) = KAT::sqrt (theView_(i));
  }

private:
  ViewType theView_;
};

template<class RV>
void
gblNormImpl (const RV& normsOut,
             const Teuchos::Comm<int>* const comm,
             const bool distributed,
             const EWhichNorm whichNorm)
{
  using Teuchos::REDUCE_MAX;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef typename RV::non_const_value_type mag_type;

  const size_t numVecs = normsOut.extent (0);

  // If the MultiVector is distributed over multiple processes, do
  // the distributed (interprocess) part of the norm.  We assume
  // that the MPI implementation can read from and write to device
  // memory.
  //
  // replaceMap() may have removed some processes.  Those processes
  // have a null Map.  They must not participate in any collective
  // operations.  We ask first whether the Map is null, because
  // isDistributed() defers that question to the Map.  We still
  // compute and return local norms for processes not participating
  // in collective operations; those probably don't make any sense,
  // but it doesn't hurt to do them, since it's illegal to call
  // norm*() on those processes anyway.
  if (distributed && comm != nullptr) {
    // The calling process only participates in the collective if
    // both the Map and its Comm on that process are nonnull.
    //
    // MPI doesn't allow aliasing of arguments, so we have to make
    // a copy of the local sum.
    RV lclNorms ("MV::normImpl lcl", numVecs);
    Kokkos::deep_copy (lclNorms, normsOut);
    const mag_type* const lclSum = lclNorms.data ();
    mag_type* const gblSum = normsOut.data ();
    const int nv = static_cast<int> (numVecs);
    if (whichNorm == NORM_INF) {
      reduceAll<int, mag_type> (*comm, REDUCE_MAX, nv, lclSum, gblSum);
    } else {
      reduceAll<int, mag_type> (*comm, REDUCE_SUM, nv, lclSum, gblSum);
    }
  }

  if (whichNorm == NORM_TWO) {
    // Replace the norm-squared results with their square roots in
    // place, to get the final output.  If the device memory and
    // the host memory are the same, it probably doesn't pay to
    // launch a parallel kernel for that, since there isn't enough
    // parallelism for the typical MultiVector case.
    const bool inHostMemory =
      Kokkos::Impl::is_same<typename RV::memory_space,
      typename RV::host_mirror_space::memory_space>::value;
    if (inHostMemory) {
      for (size_t j = 0; j < numVecs; ++j) {
        normsOut(j) = Kokkos::Details::ArithTraits<mag_type>::sqrt (normsOut(j));
      }
    }
    else {
      // There's not as much parallelism now, but that's OK.  The
      // point of doing parallel dispatch here is to keep the norm
      // results on the device, thus avoiding a copy to the host
      // and back again.
      SquareRootFunctor<RV> f (normsOut);
      typedef typename RV::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      Kokkos::parallel_for (range_type (0, numVecs), f);
    }
  }
}

} // namespace Impl

template <class ValueType,
          class ArrayLayout,
          class DeviceType,
          class MagnitudeType>
void
normImpl (MagnitudeType norms[],
          const Kokkos::View<const ValueType**, ArrayLayout, DeviceType>& X,
          const EWhichNorm whichNorm,
          const Teuchos::ArrayView<const size_t>& whichVecs,
          const bool isConstantStride,
          const bool isDistributed,
          const Teuchos::Comm<int>* comm)
{
  using RV = Kokkos::View<MagnitudeType*, Kokkos::HostSpace>;
  //using XMV = Kokkos::View<const ValueType**, ArrayLayout, DeviceType>;
  //using pair_type = std::pair<size_t, size_t>;

  const size_t numVecs = isConstantStride ?
    static_cast<size_t> (X.extent (1)) :
    static_cast<size_t> (whichVecs.size ());
  if (numVecs == 0) {
    return; // nothing to do
  }
  RV normsOut (norms, numVecs);

  Impl::lclNormImpl (normsOut, X, numVecs, whichVecs,
                     isConstantStride, whichNorm);
  Impl::gblNormImpl (normsOut, comm, isDistributed, whichNorm);
}

} // namespace Details
} // namespace Tpetra

#define TPETRA_DETAILS_NORMIMPL_INSTANT( SC, NT ) \
namespace Details { \
  template void \
  normImpl< \
    Kokkos::ArithTraits<SC>::val_type, \
    Kokkos::LayoutLeft, \
    NT::device_type, \
    Kokkos::ArithTraits<SC>::mag_type \
  > (Kokkos::ArithTraits<SC>::mag_type[], \
     const Kokkos::View< \
       const Kokkos::ArithTraits<SC>::val_type**, \
       Kokkos::LayoutLeft, \
       NT::device_type>&, \
     const EWhichNorm, \
     const Teuchos::ArrayView<const size_t>&, \
     const bool, \
     const bool, \
     const Teuchos::Comm<int>* ); \
}

#endif // TPETRA_DETAILS_NORMIMPL_DEF_HPP
