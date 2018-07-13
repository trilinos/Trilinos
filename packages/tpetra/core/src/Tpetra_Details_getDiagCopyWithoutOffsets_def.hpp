/*
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
*/

#ifndef TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DEF_HPP
#define TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DEF_HPP

/// \file Tpetra_Details_getDiagCopyWithoutOffsets_def.hpp
/// \brief Definition of
///   Tpetra::Details::getDiagCopyWithoutOffsetsNotFillComplete (an
///   implementation detail of Tpetra::CrsMatrix).
///
/// This function, and any declarations and/or definitions in it, are
/// implementation details of Tpetra::CrsMatrix.

#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {
namespace Details {

// Work-around for #499: Implementation of one-argument (no offsets)
// getLocalDiagCopy for the NOT fill-complete case.
//
// NOTE (mfh 18 Jul 2016) This calls functions that are NOT GPU device
// functions!  Thus, we do NOT use KOKKOS_INLINE_FUNCTION or
// KOKKOS_FUNCTION here, because those attempt to mark the functions
// they modify as CUDA device functions.  This functor is ONLY for
// non-CUDA execution spaces!
template<class SC, class LO, class GO, class NT>
class GetLocalDiagCopyWithoutOffsetsNotFillCompleteFunctor {
public:
  typedef ::Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef ::Tpetra::Vector<SC, LO, GO, NT> vec_type;

  typedef typename vec_type::impl_scalar_type IST;
  // The output Vector determines the execution space.
  typedef typename vec_type::device_type device_type;

private:
  typedef typename vec_type::dual_view_type::host_mirror_space::execution_space host_execution_space;
  typedef typename vec_type::map_type map_type;

  static bool
  graphIsSorted (const row_matrix_type& A)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
    typedef Tpetra::RowGraph<LO, GO, NT> row_graph_type;

    // We conservatively assume not sorted.  RowGraph lacks an
    // "isSorted" predicate, so we can't know for sure unless the cast
    // to CrsGraph succeeds.
    bool sorted = false;

    RCP<const row_graph_type> G_row = A.getGraph ();
    if (! G_row.is_null ()) {
      RCP<const crs_graph_type> G_crs =
        rcp_dynamic_cast<const crs_graph_type> (G_row);
      if (! G_crs.is_null ()) {
        sorted = G_crs->isSorted ();
      }
    }

    return sorted;
  }

public:
  // lclNumErrs [out] Total count of errors on this process.
  GetLocalDiagCopyWithoutOffsetsNotFillCompleteFunctor (LO& lclNumErrs,
                                                        vec_type& diag,
                                                        const row_matrix_type& A) :
    A_ (A),
    lclRowMap_ (A.getRowMap ()->getLocalMap ()),
    lclColMap_ (A.getColMap ()->getLocalMap ()),
    sorted_ (graphIsSorted (A))
  {
    const LO lclNumRows = static_cast<LO> (diag.getLocalLength ());
    {
      const LO matLclNumRows =
        static_cast<LO> (lclRowMap_.getNodeNumElements ());
      TEUCHOS_TEST_FOR_EXCEPTION
        (lclNumRows != matLclNumRows, std::invalid_argument,
         "diag.getLocalLength() = " << lclNumRows << " != "
         "A.getRowMap()->getNodeNumElements() = " << matLclNumRows << ".");
    }

    // Side effects start below this point.

    diag.template modify<Kokkos::HostSpace> ();
    D_lcl_ = diag.template getLocalView<Kokkos::HostSpace> ();
    D_lcl_1d_ = Kokkos::subview (D_lcl_, Kokkos::ALL (), 0);

    Kokkos::RangePolicy<host_execution_space, LO> range (0, lclNumRows);
    lclNumErrs = 0;
    Kokkos::parallel_reduce (range, *this, lclNumErrs);

    // sync changes back to device, since the user doesn't know that
    // we had to run on host.
    diag.template sync<typename device_type::memory_space> ();
  }

  void operator () (const LO& lclRowInd, LO& errCount) const {
    using KokkosSparse::findRelOffset;

    D_lcl_1d_(lclRowInd) = Kokkos::Details::ArithTraits<IST>::zero ();
    const GO gblInd = lclRowMap_.getGlobalElement (lclRowInd);
    const LO lclColInd = lclColMap_.getLocalElement (gblInd);

    if (lclColInd == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
      errCount++;
    }
    else { // row index is also in the column Map on this process
      LO numEnt;
      const LO* lclColInds;
      const SC* curVals;
      const LO err = A_.getLocalRowViewRaw (lclRowInd, numEnt, lclColInds, curVals);
      if (err != 0) {
        errCount++;
      }
      else {
        // The search hint is always zero, since we only call this
        // once per row of the matrix.
        const LO hint = 0;
        const LO offset =
          findRelOffset (lclColInds, numEnt, lclColInd, hint, sorted_);
        if (offset == numEnt) { // didn't find the diagonal column index
          errCount++;
        }
        else {
          D_lcl_1d_(lclRowInd) = curVals[offset];
        }
      }
    }
  }

private:
  const row_matrix_type& A_;
  typename map_type::local_map_type lclRowMap_;
  typename map_type::local_map_type lclColMap_;
  typename vec_type::dual_view_type::t_host D_lcl_;
  decltype (Kokkos::subview (D_lcl_, Kokkos::ALL (), 0)) D_lcl_1d_;
  const bool sorted_;
};


template<class SC, class LO, class GO, class NT>
LO
getLocalDiagCopyWithoutOffsetsNotFillComplete ( ::Tpetra::Vector<SC, LO, GO, NT>& diag,
                                                const ::Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                                                const bool debug)
{
  using ::Tpetra::Details::gathervPrint;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  typedef GetLocalDiagCopyWithoutOffsetsNotFillCompleteFunctor<SC,
    LO, GO, NT> functor_type;

  // The functor's constructor does error checking and executes the
  // thread-parallel kernel.

  LO lclNumErrs = 0;

  if (debug) {
    int lclSuccess = 1;
    int gblSuccess = 0;
    std::ostringstream errStrm;
    Teuchos::RCP<const Teuchos::Comm<int> > commPtr = A.getComm ();
    if (commPtr.is_null ()) {
      return lclNumErrs; // this process does not participate
    }
    const Teuchos::Comm<int>& comm = *commPtr;

    try {
      functor_type functor (lclNumErrs, diag, A);
    }
    catch (std::exception& e) {
      lclSuccess = -1;
      errStrm << "Process " << A.getComm ()->getRank () << ": "
              << e.what () << std::endl;
    }
    if (lclNumErrs != 0) {
      lclSuccess = 0;
    }

    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess == -1) {
      if (comm.getRank () == 0) {
        // We gather into std::cerr, rather than using an
        // std::ostringstream, because there might be a lot of MPI
        // processes.  It could take too much memory to gather all the
        // messages to Process 0 before printing.  gathervPrint gathers
        // and prints one message at a time, thus saving memory.  I
        // don't want to run out of memory while trying to print an
        // error message; that would hide the real problem.
        std::cerr << "getLocalDiagCopyWithoutOffsetsNotFillComplete threw an "
          "exception on one or more MPI processes in the matrix's comunicator."
                  << std::endl;
      }
      gathervPrint (std::cerr, errStrm.str (), comm);
      // Don't need to print anything here, since we've already
      // printed to std::cerr above.
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "");
    }
    else if (gblSuccess == 0) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (gblSuccess != 1, std::runtime_error,
         "getLocalDiagCopyWithoutOffsetsNotFillComplete failed on "
         "one or more MPI processes in the matrix's communicator.");
    }
  }
  else { // ! debug
    functor_type functor (lclNumErrs, diag, A);
  }

  return lclNumErrs;
}

} // namespace Details
} // namespace Tpetra

// Explicit template instantiation macro for
// getLocalDiagCopyWithoutOffsetsNotFillComplete.  NOT FOR USERS!!!
// Must be used inside the Tpetra namespace.
#define TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_INSTANT( SCALAR, LO, GO, NODE ) \
  template LO \
  Details::getLocalDiagCopyWithoutOffsetsNotFillComplete< SCALAR, LO, GO, NODE > \
    ( ::Tpetra::Vector< SCALAR, LO, GO, NODE >& diag, \
      const ::Tpetra::RowMatrix< SCALAR, LO, GO, NODE >& A, \
      const bool debug);

#endif // TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DEF_HPP
