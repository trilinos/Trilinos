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

#ifndef TPETRA_ROWMATRIX_DEF_HPP
#define TPETRA_ROWMATRIX_DEF_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~RowMatrix() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  add (const Scalar& alpha,
       const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
       const Scalar& beta,
       const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
       const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    using Teuchos::sublist;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Map<LO, GO, Node> map_type;
    typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

    const this_type& B = *this; // a convenient abbreviation

    // If the user didn't supply a domain or range Map, then try to
    // get one from B first (if it has them), then from A (if it has
    // them) If we don't have any domain or range Maps, scold the
    // user.
    RCP<const map_type> A_domainMap = A.getDomainMap ();
    RCP<const map_type> A_rangeMap = A.getRangeMap ();
    RCP<const map_type> B_domainMap = B.getDomainMap ();
    RCP<const map_type> B_rangeMap = B.getRangeMap ();

    RCP<const map_type> theDomainMap = domainMap;
    RCP<const map_type> theRangeMap = rangeMap;

    if (domainMap.is_null ()) {
      if (B_domainMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_domainMap.is_null (), std::invalid_argument,
          "Tpetra::RowMatrix::add: If neither A nor B have a domain Map, "
          "then you must supply a nonnull domain Map to this method.");
        theDomainMap = A_domainMap;
      } else {
        theDomainMap = B_domainMap;
      }
    }
    if (rangeMap.is_null ()) {
      if (B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_rangeMap.is_null (), std::invalid_argument,
          "Tpetra::RowMatrix::add: If neither A nor B have a range Map, "
          "then you must supply a nonnull range Map to this method.");
        theRangeMap = A_rangeMap;
      } else {
        theRangeMap = B_rangeMap;
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, check that A and B have matching domain and
    // range Maps, if they have domain and range Maps at all.  (If
    // they aren't fill complete, then they may not yet have them.)
    if (! A_domainMap.is_null () && ! A_rangeMap.is_null ()) {
      if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_domainMap->isSameAs (*A_domainMap), std::invalid_argument,
          "Tpetra::RowMatrix::add: The input RowMatrix A must have a domain Map "
          "which is the same as (isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_rangeMap->isSameAs (*A_rangeMap), std::invalid_argument,
          "Tpetra::RowMatrix::add: The input RowMatrix A must have a range Map "
          "which is the same as (isSameAs) this RowMatrix's range Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
          std::invalid_argument,
          "Tpetra::RowMatrix::add: The input domain Map must be the same as "
          "(isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
          std::invalid_argument,
          "Tpetra::RowMatrix::add: The input range Map must be the same as "
          "(isSameAs) this RowMatrix's range Map.");
      }
    }
    else if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
        std::invalid_argument,
        "Tpetra::RowMatrix::add: The input domain Map must be the same as "
        "(isSameAs) this RowMatrix's domain Map.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
        std::invalid_argument,
        "Tpetra::RowMatrix::add: The input range Map must be the same as "
        "(isSameAs) this RowMatrix's range Map.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        domainMap.is_null () || rangeMap.is_null (), std::invalid_argument,
        "Tpetra::RowMatrix::add: If neither A nor B have a domain and range "
        "Map, then you must supply a nonnull domain and range Map to this "
        "method.");
    }
#endif // HAVE_TPETRA_DEBUG

    // What parameters do we pass to C's constructor?  Do we call
    // fillComplete on C after filling it?  And if so, what parameters
    // do we pass to C's fillComplete call?
    bool callFillComplete = true;
    RCP<ParameterList> constructorSublist;
    RCP<ParameterList> fillCompleteSublist;
    if (! params.is_null ()) {
      callFillComplete = params->get ("Call fillComplete", callFillComplete);
      constructorSublist = sublist (params, "Constructor parameters");
      fillCompleteSublist = sublist (params, "fillComplete parameters");
    }

    RCP<const map_type> A_rowMap = A.getRowMap ();
    RCP<const map_type> B_rowMap = B.getRowMap ();
    RCP<const map_type> C_rowMap = B_rowMap; // see discussion in documentation
    RCP<crs_matrix_type> C; // The result matrix.

    // If A and B's row Maps are the same, we can compute an upper
    // bound on the number of entries in each row of C, before
    // actually computing the sum.  A reasonable upper bound is the
    // sum of the two entry counts in each row.  If we choose this as
    // the actual per-row upper bound, we can use static profile.
    if (A_rowMap->isSameAs (*B_rowMap)) {
      const LO localNumRows = as<LO> (A_rowMap->getNodeNumElements ());
      ArrayRCP<size_t> C_maxNumEntriesPerRow (localNumRows, 0);

      // Get the number of entries in each row of A.
      if (alpha != STS::zero ()) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += A_numEntries;
        }
      }
      // Get the number of entries in each row of B.
      if (beta != STS::zero ()) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += B_numEntries;
        }
      }
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile, constructorSublist));
      }
      // Since A and B have the same row Maps, we could add them
      // together all at once and merge values before we call
      // insertGlobalValues.  However, we don't really need to, since
      // we've already allocated enough space in each row of C for C
      // to do the merge itself.
    }
    else { // the row Maps of A and B are not the same
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile,
                                      constructorSublist));
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
      "Tpetra::RowMatrix::add: C should not be null at this point.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    //
    // Compute C = alpha*A + beta*B.
    //
    Array<GO> ind;
    Array<Scalar> val;

    if (alpha != STS::zero ()) {
      const LO A_localNumRows = as<LO> (A_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > as<size_t> (ind.size ())) {
          ind.resize (A_numEntries);
          val.resize (A_numEntries);
        }
        ArrayView<GO> indView = ind (0, A_numEntries);
        ArrayView<Scalar> valView = val (0, A_numEntries);
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != STS::one ()) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (beta != STS::zero ()) {
      const LO B_localNumRows = as<LO> (B_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > as<size_t> (ind.size ())) {
          ind.resize (B_numEntries);
          val.resize (B_numEntries);
        }
        ArrayView<GO> indView = ind (0, B_numEntries);
        ArrayView<Scalar> valView = val (0, B_numEntries);
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != STS::one ()) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (callFillComplete) {
      if (fillCompleteSublist.is_null ()) {
        C->fillComplete (theDomainMap, theRangeMap);
      } else {
        C->fillComplete (theDomainMap, theRangeMap, fillCompleteSublist);
      }
    }

    return rcp_implicit_cast<this_type> (C);
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrix< SCALAR , LO , GO , NODE >;


#endif // TPETRA_ROWMATRIX_DEF_HPP

