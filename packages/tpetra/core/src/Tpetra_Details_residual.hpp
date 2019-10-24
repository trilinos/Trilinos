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
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RESIDUAL_HPP
#define TPETRA_DETAILS_RESIDUAL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"

/// \file Tpetra_Details_residual.hpp
/// \brief Functions that allow for fused residual calculation.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.-

namespace Tpetra {
  namespace Details {

template<class SC, class LO, class GO, class NO>
void localResidual(const CrsMatrix<SC,LO,GO,NO> &  A,
                   const MultiVector<SC,LO,GO,NO> & X,
                   const MultiVector<SC,LO,GO,NO> & B,
                   MultiVector<SC,LO,GO,NO> & R) {
  using Tpetra::Details::ProfilingRegion;
  using Teuchos::NO_TRANS;
  ProfilingRegion regionLocalApply ("Tpetra::CrsMatrix::localResidual");

  auto X_lcl = X.getLocalViewDevice ();
  auto B_lcl = B.getLocalViewDevice ();
  auto R_lcl = R.getLocalViewDevice ();
  auto lclMatrix_ = A.getLocalMatrix ();

  const bool debug = ::Tpetra::Details::Behavior::debug ();
  if (debug) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != R.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X.getNumVectors () << " != "
       "R.getNumVectors() = " << R.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != R.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X.getNumVectors () << " != "
       "R.getNumVectors() = " << R.getNumVectors () << ".");

    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getLocalLength () !=
       A.getColMap ()->getNodeNumElements (), std::runtime_error,
       "X has the wrong number of local rows.  "
       "X.getLocalLength() = " << X.getLocalLength () << " != "
       "A.getColMap()->getNodeNumElements() = " <<
       A.getColMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (R.getLocalLength () !=
       A.getRowMap ()->getNodeNumElements (), std::runtime_error,
       "R has the wrong number of local rows.  "
       "R.getLocalLength() = " << R.getLocalLength () << " != "
       "A.getRowMap()->getNodeNumElements() = " <<
       A.getRowMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.getLocalLength () !=
       A.getRowMap ()->getNodeNumElements (), std::runtime_error,
       "B has the wrong number of local rows.  "
       "B.getLocalLength() = " << B.getLocalLength () << " != "
       "A.getRowMap()->getNodeNumElements() = " <<
       A.getRowMap ()->getNodeNumElements () << ".");

    TEUCHOS_TEST_FOR_EXCEPTION
      (! A.isFillComplete (), std::runtime_error, "The matrix A is not "
       "fill complete.  You must call fillComplete() (possibly with "
       "domain and range Map arguments) without an intervening "
       "resumeFill() call before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! X.isConstantStride () || ! R.isConstantStride () || ! B.isConstantStride (),
       std::runtime_error, "X, Y and B must be constant stride.");
    // If the two pointers are NULL, then they don't alias one
    // another, even though they are equal.
    TEUCHOS_TEST_FOR_EXCEPTION
      ((X_lcl.data () == R_lcl.data () && X_lcl.data () != nullptr) ||
       (X_lcl.data () == B_lcl.data () && X_lcl.data () != nullptr),
       std::runtime_error, "X, Y and R may not alias one another.");
  }
      
  // This is currently a "reference implementation" waiting until Kokkos Kernels provides
  // a residual kernel.
  SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one, zero = Teuchos::ScalarTraits<SC>::zero();    
  A.localApply(X,R,Teuchos::NO_TRANS, one, zero);
  R.update(one,B,negone);
}
    

//! Computes R = B - A * X 
template<class SC, class LO, class GO, class NO>
void residual(const Operator<SC,LO,GO,NO> &   Aop,
              const MultiVector<SC,LO,GO,NO> & X_in,
              const MultiVector<SC,LO,GO,NO> & B_in,
              MultiVector<SC,LO,GO,NO> & R_in) {
  using Tpetra::Details::ProfilingRegion;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;
  
  const CrsMatrix<SC,LO,GO,NO> * Apt  = dynamic_cast<const CrsMatrix<SC,LO,GO,NO>*>(&Aop);
  if(!Apt) {
    // If we're not a CrsMatrix, we can't do fusion, so just do apply+update
     SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one, zero = Teuchos::ScalarTraits<SC>::zero();    
     Aop.apply(X_in,R_in,Teuchos::NO_TRANS, one, zero);
     R_in.update(one,B_in,negone);
     return;
  }
  const CrsMatrix<SC,LO,GO,NO> & A = *Apt;

  using import_type = typename CrsMatrix<SC,LO,GO,NO>::import_type;
  using export_type = typename CrsMatrix<SC,LO,GO,NO>::export_type;
  using MV = MultiVector<SC,LO,GO,NO>;

  // We treat the case of a replicated MV output specially.
  const bool R_is_replicated =
    (! R_in.isDistributed () && A.getComm ()->getSize () != 1);

  // It's possible that R is a view of X or B.  
  // We don't try to to detect the more subtle cases (e.g., one is a
  // subview of the other, but their initial pointers differ).  We
  // only need to do this if this matrix's Import is trivial;
  // otherwise, we don't actually apply the operator from X into Y.
  
  RCP<const import_type> importer = A.getGraph ()->getImporter ();
  RCP<const export_type> exporter = A.getGraph ()->getExporter ();

  // Temporary MV for Import operation.  After the block of code
  // below, this will be an (Imported if necessary) column Map MV
  // ready to give to localApply(...).
  RCP<const MV> X_colMap;
  if (importer.is_null ()) {
    if (! X_in.isConstantStride ()) {
      // Not all sparse mat-vec kernels can handle an input MV with
      // nonconstant stride correctly, so we have to copy it in that
      // case into a constant stride MV.  To make a constant stride
      // copy of X_in, we force creation of the column (== domain)
      // Map MV (if it hasn't already been created, else fetch the
      // cached copy).  This avoids creating a new MV each time.
      RCP<MV> X_colMapNonConst = A.getColumnMapMultiVector (X_in, true);
      Tpetra::deep_copy (*X_colMapNonConst, X_in);
      X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
    }
    else {
      // The domain and column Maps are the same, so do the local
      // multiply using the domain Map input MV X_in.
      X_colMap = rcpFromRef (X_in);
    }
  }
  else { // need to Import source (multi)vector
    ProfilingRegion regionImport ("Tpetra::CrsMatrix::residual: Import");
    // We're doing an Import anyway, which will copy the relevant
    // elements of the domain Map MV X_in into a separate column Map
    // MV.  Thus, we don't have to worry whether X_in is constant
    // stride.
    RCP<MV> X_colMapNonConst = A.getColumnMapMultiVector (X_in);
    
    // Import from the domain Map MV to the column Map MV.
    X_colMapNonConst->doImport (X_in, *importer, INSERT);
    X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
  }

  // Get a vector for the R_rowMap output residual, handling the 
  // non-constant stride and exporter cases.  Since R gets clobbered
  // we don't need to worry about the data in it
  RCP<MV> R_rowMap;
  if(exporter.is_null()) {
    if (! R_in.isConstantStride ()) {
      R_rowMap = A.getRowMapMultiVector(R_in);
    }
    else {
      R_rowMap = rcpFromRef (R_in);
    }
  }
  else {
    R_rowMap = A.getRowMapMultiVector (R_in);    
  }
  
  // Get a vector for the B_rowMap output residual, handling the 
  // non-constant stride and exporter cases
  RCP<const MV> B_rowMap;
  if(exporter.is_null()) {
    if (! B_in.isConstantStride ()) {
      // Do an allocation here.  If we need to optimize this later, we can have the matrix 
      // cache this.
      RCP<MV> B_rowMapNonConst = rcp(new MV(A.getRowMap(),B_in.getNumVectors()));
      Tpetra::deep_copy (*B_rowMapNonConst, B_in);
      B_rowMap = rcp_const_cast<const MV> (B_rowMapNonConst);
    }
    else {
      B_rowMap = rcpFromRef (B_in);
    }
  }
  else {
    // Do an allocation here.  If we need to optimize this later, we can have the matrix 
    // cache this.
    ProfilingRegion regionExport ("Tpetra::CrsMatrix::residual: B Import");
    RCP<MV> B_rowMapNonConst = rcp(new MV(A.getRowMap(),B_in.getNumVectors()));
    B_rowMapNonConst->doImport(B_in, *exporter, ADD);
    B_rowMap = rcp_const_cast<const MV> (B_rowMapNonConst);
  }

  // If we have a nontrivial Export object, we must perform an
  // Export.  In that case, the local multiply result will go into
  // the row Map multivector.  We don't have to make a
  // constant-stride version of R_in in this case, because we had to
  // make a constant stride R_rowMap MV and do an Export anyway.
  if (! exporter.is_null ()) {

    localResidual (A, *X_colMap, *B_rowMap, *R_rowMap);
    
    {
      ProfilingRegion regionExport ("Tpetra::CrsMatrix::residual: R Export");
      
      // Do the Export operation.
      R_in.doExport (*R_rowMap, *exporter, ADD);
    }
  }
  else { // Don't do an Export: row Map and range Map are the same.
    //
    // If R_in does not have constant stride,
    // then we can't let the kernel write directly to R_in.  
    // Instead, we have to use the cached row (== range)
    // Map MV as temporary storage.
    //
    if (! R_in.isConstantStride () ) {
      // We need to be sure to do a copy out in this case.
      localResidual (A, *X_colMap, *B_rowMap, *R_rowMap);
      Tpetra::deep_copy (R_in, *R_rowMap);
    }
    else {
      localResidual (A, *X_colMap, *B_rowMap, *R_rowMap);
    }
  }

  // If the range Map is a locally replicated Map, sum up
  // contributions from each process. 
  if (R_is_replicated) {
    ProfilingRegion regionReduce ("Tpetra::CrsMatrix::residual: Reduce Y");
    R_in.reduce ();
  }
}





  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_RESIDUAL_HPP
