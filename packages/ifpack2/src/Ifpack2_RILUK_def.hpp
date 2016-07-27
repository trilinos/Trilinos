/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_CRSRILUK_DEF_HPP
#define IFPACK2_CRSRILUK_DEF_HPP

#include "Ifpack2_LocalFilter.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"

#ifdef IFPACK2_ILUK_EXPERIMENTAL
#include <shylubasker_def.hpp>
#include <Kokkos_CrsMatrix.hpp>
# ifdef IFPACK2_HTS_EXPERIMENTAL
#  include <shylu_hts.hpp>
# endif
#endif

namespace Ifpack2 {

template<class MatrixType>
RILUK<MatrixType>::RILUK (const Teuchos::RCP<const row_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    LevelOfFill_ (0),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    isExperimental_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    initializeTime_ (0.0),
    computeTime_ (0.0),
    applyTime_ (0.0),
    RelaxValue_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ())
{}


template<class MatrixType>
RILUK<MatrixType>::RILUK (const Teuchos::RCP<const crs_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    LevelOfFill_ (0),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    isExperimental_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    initializeTime_ (0.0),
    computeTime_ (0.0),
    applyTime_ (0.0),
    RelaxValue_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ())
{}


template<class MatrixType>
RILUK<MatrixType>::~RILUK() {}


template<class MatrixType>
void
RILUK<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  if (A.getRawPtr () != A_.getRawPtr ()) {
    isAllocated_ = false;
    isInitialized_ = false;
    isComputed_ = false;
    A_local_crs_ = Teuchos::null;
    Graph_ = Teuchos::null;

    // The sparse triangular solvers get a triangular factor as their
    // input matrix.  The triangular factors L_ and U_ are getting
    // reset, so we reset the solvers' matrices to null.  Do that
    // before setting L_ and U_ to null, so that latter step actually
    // frees the factors.
    if (! L_solver_.is_null ()) {
      L_solver_->setMatrix (Teuchos::null);
    }
    if (! U_solver_.is_null ()) {
      U_solver_->setMatrix (Teuchos::null);
    }

    L_ = Teuchos::null;
    U_ = Teuchos::null;
    D_ = Teuchos::null;
    A_ = A;
  }
}



template<class MatrixType>
const typename RILUK<MatrixType>::crs_matrix_type&
RILUK<MatrixType>::getL () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    L_.is_null (), std::runtime_error, "Ifpack2::RILUK::getL: The L factor "
    "is null.  Please call initialize() (and preferably also compute()) "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *L_;
}


template<class MatrixType>
const Tpetra::Vector<typename RILUK<MatrixType>::scalar_type,
                     typename RILUK<MatrixType>::local_ordinal_type,
                     typename RILUK<MatrixType>::global_ordinal_type,
                     typename RILUK<MatrixType>::node_type>&
RILUK<MatrixType>::getD () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_.is_null (), std::runtime_error, "Ifpack2::RILUK::getD: The D factor "
    "(of diagonal entries) is null.  Please call initialize() (and "
    "preferably also compute()) before calling this method.  If the input "
    "matrix has not yet been set, you must first call setMatrix() with a "
    "nonnull input matrix before you may call initialize() or compute().");
  return *D_;
}


template<class MatrixType>
const typename RILUK<MatrixType>::crs_matrix_type&
RILUK<MatrixType>::getU () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    U_.is_null (), std::runtime_error, "Ifpack2::RILUK::getU: The U factor "
    "is null.  Please call initialize() (and preferably also compute()) "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *U_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename RILUK<MatrixType>::local_ordinal_type,
                               typename RILUK<MatrixType>::global_ordinal_type,
                               typename RILUK<MatrixType>::node_type> >
RILUK<MatrixType>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::RILUK::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  // FIXME (mfh 25 Jan 2014) Shouldn't this just come from A_?
  TEUCHOS_TEST_FOR_EXCEPTION(
    Graph_.is_null (), std::runtime_error, "Ifpack2::RILUK::getDomainMap: "
    "The computed graph is null.  Please call initialize() before calling "
    "this method.");
  return Graph_->getL_Graph ()->getDomainMap ();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename RILUK<MatrixType>::local_ordinal_type,
                               typename RILUK<MatrixType>::global_ordinal_type,
                               typename RILUK<MatrixType>::node_type> >
RILUK<MatrixType>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::RILUK::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  // FIXME (mfh 25 Jan 2014) Shouldn't this just come from A_?
  TEUCHOS_TEST_FOR_EXCEPTION(
    Graph_.is_null (), std::runtime_error, "Ifpack2::RILUK::getRangeMap: "
    "The computed graph is null.  Please call initialize() before calling "
    "this method.");
  return Graph_->getL_Graph ()->getRangeMap ();
}


template<class MatrixType>
void RILUK<MatrixType>::allocate_L_and_U ()
{
  using Teuchos::null;
  using Teuchos::rcp;

  if (! isAllocated_) {
    // Deallocate any existing storage.  This avoids storing 2x
    // memory, since RCP op= does not deallocate until after the
    // assignment.
    L_ = null;
    U_ = null;
    D_ = null;

    // Allocate Matrix using ILUK graphs
    L_ = rcp (new crs_matrix_type (Graph_->getL_Graph ()));
    U_ = rcp (new crs_matrix_type (Graph_->getU_Graph ()));
    L_->setAllToScalar (STS::zero ()); // Zero out L and U matrices
    U_->setAllToScalar (STS::zero ());

    // FIXME (mfh 24 Jan 2014) This assumes domain == range Map for L and U.
    L_->fillComplete ();
    U_->fillComplete ();
    D_ = rcp (new vec_type (Graph_->getL_Graph ()->getRowMap ()));
  }
  isAllocated_ = true;
}


template<class MatrixType>
void
RILUK<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  int fillLevel = 0;
  magnitude_type absThresh = STM::zero ();
  magnitude_type relThresh = STM::one ();
  magnitude_type relaxValue = STM::zero ();

  //
  // "fact: iluk level-of-fill" parsing is more complicated, because
  // we want to allow as many types as make sense.  int is the native
  // type, but we also want to accept magnitude_type (for
  // compatibility with ILUT) and double (for backwards compatibilty
  // with ILUT).
  //

  bool gotFillLevel = false;
  try {
    fillLevel = params.get<int> ("fact: iluk level-of-fill");
    gotFillLevel = true;
  }
  catch (InvalidParameterType&) {
    // Throwing again in the catch block would just unwind the stack.
    // Instead, we do nothing here, and check the Boolean outside to
    // see if we got the value.
  }
  catch (InvalidParameterName&) {
    gotFillLevel = true; // Accept the default value.
  }

  if (! gotFillLevel) {
    try {
      // Try global_ordinal_type.  The cast to int must succeed.
      fillLevel = as<int> (params.get<global_ordinal_type> ("fact: iluk level-of-fill"));
      gotFillLevel = true;
    }
    catch (InvalidParameterType&) {
      // Try the next type.
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  if (! gotFillLevel) {
    try {
      // Try magnitude_type, for compatibility with ILUT.
      // The cast from magnitude_type to int must succeed.
      fillLevel = as<int> (params.get<magnitude_type> ("fact: iluk level-of-fill"));
      gotFillLevel = true;
    }
    catch (InvalidParameterType&) {
      // Try the next type.
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  if (! gotFillLevel) {
    try {
      // Try double, for compatibility with ILUT.
      // The cast from double to int must succeed.
      fillLevel = as<int> (params.get<double> ("fact: iluk level-of-fill"));
      gotFillLevel = true;
    }
    catch (InvalidParameterType& e) {
      // We're out of options.  The user gave us the parameter, but it
      // doesn't have the right type.  The best thing for us to do in
      // that case is to throw, telling the user to use the right
      // type.
      throw e;
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! gotFillLevel,
    std::logic_error,
    "Ifpack2::RILUK::setParameters: We should never get here!  "
    "The method should either have read the \"fact: iluk level-of-fill\"  "
    "parameter by this point, or have thrown an exception.  "
    "Please let the Ifpack2 developers know about this bug.");

  //
  // For the other parameters, we prefer magnitude_type, but allow
  // double for backwards compatibility.
  //

  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double> ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double> ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relaxValue = params.get<magnitude_type> ("fact: relax value");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relaxValue = as<magnitude_type> (params.get<double> ("fact: relax value"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  //Experimental
  //Note: just following RILUK original style.
  //Do not think catch is the method for this. JDB
#ifdef IFPACK2_ILUK_EXPERIMENTAL
  try {
    isExperimental_ = params.get<bool> ("fact: iluk experimental basker");
  }
  catch (InvalidParameterType&) {
    //Use default
  }
  catch (InvalidParameterName&) {
    //Use default
  }

  try {
    basker_threads = params.get<int> ("fact: iluk experimental basker threads");
  }
  catch (InvalidParameterType&) {
    basker_threads = 1;
  }
  catch (InvalidParameterName&) {
    basker_threads = 1;
  }

  try {
    basker_user_fill = (scalar_type) params.get<double>("fact: iluk experimental basker user_fill");
  }
  catch (InvalidParameterType&) {
    basker_user_fill = (scalar_type) 0;
  }
  catch (InvalidParameterName&) {
    basker_user_fill = (scalar_type) 0;
  }

# ifdef IFPACK2_HTS_EXPERIMENTAL
  use_hts_ = false;
  hts_nthreads_ = basker_threads;
  {
    const std::string name("fact: iluk experimental HTS");
    if (params.isType<bool> (name))
      use_hts_ = params.get<bool> (name);
  }
  {
    const std::string name("fact: iluk experimental HTS threads");
    if (params.isType<int> (name))
      hts_nthreads_ = params.get<int> (name);
  }
# endif
#endif

  // "Commit" the values only after validating all of them.  This
  // ensures that there are no side effects if this routine throws an
  // exception.

  LevelOfFill_ = fillLevel;

  // mfh 28 Nov 2012: The previous code would not assign Athresh_,
  // Rthresh_, or RelaxValue_, if the read-in value was -1.  I don't
  // know if keeping this behavior is correct, but I'll keep it just
  // so as not to change previous behavior.

  if (absThresh != -STM::one ()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one ()) {
    Rthresh_ = relThresh;
  }
  if (relaxValue != -STM::one ()) {
    RelaxValue_ = relaxValue;
  }
}


template<class MatrixType>
Teuchos::RCP<const typename RILUK<MatrixType>::row_matrix_type>
RILUK<MatrixType>::getMatrix () const {
  return Teuchos::rcp_implicit_cast<const row_matrix_type> (A_);
}


template<class MatrixType>
Teuchos::RCP<const typename RILUK<MatrixType>::crs_matrix_type>
RILUK<MatrixType>::getCrsMatrix () const {
  return Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_, true);
}


template<class MatrixType>
Teuchos::RCP<const typename RILUK<MatrixType>::row_matrix_type>
RILUK<MatrixType>::makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;

  // If A_'s communicator only has one process, or if its column and
  // row Maps are the same, then it is already local, so use it
  // directly.
  if (A->getRowMap ()->getComm ()->getSize () == 1 ||
      A->getRowMap ()->isSameAs (* (A->getColMap ()))) {
    return A;
  }

  // If A_ is already a LocalFilter, then use it directly.  This
  // should be the case if RILUK is being used through
  // AdditiveSchwarz, for example.
  RCP<const LocalFilter<row_matrix_type> > A_lf_r =
    rcp_dynamic_cast<const LocalFilter<row_matrix_type> > (A);
  if (! A_lf_r.is_null ()) {
    return rcp_implicit_cast<const row_matrix_type> (A_lf_r);
  }
  else {
    // A_'s communicator has more than one process, its row Map and
    // its column Map differ, and A_ is not a LocalFilter.  Thus, we
    // have to wrap it in a LocalFilter.
    return rcp (new LocalFilter<row_matrix_type> (A));
  }
}


template<class MatrixType>
void RILUK<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  typedef Tpetra::CrsGraph<local_ordinal_type,
                           global_ordinal_type,
                           node_type> crs_graph_type;
  const char prefix[] = "Ifpack2::RILUK::initialize: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix is null.  Please "
     "call setMatrix() with a nonnull input before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_->isFillComplete (), std::runtime_error, prefix << "The matrix is not "
     "fill complete.  You may not invoke initialize() or compute() with this "
     "matrix until the matrix is fill complete.  If your matrix is a "
     "Tpetra::CrsMatrix, please call fillComplete on it (with the domain and "
     "range Maps, if appropriate) before calling this method.");

  Teuchos::Time timer ("RILUK::initialize");
  { // Start timing
    Teuchos::TimeMonitor timeMon (timer);

    // Calling initialize() means that the user asserts that the graph
    // of the sparse matrix may have changed.  We must not just reuse
    // the previous graph in that case.
    //
    // Regarding setting isAllocated_ to false: Eventually, we may want
    // some kind of clever memory reuse strategy, but it's always
    // correct just to blow everything away and start over.
    isInitialized_ = false;
    isAllocated_ = false;
    isComputed_ = false;
    Graph_ = Teuchos::null;

    RCP<const row_matrix_type> A_local = makeLocalFilter (A_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local.is_null (), std::logic_error, "Ifpack2::RILUK::initialize: "
      "makeLocalFilter returned null; it failed to compute A_local.  "
      "Please report this bug to the Ifpack2 developers.");

    // FIXME (mfh 24 Jan 2014, 26 Mar 2014) It would be more efficient
    // to rewrite RILUK so that it works with any RowMatrix input, not
    // just CrsMatrix.  (That would require rewriting IlukGraph to
    // handle a Tpetra::RowGraph.)  However, to make it work for now,
    // we just copy the input matrix if it's not a CrsMatrix.
    {
      RCP<const crs_matrix_type> A_local_crs =
        rcp_dynamic_cast<const crs_matrix_type> (A_local);
      if (A_local_crs.is_null ()) {
        // FIXME (mfh 24 Jan 2014) It would be smarter to count up the
        // number of elements in each row of A_local, so that we can
        // create A_local_crs_nc using static profile.  The code below is
        // correct but potentially slow.
        RCP<crs_matrix_type> A_local_crs_nc =
          rcp (new crs_matrix_type (A_local->getRowMap (),
                                    A_local->getColMap (), 0));
        // FIXME (mfh 24 Jan 2014) This Import approach will only work
        // if A_ has a one-to-one row Map.  This is generally the case
        // with matrices given to Ifpack2.
        //
        // Source and destination Maps are the same in this case.
        // That way, the Import just implements a copy.
        typedef Tpetra::Import<local_ordinal_type, global_ordinal_type,
          node_type> import_type;
        import_type import (A_local->getRowMap (), A_local->getRowMap ());
        A_local_crs_nc->doImport (*A_local, import, Tpetra::REPLACE);
        A_local_crs_nc->fillComplete (A_local->getDomainMap (), A_local->getRangeMap ());
        A_local_crs = rcp_const_cast<const crs_matrix_type> (A_local_crs_nc);
      }
      A_local_crs_ = A_local_crs;
      Graph_ = rcp (new Ifpack2::IlukGraph<crs_graph_type> (A_local_crs->getCrsGraph (),
                                                            LevelOfFill_, 0));
    }

    if(!isExperimental_)
      {
        Graph_->initialize ();
        allocate_L_and_U ();
        initAllValues (*A_local_crs_);
      }
    else
      {
#ifdef IFPACK2_ILUK_EXPERIMENTAL
        typedef typename node_type::device_type    kokkos_device;
        typedef typename kokkos_device::execution_space kokkos_exe;
        typedef typename crs_matrix_type::local_matrix_type kokkos_csr_matrix;

        static_assert( std::is_same< kokkos_exe,
                       Kokkos::OpenMP>::value,
                       "Kokkos node type not supported by experimental thread basker RILUK");


        myBasker = rcp( new BaskerNS::Basker<local_ordinal_type, scalar_type, Kokkos::OpenMP>);
        myBasker->Options.no_pivot     = true;
        myBasker->Options.transpose    = true; //CRS not CCS
        myBasker->Options.symmetric    = false;
        myBasker->Options.realloc      = true;
        myBasker->Options.btf          = false;
        myBasker->Options.incomplete   = true;
        myBasker->Options.inc_lvl      = LevelOfFill_;
        myBasker->Options.user_fill    = basker_user_fill;
        myBasker->Options.same_pattern = false;
        myBasker->SetThreads(basker_threads);
        basker_reuse = false;

        kokkos_csr_matrix kcsr = A_local_crs_->getLocalMatrix();

        local_ordinal_type* r_ptr;
        r_ptr = new local_ordinal_type[(local_ordinal_type)A_local_crs_->getNodeNumRows()+1];

        //Still need to convert
        //Lost on why Trilinos uses such odd types for row_ptrs
        for(local_ordinal_type bi = 0;
            bi < A_local_crs_->getNodeNumRows()+1; bi++)
          {
            r_ptr[bi] = (local_ordinal_type)kcsr.graph.row_map(bi);
          }

        int basker_error =
        myBasker->Symbolic(
         ((local_ordinal_type)A_local_crs_->getNodeNumRows()),
         ((local_ordinal_type)A_local_crs_->getNodeNumCols()),
         ((local_ordinal_type)A_local_crs_->getNodeNumEntries()),
         r_ptr,
         &(kcsr.graph.entries(0)),
         &(kcsr.values(0)));

        TEUCHOS_TEST_FOR_EXCEPTION(
        basker_error != 0, std::logic_error, "Ifpack2::RILUK::initialize:"
         "Error in basker Symbolic");


#else
      TEUCHOS_TEST_FOR_EXCEPTION(
      0==0, std::logic_error, "Ifpack2::RILUK::initialize: "
      "Using experimental ILUK without compiling experimental "
      "Try again with -DIFPACK2_ILUK_EXPERIMENAL.");
#endif
      }

  } // Stop timing

  isInitialized_ = true;
  ++numInitialize_;
  initializeTime_ += timer.totalElapsedTime ();
}


template<class MatrixType>
void
RILUK<MatrixType>::
initAllValues (const row_matrix_type& A)
{
  using Teuchos::ArrayRCP;
  using Teuchos::Comm;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> map_type;

  size_t NumIn = 0, NumL = 0, NumU = 0;
  bool DiagFound = false;
  size_t NumNonzeroDiags = 0;
  size_t MaxNumEntries = A.getGlobalMaxNumRowEntries();

  // First check that the local row map ordering is the same as the local portion of the column map.
  // The extraction of the strictly lower/upper parts of A, as well as the factorization,
  // implicitly assume that this is the case.
  Teuchos::ArrayView<const global_ordinal_type> rowGIDs = A.getRowMap()->getNodeElementList();
  Teuchos::ArrayView<const global_ordinal_type> colGIDs = A.getColMap()->getNodeElementList();
  bool gidsAreConsistentlyOrdered=true;
  global_ordinal_type indexOfInconsistentGID=0;
  for (global_ordinal_type i=0; i<rowGIDs.size(); ++i) {
    if (rowGIDs[i] != colGIDs[i]) {
      gidsAreConsistentlyOrdered=false;
      indexOfInconsistentGID=i;
      break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(gidsAreConsistentlyOrdered==false, std::runtime_error,
                             "The ordering of the local GIDs in the row and column maps is not the same"
                             << std::endl << "at index " << indexOfInconsistentGID
                             << ".  Consistency is required, as all calculations are done with"
                             << std::endl << "local indexing.");

  // Allocate temporary space for extracting the strictly
  // lower and upper parts of the matrix A.
  Teuchos::Array<local_ordinal_type> InI(MaxNumEntries);
  Teuchos::Array<local_ordinal_type> LI(MaxNumEntries);
  Teuchos::Array<local_ordinal_type> UI(MaxNumEntries);
  Teuchos::Array<scalar_type> InV(MaxNumEntries);
  Teuchos::Array<scalar_type> LV(MaxNumEntries);
  Teuchos::Array<scalar_type> UV(MaxNumEntries);

  // Check if values should be inserted or replaced
  const bool ReplaceValues = L_->isStaticGraph () || L_->isLocallyIndexed ();

  L_->resumeFill ();
  U_->resumeFill ();
  if (ReplaceValues) {
    L_->setAllToScalar (STS::zero ()); // Zero out L and U matrices
    U_->setAllToScalar (STS::zero ());
  }

  D_->putScalar (STS::zero ()); // Set diagonal values to zero
  ArrayRCP<scalar_type> DV = D_->get1dViewNonConst (); // Get view of diagonal

  RCP<const map_type> rowMap = L_->getRowMap ();

  // First we copy the user's matrix into L and U, regardless of fill level.
  // It is important to note that L and U are populated using local indices.
  // This means that if the row map GIDs are not monotonically increasing
  // (i.e., permuted or gappy), then the strictly lower (upper) part of the
  // matrix is not the one that you would get if you based L (U) on GIDs.
  // This is ok, as the *order* of the GIDs in the rowmap is a better
  // expression of the user's intent than the GIDs themselves.

  Teuchos::ArrayView<const global_ordinal_type> nodeGIDs = rowMap->getNodeElementList();
  for (size_t myRow=0; myRow<A.getNodeNumRows(); ++myRow) {
    local_ordinal_type local_row = myRow;

    //TODO JJH 4April2014 An optimization is to use getLocalRowView.  Not all matrices support this,
    //                    we'd need to check via the Tpetra::RowMatrix method supportsRowViews().
    A.getLocalRowCopy (local_row, InI(), InV(), NumIn); // Get Values and Indices

    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;

    for (size_t j = 0; j < NumIn; ++j) {
      const local_ordinal_type k = InI[j];

      if (k == local_row) {
        DiagFound = true;
        // Store perturbed diagonal in Tpetra::Vector D_
        DV[local_row] += Rthresh_ * InV[j] + IFPACK2_SGN(InV[j]) * Athresh_;
      }
      else if (k < 0) { // Out of range
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Ifpack2::RILUK::initAllValues: current "
          "GID k = " << k << " < 0.  I'm not sure why this is an error; it is "
          "probably an artifact of the undocumented assumptions of the "
          "original implementation (likely copied and pasted from Ifpack).  "
          "Nevertheless, the code I found here insisted on this being an error "
          "state, so I will throw an exception here.");
      }
      else if (k < local_row) {
        LI[NumL] = k;
        LV[NumL] = InV[j];
        NumL++;
      }
      else if (Teuchos::as<size_t>(k) <= rowMap->getNodeNumElements()) {
        UI[NumU] = k;
        UV[NumU] = InV[j];
        NumU++;
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) {
      ++NumNonzeroDiags;
    } else {
      DV[local_row] = Athresh_;
    }

    if (NumL) {
      if (ReplaceValues) {
        L_->replaceLocalValues(local_row, LI(0, NumL), LV(0,NumL));
      } else {
        //FIXME JJH 24April2014 Is this correct?  I believe this case is when there aren't already values
        //FIXME in this row in the column locations corresponding to UI.
        L_->insertLocalValues(local_row, LI(0,NumL), LV(0,NumL));
      }
    }

    if (NumU) {
      if (ReplaceValues) {
        U_->replaceLocalValues(local_row, UI(0,NumU), UV(0,NumU));
      } else {
        //FIXME JJH 24April2014 Is this correct?  I believe this case is when there aren't already values
        //FIXME in this row in the column locations corresponding to UI.
        U_->insertLocalValues(local_row, UI(0,NumU), UV(0,NumU));
      }
    }
  }

  // The domain of L and the range of U are exactly their own row maps
  // (there is no communication).  The domain of U and the range of L
  // must be the same as those of the original matrix, However if the
  // original matrix is a VbrMatrix, these two latter maps are
  // translation from a block map to a point map.
  L_->fillComplete (L_->getColMap (), A_local_crs_->getRangeMap ());
  U_->fillComplete (A_local_crs_->getDomainMap (), U_->getRowMap ());

  // At this point L and U have the values of A in the structure of L
  // and U, and diagonal vector D

  isInitialized_ = true;
}


template<class MatrixType>
void RILUK<MatrixType>::compute ()
{
  using Teuchos::rcp;
  const char prefix[] = "Ifpack2::RILUK::compute: ";

  // initialize() checks this too, but it's easier for users if the
  // error shows them the name of the method that they actually
  // called, rather than the name of some internally called method.
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix is null.  Please "
     "call setMatrix() with a nonnull input before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_->isFillComplete (), std::runtime_error, prefix << "The matrix is not "
     "fill complete.  You may not invoke initialize() or compute() with this "
     "matrix until the matrix is fill complete.  If your matrix is a "
     "Tpetra::CrsMatrix, please call fillComplete on it (with the domain and "
     "range Maps, if appropriate) before calling this method.");

  if (! isInitialized ()) {
    initialize (); // Don't count this in the compute() time
  }

  Teuchos::Time timer ("RILUK::compute");
  { //Start timing
    if(! isExperimental_) {
    isComputed_ = false;

    L_->resumeFill ();
    U_->resumeFill ();

    // MinMachNum should be officially defined, for now pick something a little
    // bigger than IEEE underflow value

    const scalar_type MinDiagonalValue = STS::rmin ();
    const scalar_type MaxDiagonalValue = STS::one () / MinDiagonalValue;

    size_t NumIn, NumL, NumU;

    // Get Maximum Row length
    const size_t MaxNumEntries =
      L_->getNodeMaxNumRowEntries () + U_->getNodeMaxNumRowEntries () + 1;

    Teuchos::Array<local_ordinal_type> InI(MaxNumEntries); // Allocate temp space
    Teuchos::Array<scalar_type> InV(MaxNumEntries);
    size_t num_cols = U_->getColMap()->getNodeNumElements();
    Teuchos::Array<int> colflag(num_cols);

    Teuchos::ArrayRCP<scalar_type> DV = D_->get1dViewNonConst(); // Get view of diagonal

    // Now start the factorization.

    // Need some integer workspace and pointers
    size_t NumUU;
    Teuchos::ArrayView<const local_ordinal_type> UUI;
    Teuchos::ArrayView<const scalar_type> UUV;
    for (size_t j = 0; j < num_cols; ++j) {
      colflag[j] = -1;
    }

    for (size_t i = 0; i < L_->getNodeNumRows (); ++i) {
      local_ordinal_type local_row = i;

      // Fill InV, InI with current row of L, D and U combined

      NumIn = MaxNumEntries;
      L_->getLocalRowCopy (local_row, InI (), InV (), NumL);

      InV[NumL] = DV[i]; // Put in diagonal
      InI[NumL] = local_row;

      U_->getLocalRowCopy (local_row, InI (NumL+1, MaxNumEntries-NumL-1),
                           InV (NumL+1, MaxNumEntries-NumL-1), NumU);
      NumIn = NumL+NumU+1;

      // Set column flags
      for (size_t j = 0; j < NumIn; ++j) {
        colflag[InI[j]] = j;
      }

      scalar_type diagmod = STS::zero (); // Off-diagonal accumulator

      for (size_t jj = 0; jj < NumL; ++jj) {
        local_ordinal_type j = InI[jj];
        scalar_type multiplier = InV[jj]; // current_mults++;

        InV[jj] *= DV[j];

        U_->getLocalRowView(j, UUI, UUV); // View of row above
        NumUU = UUI.size();

        if (RelaxValue_ == STM::zero ()) {
          for (size_t k = 0; k < NumUU; ++k) {
            const int kk = colflag[UUI[k]];
            // FIXME (mfh 23 Dec 2013) Wait a second, we just set
            // colflag above using size_t (which is generally unsigned),
            // but now we're querying it using int (which is signed).
            if (kk > -1) {
              InV[kk] -= multiplier * UUV[k];
            }
          }
        }
        else {
          for (size_t k = 0; k < NumUU; ++k) {
            // FIXME (mfh 23 Dec 2013) Wait a second, we just set
            // colflag above using size_t (which is generally unsigned),
            // but now we're querying it using int (which is signed).
            const int kk = colflag[UUI[k]];
            if (kk > -1) {
              InV[kk] -= multiplier*UUV[k];
            }
            else {
              diagmod -= multiplier*UUV[k];
            }
          }
        }
      }
      if (NumL) {
        // Replace current row of L
        L_->replaceLocalValues (local_row, InI (0, NumL), InV (0, NumL));
      }

      DV[i] = InV[NumL]; // Extract Diagonal value

      if (RelaxValue_ != STM::zero ()) {
        DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
      }

      if (STS::magnitude (DV[i]) > STS::magnitude (MaxDiagonalValue)) {
        if (STS::real (DV[i]) < STM::zero ()) {
          DV[i] = -MinDiagonalValue;
        }
        else {
          DV[i] = MinDiagonalValue;
        }
      }
      else {
        DV[i] = STS::one () / DV[i]; // Invert diagonal value
      }

      for (size_t j = 0; j < NumU; ++j) {
        InV[NumL+1+j] *= DV[i]; // Scale U by inverse of diagonal
      }

      if (NumU) {
        // Replace current row of L and U
        U_->replaceLocalValues (local_row, InI (NumL+1, NumU), InV (NumL+1, NumU));
      }

      // Reset column flags
      for (size_t j = 0; j < NumIn; ++j) {
        colflag[InI[j]] = -1;
      }
    }

    // FIXME (mfh 23 Dec 2013) Do we know that the column Map of L_ is
    // always one-to-one?
    L_->fillComplete (L_->getColMap (), A_local_crs_->getRangeMap ());
    U_->fillComplete (A_local_crs_->getDomainMap (), U_->getRowMap ());

    // Validate that the L and U factors are actually lower and upper triangular
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! L_->isLowerTriangular (), std::runtime_error,
      "Ifpack2::RILUK::compute: L isn't lower triangular.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! U_->isUpperTriangular (), std::runtime_error,
      "Ifpack2::RILUK::compute: U isn't lower triangular.");

    L_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> (L_));
    L_solver_->initialize ();
    L_solver_->compute ();

    U_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> (U_));
    U_solver_->initialize ();
    U_solver_->compute ();
  }
  else {
#ifdef IFPACK2_ILUK_EXPERIMENTAL
    Teuchos::Time timerbasker ("RILUK::basercompute");
    {
      //
      if(basker_reuse == false)
        {
          //We don't have to reimport Matrix because same
          {
            //Teuchos::TimeMonitor timeMon(timerbasker);
            myBasker->Factor_Inc(0);
            basker_reuse = true;
          }
        }
      else
        {
          //Do we need to import Matrix with file again?
          myBasker->Options.same_pattern = true;

          typedef typename crs_matrix_type::local_matrix_type kokkos_csr_matrix;
          kokkos_csr_matrix kcsr = A_local_crs_->getLocalMatrix();


         local_ordinal_type* r_ptr;
        r_ptr = new local_ordinal_type[(local_ordinal_type)A_local_crs_->getNodeNumRows()+1];

        //Still need to convert
        //Lost on why Trilinos uses such odd types for row_ptrs
        for(local_ordinal_type bi = 0; bi < A_local_crs_->getNodeNumRows()+1; bi++)
          {
            r_ptr[bi] = (local_ordinal_type)kcsr.graph.row_map(bi);
          }

        int basker_error =
        myBasker->Factor(
         ((local_ordinal_type)A_local_crs_->getNodeNumRows()),
         ((local_ordinal_type)A_local_crs_->getNodeNumCols()),
         ((local_ordinal_type)A_local_crs_->getNodeNumEntries()),
         r_ptr,
         &(kcsr.graph.entries(0)),
         &(kcsr.values(0)));

        //myBasker->DEBUG_PRINT();

        TEUCHOS_TEST_FOR_EXCEPTION(
        basker_error != 0, std::logic_error, "Ifpack2::RILUK::initialize:"
         "Error in basker compute");


        }
    }
# ifdef IFPACK2_HTS_EXPERIMENTAL
    //TODO Reuse symbolic information.
    if (use_hts_) {

      Teuchos::Time basker_getL("basker_getL");
      Teuchos::Time hts_buildL ("hts_buildL");
      Teuchos::Time basker_getU("basker_getU");
      Teuchos::Time hts_buildU ("hts_buildU");


      hts_mgr_ = Teuchos::rcp(new HTSManager());
      local_ordinal_type* p, * q;
      local_ordinal_type nnz;



      myBasker->GetPerm(&p, &q);

      {

        HTSData d;
        myBasker->GetL(d.n, nnz, &d.ir, &d.jc, &d.v);
        d.sort();
        typename HTST::CrsMatrix* T = HTST::make_CrsMatrix(d.n, d.ir, d.jc, d.v, true);
        hts_mgr_->Limpl = HTST::preprocess(T, 1, hts_nthreads_, true, p, 0);
        delete[] p;

      }

      {
        HTSData d;
        myBasker->GetU(d.n, nnz, &d.ir, &d.jc, &d.v);
        d.sort();

        typename HTST::CrsMatrix* T = HTST::make_CrsMatrix(d.n, d.ir, d.jc, d.v, true);
        hts_mgr_->Uimpl = HTST::preprocess(T, 1, hts_nthreads_, true, 0, q);


      }


    }
# endif

#else
    TEUCHOS_TEST_FOR_EXCEPTION(
       0==1, std::runtime_error,
       "Ifpack2::RILUK::compute: experimental not enabled");
#endif
  }//end -- if experimental
  }//end timing

  isComputed_ = true;
  ++numCompute_;
  computeTime_ += timer.totalElapsedTime ();
}


template<class MatrixType>
void
RILUK<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::RILUK::apply: The matrix is "
    "null.  Please call setMatrix() with a nonnull input, then initialize() "
    "and compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::RILUK::apply: If you have not yet called compute(), "
    "you must call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::RILUK::apply: X and Y do not have the same number of columns.  "
    "X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isComplex && mode == Teuchos::CONJ_TRANS, std::logic_error,
    "Ifpack2::RILUK::apply: mode = Teuchos::CONJ_TRANS is not implemented for "
    "complex Scalar type.  Please talk to the Ifpack2 developers to get this "
    "fixed.  There is a FIXME in this file about this very issue.");
#ifdef HAVE_IFPACK2_DEBUG
  {
    const magnitude_type D_nrm1 = D_->norm1 ();
    TEUCHOS_TEST_FOR_EXCEPTION( STM::isnaninf (D_nrm1), std::runtime_error, "Ifpack2::RILUK::apply: The 1-norm of the stored diagonal is NaN or Inf.");
    Teuchos::Array<magnitude_type> norms (X.getNumVectors ());
    X.norm1 (norms ());
    bool good = true;
    for (size_t j = 0; j < X.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION( ! good, std::runtime_error, "Ifpack2::RILUK::apply: The 1-norm of the input X is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  const scalar_type one = STS::one ();
  const scalar_type zero = STS::zero ();

  Teuchos::Time timer ("RILUK::apply");
  { // Start timing
    if(!isExperimental_){
    Teuchos::TimeMonitor timeMon (timer);
    if (alpha == one && beta == zero) {
      if (mode == Teuchos::NO_TRANS) { // Solve L (D (U Y)) = X for Y.
        // Start by solving L C = X for C.  C must have the same Map
        // as D.  We have to use a temp multivector, since our
        // implementation of triangular solves does not allow its
        // input and output to alias one another.
        //
        // FIXME (mfh 24 Jan 2014) Cache this temp multivector.
        MV C (D_->getMap (), X.getNumVectors ());
        L_solver_->apply (X, C, mode);

        // Solve D Y_tmp = C.  Y_tmp must have the same Map as C, and
        // the operation lets us do this in place in C, so we can
        // write "solve D C = C for C."
        C.elementWiseMultiply (one, *D_, C, zero);

        U_solver_->apply (C, Y, mode); // Solve U Y = C.
      }
      else { // Solve U^P (D^P (U^P Y)) = X for Y (where P is * or T).

        // Start by solving U^P C = X for C.  C must have the same Map
        // as D.  We have to use a temp multivector, since our
        // implementation of triangular solves does not allow its
        // input and output to alias one another.
        //
        // FIXME (mfh 24 Jan 2014) Cache this temp multivector.
        MV C (D_->getMap (), X.getNumVectors ());
        U_solver_->apply (X, C, mode);

        // Solve D^P Y_tmp = C.  Y_tmp must have the same Map as C,
        // and the operation lets us do this in place in C, so we can
        // write "solve D^P C = C for C."
        //
        // FIXME (mfh 24 Jan 2014) If mode = Teuchos::CONJ_TRANS, we
        // need to do an elementwise multiply with the conjugate of
        // D_, not just with D_ itself.
        C.elementWiseMultiply (one, *D_, C, zero);

        L_solver_->apply (C, Y, mode); // Solve L^P Y = C.
      }
    }
    else { // alpha != 1 or beta != 0
      if (alpha == zero) {
        if (beta == zero) {
          Y.putScalar (zero);
        } else {
          Y.scale (beta);
        }
      } else { // alpha != zero
        MV Y_tmp (Y.getMap (), Y.getNumVectors ());
        apply (X, Y_tmp, mode);
        Y.update (alpha, Y_tmp, beta);
      }
    }
  }
  else
    {
#ifdef IFPACK2_ILUK_EXPERIMENTAL
      Teuchos::ArrayRCP<const scalar_type> XX;
      Teuchos::ArrayRCP<const scalar_type> YY;
      XX = X.get1dView();
      YY = Y.get1dView();

# ifdef IFPACK2_HTS_EXPERIMENTAL
      if (use_hts_) {
        HTST::solve_omp(hts_mgr_->Limpl, const_cast<scalar_type*>(XX.getRawPtr()),
                        X.getNumVectors(), const_cast<scalar_type*>(YY.getRawPtr()));
        HTST::solve_omp(hts_mgr_->Uimpl, const_cast<scalar_type*>(YY.getRawPtr()),
                        X.getNumVectors());
      } else
# endif
      {
        myBasker->Solve(((local_ordinal_type)X.getNumVectors()),
                        (const_cast<scalar_type *>(XX.getRawPtr())),
                        (const_cast<scalar_type *>(YY.getRawPtr())));
      }
#else
      TEUCHOS_TEST_FOR_EXCEPTION(
      0==1, std::runtime_error,
      "Ifpack2::RILUK::apply: Experimental no enabled");
#endif
    }//end isExperimental
  }//end timing

#ifdef HAVE_IFPACK2_DEBUG
  {
    Teuchos::Array<magnitude_type> norms (Y.getNumVectors ());
    Y.norm1 (norms ());
    bool good = true;
    for (size_t j = 0; j < Y.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION( ! good, std::runtime_error, "Ifpack2::RILUK::apply: The 1-norm of the output Y is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  ++numApply_;
  applyTime_ = timer.totalElapsedTime ();
}


template<class MatrixType>
void RILUK<MatrixType>::
multiply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
          const Teuchos::ETransp mode) const
{
  const scalar_type zero = STS::zero ();
  const scalar_type one = STS::one ();

  if (mode != Teuchos::NO_TRANS) {
    U_->apply (X, Y, mode); //
    Y.update (one, X, one); // Y = Y + X (account for implicit unit diagonal)

    // FIXME (mfh 24 Jan 2014) If mode = Teuchos::CONJ_TRANS, we need
    // to do an elementwise multiply with the conjugate of D_, not
    // just with D_ itself.
    Y.elementWiseMultiply (one, *D_, Y, zero); // y = D*y (D_ has inverse of diagonal)

    MV Y_tmp (Y, Teuchos::Copy); // Need a temp copy of Y
    L_->apply (Y_tmp, Y, mode);
    Y.update (one, Y_tmp, one); // (account for implicit unit diagonal)
  }
  else {
    L_->apply (X, Y, mode);
    Y.update (one, X, one); // Y = Y + X (account for implicit unit diagonal)
    Y.elementWiseMultiply (one, *D_, Y, zero); // y = D*y (D_ has inverse of diagonal)
    MV Y_tmp (Y, Teuchos::Copy); // Need a temp copy of Y1
    U_->apply (Y_tmp, Y, mode);
    Y.update (one, Y_tmp, one); // (account for implicit unit diagonal)
  }
}


template<class MatrixType>
std::string RILUK<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::RILUK\": {";
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  os << "Level-of-fill: " << getLevelOfFill() << ", ";

  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]"
       << ", Global nnz: " << A_->getGlobalNumEntries();
  }

  os << "}";
  return os.str ();
}

#if defined IFPACK2_ILUK_EXPERIMENTAL && defined IFPACK2_HTS_EXPERIMENTAL
template<class MatrixType>
RILUK<MatrixType>::HTSManager::HTSManager ()
  : Limpl(0), Uimpl(0)
{}

template<class MatrixType>
RILUK<MatrixType>::HTSManager::~HTSManager ()
{
  HTST::delete_Impl(Limpl);
  HTST::delete_Impl(Uimpl);
}

template<class MatrixType>
RILUK<MatrixType>::HTSData::HTSData ()
  : jc(0), ir(0), v(0)
{}

template<class MatrixType>
RILUK<MatrixType>::HTSData::~HTSData ()
{
  free_CrsMatrix_data();
}

template<class MatrixType>
void RILUK<MatrixType>::HTSData::free_CrsMatrix_data ()
{
  if (jc) delete[] jc;
  if (ir) delete[] ir;
  if (v) delete[] v;
  jc = ir = 0;
  v = 0;
}

template<class MatrixType>
void RILUK<MatrixType>::HTSData::sort ()
{
  if ( ! ir || ! jc) return;
  std::vector<Entry> es;
  for (local_ordinal_type i = 0; i < n; ++i) {
    es.resize(ir[i+1] - ir[i]);
    const local_ordinal_type os = ir[i];
    for (local_ordinal_type j = 0; j < ir[i+1] - os; ++j) {
      es[j].j = jc[os+j];
      es[j].v = v[os+j];
    }
    std::sort(es.begin(), es.end());
    for (local_ordinal_type j = 0; j < ir[i+1] - os; ++j) {
      jc[os+j] = es[j].j;
      v[os+j] = es[j].v;
    }
  }
}
#endif

} // namespace Ifpack2

#define IFPACK2_RILUK_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::RILUK< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
