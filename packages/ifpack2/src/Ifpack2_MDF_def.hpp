// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_MDF_DEF_HPP
#define IFPACK2_MDF_DEF_HPP

#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Ifpack2_Details_getParamTryingTypes.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "KokkosKernels_Sorting.hpp"
#include <exception>
#include <type_traits>

namespace Ifpack2 {

namespace Details {

namespace MDFImpl
{

template<class dev_view_t>
auto copy_view(const dev_view_t & vals)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  typename dev_view_t::non_const_type newvals (view_alloc (vals.label(), WithoutInitializing), vals.extent (0));
  Kokkos::deep_copy(newvals,vals);
  return newvals;
}

template<class array_t,class dev_view_t>
void copy_dev_view_to_host_array(array_t & array, const dev_view_t & dev_view)
{
  using host_view_t = typename dev_view_t::HostMirror;

  // Clear out existing and allocate
  const auto ext = dev_view.extent(0);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ext != size_t(array.size()), std::logic_error, "Ifpack2::MDF::copy_dev_view_to_host_array: "
    "Size of permuations on host and device do not match.  "
    "Please report this bug to the Ifpack2 developers.");

  //Wrap array data in view and copy
  Kokkos::deep_copy(host_view_t(array.get(),ext),dev_view);
}

template<class scalar_type,class local_ordinal_type,class global_ordinal_type,class node_type>
void applyReorderingPermutations(
  const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
  Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
  const Teuchos::ArrayRCP<local_ordinal_type> & perm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
                             "Ifpack2::MDF::applyReorderingPermuations ERROR: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> >       y_ptr = Y.get2dViewNonConst();

  for(size_t k=0; k < X.getNumVectors(); k++)
    for(local_ordinal_type i=0; (size_t)i< X.getLocalLength(); i++)
      y_ptr[k][perm[i]] = x_ptr[k][i];
}


template<class scalar_type,class local_ordinal_type,class global_ordinal_type,class node_type>
auto get_local_crs_row_matrix(
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type>> A_local)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;

  using crs_matrix_type = Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;

  using nonconst_local_inds_host_view_type = typename crs_matrix_type::nonconst_local_inds_host_view_type;
  using nonconst_values_host_view_type =  typename crs_matrix_type::nonconst_values_host_view_type;

  RCP<const crs_matrix_type> A_local_crs = rcp_dynamic_cast<const crs_matrix_type>(A_local);
  if (A_local_crs.is_null ()) {
    local_ordinal_type numRows = A_local->getLocalNumRows();
    Array<size_t> entriesPerRow(numRows);
    for(local_ordinal_type i = 0; i < numRows; i++) {
      entriesPerRow[i] = A_local->getNumEntriesInLocalRow(i);
    }
    RCP<crs_matrix_type> A_local_crs_nc =
      rcp (new crs_matrix_type (A_local->getRowMap (),
                                A_local->getColMap (),
                                entriesPerRow()));
    // copy entries into A_local_crs
    nonconst_local_inds_host_view_type indices("indices",A_local->getLocalMaxNumRowEntries());
    nonconst_values_host_view_type values("values",A_local->getLocalMaxNumRowEntries());
    for(local_ordinal_type i = 0; i < numRows; i++) {
      size_t numEntries = 0;
      A_local->getLocalRowCopy(i, indices, values, numEntries);
      A_local_crs_nc->insertLocalValues(i, numEntries, reinterpret_cast<scalar_type*>(values.data()), indices.data());
    }
    A_local_crs_nc->fillComplete (A_local->getDomainMap (), A_local->getRangeMap ());
    A_local_crs = rcp_const_cast<const crs_matrix_type> (A_local_crs_nc);
  }

  return A_local_crs;
}



}

}

template<class MatrixType>
MDF<MatrixType>::MDF (const Teuchos::RCP<const row_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    Verbosity_(0),
    LevelOfFill_ (0),
    Overalloc_ (2.),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    initializeTime_ (0.0),
    computeTime_ (0.0),
    applyTime_ (0.0)
{
  allocateSolvers();
  allocatePermutations();
}


template<class MatrixType>
MDF<MatrixType>::MDF (const Teuchos::RCP<const crs_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    Verbosity_(0),
    LevelOfFill_ (0),
    Overalloc_ (2.),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    initializeTime_ (0.0),
    computeTime_ (0.0),
    applyTime_ (0.0)
{
  allocateSolvers();
  allocatePermutations();
}

template<class MatrixType>
void MDF<MatrixType>::allocatePermutations (bool force)
{
  if (A_.is_null()) return;

  // Allocate arrays as soon as size as known so their pointer is availabe
  if (force || permutations_.is_null() || A_->getLocalNumRows() != size_t(permutations_.size()))
  {
    permutations_ = Teuchos::null;
    reversePermutations_ = Teuchos::null;
    permutations_ = permutations_type(A_->getLocalNumRows());
    reversePermutations_ = permutations_type(A_->getLocalNumRows());
  }
}

template<class MatrixType>
void MDF<MatrixType>::allocateSolvers ()
{
  L_solver_ = Teuchos::null;
  U_solver_ = Teuchos::null;
  L_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> ());
  L_solver_->setObjectLabel("lower");
  U_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> ());
  U_solver_->setObjectLabel("upper");
}

template<class MatrixType>
void
MDF<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  if (A.getRawPtr () != A_.getRawPtr ()) {
    isAllocated_ = false;
    isInitialized_ = false;
    isComputed_ = false;
    A_local_ = Teuchos::null;
    MDF_handle_ = Teuchos::null;

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
    A_ = A;

    allocatePermutations(true);
  }
}



template<class MatrixType>
const typename MDF<MatrixType>::crs_matrix_type&
MDF<MatrixType>::getL () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    L_.is_null (), std::runtime_error, "Ifpack2::MDF::getL: The L factor "
    "is null.  Please call initialize() and compute() "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *L_;
}

template<class MatrixType>
typename MDF<MatrixType>::permutations_type &
MDF<MatrixType>::getPermutations () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    permutations_.is_null (), std::runtime_error, "Ifpack2::MDF::getPermutations: "
    "The permulations are null.  Please call initialize() and compute() "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return const_cast<permutations_type &>(permutations_);
}
template<class MatrixType>
typename MDF<MatrixType>::permutations_type &
MDF<MatrixType>::getReversePermutations () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    reversePermutations_.is_null (), std::runtime_error, "Ifpack2::MDF::getReversePermutations: "
    "The permulations are null.  Please call initialize() and compute() "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return const_cast<permutations_type &>(reversePermutations_);
}

template<class MatrixType>
const typename MDF<MatrixType>::crs_matrix_type&
MDF<MatrixType>::getU () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    U_.is_null (), std::runtime_error, "Ifpack2::MDF::getU: The U factor "
    "is null.  Please call initialize() and compute() "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *U_;
}


template<class MatrixType>
size_t MDF<MatrixType>::getNodeSmootherComplexity() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::MDF::getNodeSmootherComplexity: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  // MDF methods cost roughly one apply + the nnz in the upper+lower triangles
  if(!L_.is_null() && !U_.is_null())
    return A_->getLocalNumEntries() + L_->getLocalNumEntries() + U_->getLocalNumEntries();
  else
    return 0;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MDF<MatrixType>::local_ordinal_type,
                               typename MDF<MatrixType>::global_ordinal_type,
                               typename MDF<MatrixType>::node_type> >
MDF<MatrixType>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::MDF::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  // FIXME (mfh 25 Jan 2014) Shouldn't this just come from A_?
  TEUCHOS_TEST_FOR_EXCEPTION(
    L_.is_null (), std::runtime_error, "Ifpack2::MDF::getDomainMap: "
    "The computed graph is null.  Please call initialize() and compute() before calling "
    "this method.");
  return L_->getDomainMap ();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MDF<MatrixType>::local_ordinal_type,
                               typename MDF<MatrixType>::global_ordinal_type,
                               typename MDF<MatrixType>::node_type> >
MDF<MatrixType>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::MDF::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  // FIXME (mfh 25 Jan 2014) Shouldn't this just come from A_?
  TEUCHOS_TEST_FOR_EXCEPTION(
    L_.is_null (), std::runtime_error, "Ifpack2::MDF::getRangeMap: "
    "The computed graph is null.  Please call initialize() abd compute() before calling "
    "this method.");
  return L_->getRangeMap ();
}

template<class MatrixType>
void
MDF<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Details::getParamTryingTypes;
  const char prefix[] = "Ifpack2::MDF: ";

  // Default values of the various parameters.
  int fillLevel = 0;
  double overalloc = 2.;
  int verbosity = 0;

  // "fact: mdf level-of-fill" parsing is more complicated, because
  // we want to allow as many types as make sense.  int is the native
  // type, but we also want to accept double (for backwards
  // compatibilty with ILUT).  You can't cast arbitrary magnitude_type
  // (e.g., Sacado::MP::Vector) to int, so we use float instead, to
  // get coverage of the most common magnitude_type cases.  Weirdly,
  // there's an Ifpack2 test that sets the fill level as a
  // global_ordinal_type.
  {
    const std::string paramName ("fact: mdf level-of-fill");
    getParamTryingTypes<int, int, global_ordinal_type, double, float>
      (fillLevel, params, paramName, prefix);

    TEUCHOS_TEST_FOR_EXCEPTION
      (fillLevel != 0, std::runtime_error, prefix << "MDF with level of fill != 0 is not yet implemented.");
  }
  {
    const std::string paramName ("Verbosity");
    getParamTryingTypes<int, int, global_ordinal_type, double, float>
      (verbosity, params, paramName, prefix);
  }
  {
    const std::string paramName ("fact: mdf overalloc");
    getParamTryingTypes<double, double>
      (overalloc, params, paramName, prefix);
  }

  // Forward to trisolvers.
  L_solver_->setParameters(params);
  U_solver_->setParameters(params);

  // "Commit" the values only after validating all of them.  This
  // ensures that there are no side effects if this routine throws an
  // exception.

  LevelOfFill_ = fillLevel;
  Overalloc_ = overalloc;
  Verbosity_ = verbosity;
}


template<class MatrixType>
Teuchos::RCP<const typename MDF<MatrixType>::row_matrix_type>
MDF<MatrixType>::getMatrix () const {
  return Teuchos::rcp_implicit_cast<const row_matrix_type> (A_);
}


template<class MatrixType>
Teuchos::RCP<const typename MDF<MatrixType>::crs_matrix_type>
MDF<MatrixType>::getCrsMatrix () const {
  return Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_, true);
}


template<class MatrixType>
Teuchos::RCP<const typename MDF<MatrixType>::row_matrix_type>
MDF<MatrixType>::makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
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
  // should be the case if MDF is being used through
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
void MDF<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  const char prefix[] = "Ifpack2::MDF::initialize: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix is null.  Please "
     "call setMatrix() with a nonnull input before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_->isFillComplete (), std::runtime_error, prefix << "The matrix is not "
     "fill complete.  You may not invoke initialize() or compute() with this "
     "matrix until the matrix is fill complete.  If your matrix is a "
     "Tpetra::CrsMatrix, please call fillComplete on it (with the domain and "
     "range Maps, if appropriate) before calling this method.");

  Teuchos::Time timer ("MDF::initialize");
  double startTime = timer.wallTime();
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
    isAllocated_   = false;
    isComputed_    = false;
    MDF_handle_ = Teuchos::null;

    A_local_ = makeLocalFilter (A_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local_.is_null (), std::logic_error, "Ifpack2::MDF::initialize: "
      "makeLocalFilter returned null; it failed to compute A_local.  "
      "Please report this bug to the Ifpack2 developers.");

    // FIXME (mfh 24 Jan 2014, 26 Mar 2014) It would be more efficient
    // to rewrite MDF so that it works with any RowMatrix input, not
    // just CrsMatrix.  (That would require rewriting mdfGraph to
    // handle a Tpetra::RowGraph.)  However, to make it work for now,
    // we just copy the input matrix if it's not a CrsMatrix.
    {
      RCP<const crs_matrix_type> A_local_crs = Details::MDFImpl::get_local_crs_row_matrix(A_local_);

      auto A_local_device = A_local_crs->getLocalMatrixDevice();
      MDF_handle_ = rcp( new MDF_handle_device_type(A_local_device) );
      MDF_handle_->set_verbosity(Verbosity_);

      KokkosSparse::Experimental::mdf_symbolic(A_local_device,*MDF_handle_);

      isAllocated_ = true;
    }

    checkOrderingConsistency (*A_local_);
  } // Stop timing

  isInitialized_ = true;
  ++numInitialize_;
  initializeTime_ += (timer.wallTime() - startTime);
}

template<class MatrixType>
void
MDF<MatrixType>::
checkOrderingConsistency (const row_matrix_type& A)
{
  // First check that the local row map ordering is the same as the local portion of the column map.
  // The extraction of the strictly lower/upper parts of A, as well as the factorization,
  // implicitly assume that this is the case.
  Teuchos::ArrayView<const global_ordinal_type> rowGIDs = A.getRowMap()->getLocalElementList();
  Teuchos::ArrayView<const global_ordinal_type> colGIDs = A.getColMap()->getLocalElementList();
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
}

template<class MatrixType>
void MDF<MatrixType>::compute ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  const char prefix[] = "Ifpack2::MDF::compute: ";

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

  Teuchos::Time timer ("MDF::compute");

  // Start timing
  Teuchos::TimeMonitor timeMon (timer);
  double startTime = timer.wallTime();

  isComputed_ = false;

  {//Make sure values in A is picked up even in case of pattern reuse
    RCP<const crs_matrix_type> A_local_crs = Details::MDFImpl::get_local_crs_row_matrix(A_local_);

    // Compute the ordering and factorize
    auto A_local_device = A_local_crs->getLocalMatrixDevice();

    KokkosSparse::Experimental::mdf_numeric(A_local_device,*MDF_handle_);
  }

  // Ordering convention for MDF impl and here are reversed. Do reverse here to avoid confusion
  Details::MDFImpl::copy_dev_view_to_host_array(reversePermutations_, MDF_handle_->permutation);
  Details::MDFImpl::copy_dev_view_to_host_array(permutations_, MDF_handle_->permutation_inv);

  // TMR: Need to COPY the values held by the MDF handle because the CRS matrix needs to
  // exclusively own them and the MDF_handles use_count contribution throws that off
  {
    auto L_mdf = MDF_handle_->getL();
    L_ = rcp(new crs_matrix_type(
      A_local_->getRowMap (),
      A_local_->getColMap (),
      Details::MDFImpl::copy_view(L_mdf.graph.row_map),
      Details::MDFImpl::copy_view(L_mdf.graph.entries),
      Details::MDFImpl::copy_view(L_mdf.values)
    ));
  }
  {
    auto U_mdf = MDF_handle_->getU();
    U_ = rcp(new crs_matrix_type(
      A_local_->getRowMap (),
      A_local_->getColMap (),
      Details::MDFImpl::copy_view(U_mdf.graph.row_map),
      Details::MDFImpl::copy_view(U_mdf.graph.entries),
      Details::MDFImpl::copy_view(U_mdf.values)
    ));
  }
  L_->fillComplete ();
  U_->fillComplete ();
  L_solver_->setMatrix (L_);
  L_solver_->initialize ();
  L_solver_->compute ();
  U_solver_->setMatrix (U_);
  U_solver_->initialize ();
  U_solver_->compute ();

  isComputed_ = true;
  ++numCompute_;
  computeTime_ += (timer.wallTime() - startTime);
}

template<class MatrixType>
void
MDF<MatrixType>::
apply_impl (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  const scalar_type one = STS::one ();
  const scalar_type zero = STS::zero ();

  if (alpha == one && beta == zero) {
    MV tmp (Y.getMap (), Y.getNumVectors ());
    Details::MDFImpl::applyReorderingPermutations(X,tmp,permutations_);
    if (mode == Teuchos::NO_TRANS) { // Solve L (D (U Y)) = X for Y.
      // Start by solving L Y = X for Y.
      L_solver_->apply (tmp, Y, mode);
      U_solver_->apply (Y, tmp, mode); // Solve U Y = Y.
    }
    else { // Solve U^P (D^P (L^P Y)) = X for Y (where P is * or T).
      // Start by solving U^P Y = X for Y.
      U_solver_->apply (tmp, Y, mode);
      L_solver_->apply (Y, tmp, mode); // Solve L^P Y = Y.
    }
    Details::MDFImpl::applyReorderingPermutations(tmp,Y,reversePermutations_);
  }
  else { // alpha != 1 or beta != 0
    if (alpha == zero) {
      // The special case for beta == 0 ensures that if Y contains Inf
      // or NaN values, we replace them with 0 (following BLAS
      // convention), rather than multiplying them by 0 to get NaN.
      if (beta == zero) {
        Y.putScalar (zero);
      } else {
        Y.scale (beta);
      }
    } else { // alpha != zero
      MV Y_tmp (Y.getMap (), Y.getNumVectors ());
      apply_impl (X, Y_tmp, mode);
      Y.update (alpha, Y_tmp, beta);
    }
  }
}

template<class MatrixType>
void
MDF<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::MDF::apply: The matrix is "
    "null.  Please call setMatrix() with a nonnull input, then initialize() "
    "and compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::MDF::apply: If you have not yet called compute(), "
    "you must call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::MDF::apply: X and Y do not have the same number of columns.  "
    "X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isComplex && mode == Teuchos::CONJ_TRANS, std::logic_error,
    "Ifpack2::MDF::apply: mode = Teuchos::CONJ_TRANS is not implemented for "
    "complex Scalar type.  Please talk to the Ifpack2 developers to get this "
    "fixed.  There is a FIXME in this file about this very issue.");
#ifdef HAVE_IFPACK2_DEBUG
  {
    Teuchos::Array<magnitude_type> norms (X.getNumVectors ());
    X.norm1 (norms ());
    bool good = true;
    for (size_t j = 0; j < X.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION( ! good, std::runtime_error, "Ifpack2::MDF::apply: The 1-norm of the input X is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  Teuchos::Time timer ("MDF::apply");
  double startTime = timer.wallTime();
  { // Start timing
    Teuchos::TimeMonitor timeMon (timer);
    apply_impl(X,Y,mode,alpha,beta);
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
    TEUCHOS_TEST_FOR_EXCEPTION( ! good, std::runtime_error, "Ifpack2::MDF::apply: The 1-norm of the output Y is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  ++numApply_;
  applyTime_ += (timer.wallTime() - startTime);
}

template<class MatrixType>
std::string MDF<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::MDF\": {";
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

  if (! L_solver_.is_null ()) os << ", " << L_solver_->description ();
  if (! U_solver_.is_null ()) os << ", " << U_solver_->description ();

  os << "}";
  return os.str ();
}

} // namespace Ifpack2

#define IFPACK2_MDF_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::MDF< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif /* IFPACK2_MDF_DEF_HPP */
