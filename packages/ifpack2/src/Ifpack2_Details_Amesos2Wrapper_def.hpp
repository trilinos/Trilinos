// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP
#define IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP

#ifdef HAVE_IFPACK2_AMESOS2

#include "Ifpack2_LocalFilter.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// FIXME (mfh 25 Aug 2015) Work-around for Bug 6392.  This doesn't
// need to be a weak symbol as long as the Ifpack2 package depends on
// Amesos2.
namespace Amesos2 {
namespace Details {
  extern void registerLinearSolverFactory ();
} // namespace Details
} // namespace Amesos2

namespace Ifpack2 {
namespace Details {

template <class MatrixType>
Amesos2Wrapper<MatrixType>::
Amesos2Wrapper (const Teuchos::RCP<const row_matrix_type>& A) :
  A_(A),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false),
  SolverName_ ("")
{}

template <class MatrixType>
Amesos2Wrapper<MatrixType>::~Amesos2Wrapper()
{}

template <class MatrixType>
void Amesos2Wrapper<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // FIXME (mfh 12 Sep 2014) Why does this code make a deep copy of
  // the input ParameterList?  Does Amesos2 want a deep copy?

  // Extract the list called "Amesos2" that contains the Amesos2
  // solver's options.
  RCP<ParameterList> theList;
  if (params.name () == "Amesos2") {
    theList = rcp (new ParameterList (params));
  } else if (params.isSublist ("Amesos2")) {
    // FIXME (mfh 12 Sep 2014) This code actually makes _two_ deep copies.
    ParameterList subpl = params.sublist ("Amesos2");
    theList = rcp (new ParameterList (subpl));
    theList->setName ("Amesos2"); //FIXME hack until Teuchos sublist name bug is fixed
    if (params.isParameter ("Amesos2 solver name")) {
      SolverName_ = params.get<std::string>("Amesos2 solver name");
    }
  } else {
    // Amesos2 silently ignores any list not called "Amesos2".  We'll
    // throw an exception.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "The ParameterList passed to Amesos2 must be "
      "called \"Amesos2\".");
  }

  // If solver_ hasn't been allocated yet, cache the parameters and set them
  // once the concrete solver does exist.
  if (solver_.is_null ()) {
    parameterList_ = theList;
    return;
  }
  // FIXME (mfh 25 Aug 2015) Why doesn't this code set parameterList_
  // when the solver is NOT null?

  solver_->setParameters(theList);
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Amesos2Wrapper<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::row_matrix_type>
Amesos2Wrapper<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::map_type>
Amesos2Wrapper<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::map_type>
Amesos2Wrapper<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool Amesos2Wrapper<MatrixType>::hasTransposeApply () const {
  return true;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double Amesos2Wrapper<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template<class MatrixType>
double Amesos2Wrapper<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template<class MatrixType>
double Amesos2Wrapper<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}

template<class MatrixType>
void Amesos2Wrapper<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  IsInitialized_ = false;
  IsComputed_ = false;

  if (A.is_null ()) {
    A_ = Teuchos::null;
  }
  else {
    A_ = A;
  }

  // FIXME (mfh 10 Dec 2013) Currently, initialize() recreates
  // solver_ unconditionally, so this code won't have any
  // effect.  Once we fix initialize() so that it keeps
  // solver_, the code below will be effective.
  //if (! solver_.is_null ()) {
  //  solver_->setA (A_);
  //}
  // FIXME JJH 2014-July18 A_ might not be a locally filtered CRS matrix, which
  // means we have to do that dance all over again before calling solver_->setA ....
}

template<class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::row_matrix_type>
Amesos2Wrapper<MatrixType>::makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
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
void Amesos2Wrapper<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::Array;
  using Teuchos::ArrayView;

  const std::string timerName ("Ifpack2::Amesos2Wrapper::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    // Check that the matrix is nonnull.
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;

    RCP<const row_matrix_type> A_local = makeLocalFilter (A_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local.is_null (), std::logic_error, "Ifpack2::AmesosWrapper::initialize: "
      "makeLocalFilter returned null; it failed to compute A_local.  "
      "Please report this bug to the Ifpack2 developers.");

    {
      // The matrix that Amesos2 will build the preconditioner on must be a Tpetra::Crs matrix.
      // If A_local isn't, then we build one.
      A_local_crs_ = rcp_dynamic_cast<const crs_matrix_type> (A_local);

      if (A_local_crs_.is_null ()) {
        local_ordinal_type numRows = A_local->getLocalNumRows();
        Array<size_t> entriesPerRow(numRows);
        for(local_ordinal_type i = 0; i < numRows; i++)
        {
          entriesPerRow[i] = A_local->getNumEntriesInLocalRow(i);
        }
        RCP<crs_matrix_type> A_local_crs_nc =
          rcp (new crs_matrix_type (A_local->getRowMap (),
                                    A_local->getColMap (),
                                    entriesPerRow()));
        // copy entries into A_local_crs
	typename crs_matrix_type::nonconst_local_inds_host_view_type indices("Indices",A_local->getLocalMaxNumRowEntries() );
	typename crs_matrix_type::nonconst_values_host_view_type values("Values", A_local->getLocalMaxNumRowEntries());
        for(local_ordinal_type i = 0; i < numRows; i++)
        {
          size_t numEntries = 0;
          A_local->getLocalRowCopy(i, indices, values, numEntries);
          ArrayView<const local_ordinal_type> indicesInsert(indices.data(), numEntries);
          ArrayView<const scalar_type> valuesInsert((const scalar_type*)values.data(), numEntries);
          A_local_crs_nc->insertLocalValues(i, indicesInsert, valuesInsert);
        }
        A_local_crs_nc->fillComplete (A_local->getDomainMap (), A_local->getRangeMap ());
        A_local_crs_ = rcp_const_cast<const crs_matrix_type> (A_local_crs_nc);
      }
    }

    // FIXME (10 Dec 2013, 25 Aug 2015) It shouldn't be necessary to
    // recreate the solver each time, since
    // Trilinos::Details::LinearSolver has a setA() method.  See the
    // implementation of setMatrix().  I don't want to break anything
    // so I will leave the code as it is, possibly inefficient.


    // FIXME (mfh 25 Aug 2015) This is a work-around for Bug 6392.
    if (! Trilinos::Details::Impl::rememberRegisteredSomeLinearSolverFactory ("Amesos2")) {
      Amesos2::Details::registerLinearSolverFactory ();
    }

    solver_ = Trilinos::Details::getLinearSolver<MV, OP, typename MV::mag_type> ("Amesos2", SolverName_);
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver_.is_null (), std::runtime_error, "Ifpack2::Details::"
       "Amesos2Wrapper::initialize: Failed to create Amesos2 solver!");

    solver_->setMatrix (A_local_crs_);
    // If parameters have been already been cached via setParameters, set them now.
    if (parameterList_ != Teuchos::null) {
      setParameters (*parameterList_);
      parameterList_ = Teuchos::null;
    }
    // The symbolic factorization properly belongs to initialize(),
    // since initialize() is concerned with the matrix's structure
    // (and compute() with the matrix's values).
    solver_->symbolic ();
  } // Stop timing here.

  IsInitialized_ = true;
  ++NumInitialize_;

  InitializeTime_ += (timer->wallTime() - startTime);
}

template<class MatrixType>
void Amesos2Wrapper<MatrixType>::compute ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  const std::string timerName ("Ifpack2::Details::Amesos2Wrapper::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  { // Start timing here.
    TimeMonitor timeMon (*timer);
    solver_->numeric ();
  } // Stop timing here.

  IsComputed_ = true;
  ++NumCompute_;

  ComputeTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
void Amesos2Wrapper<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // X = RHS
  // Y = solution

  const std::string timerName ("Ifpack2::Amesos2Wrapper::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::Amesos2Wrapper::apply: You must call compute() to compute the "
      "incomplete factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::Amesos2Wrapper::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      mode != Teuchos::NO_TRANS, std::logic_error,
      "Ifpack2::Amesos2Wrapper::apply: Solving with the transpose (mode == "
      "Teuchos::TRANS) or conjugate transpose (Teuchos::CONJ_TRANS) is not "
      "implemented.");

    // If alpha != 1 or beta != 0, create a temporary multivector
    // Y_temp to hold the contents of alpha*M^{-1}*X.  Otherwise,
    // alias Y_temp to Y.
    RCP<MV> Y_temp = (alpha != STS::one () || beta != STS::zero ()) ?
      rcp (new MV (Y.getMap (), Y.getNumVectors ())) :
      rcpFromRef (Y);

    // If X and Y are pointing to the same memory location, create an
    // auxiliary vector, X_temp, so that we don't clobber the input
    // when computing the output.  Otherwise, alias X_temp to X.
    RCP<const MV> X_temp;
    {
      if (X.aliases(Y)) {
        X_temp = rcp (new MV (X, Teuchos::Copy));
      } else {
        X_temp = rcpFromRef (X);
      }
    }

    // Set up "local" views of X and Y.
    RCP<const MV> X_local;
    RCP<MV> Y_local;
    //JJH 15-Apr-2016 I changed this from ">=" to ">".  Otherwise the else block
    //is never hit.
    //bmk 6-19-17: previously, the next line only set multipleProcs if A_ was distributed
    //  This doesn't work if A_ is local but X/Y are distributed, as in AdditiveSchwarz.
    const bool multipleProcs = (A_->getRowMap ()->getComm ()->getSize () > 1) || (X.getMap ()->getComm ()->getSize () > 1);
    if (multipleProcs) {
      // Interpret X and Y as "local" multivectors, that is, in the
      // local filter's domain resp. range Maps.  "Interpret" means that
      // we create views with different Maps; we don't have to copy.
      X_local = X_temp->offsetView (A_local_crs_->getDomainMap (), 0);
      Y_local = Y_temp->offsetViewNonConst (A_local_crs_->getRangeMap (), 0);
    }
    else { // only one process in A_'s communicator
      // X and Y are already "local"; no need to set up local views.
      X_local = X_temp;
      Y_local = Y_temp;
    }

    // Use the precomputed factorization to solve.
    solver_->solve (*Y_local, *X_local);

    if (alpha != STS::one () || beta != STS::zero ()) {
      Y.update (alpha, *Y_temp, beta);
    }
  } // Stop timing here.

  ++NumApply_;

  ApplyTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
std::string Amesos2Wrapper<MatrixType>::description () const {
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Amesos2Wrapper\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false")
     << ", Computed: " << (isComputed () ? "true" : "false");

  if (A_local_crs_.is_null ()) {
    os << ", Matrix: null";
  }
  else {
    os << ", Global matrix dimensions: ["
       << A_local_crs_->getGlobalNumRows () << ", " << A_local_crs_->getGlobalNumCols () << "]";
  }

  // If the actual solver happens to implement Describable, have it
  // describe itself.  Otherwise, don't print anything.
  if (! solver_.is_null ()) {
    Teuchos::Describable* d = dynamic_cast<Teuchos::Describable*> (solver_.getRawPtr ());
    if (d != NULL) {
      os << ", {";
      os << d->description ();
      os << "}";
    }
  }
  os << "}";
  return os.str ();
}


template <class MatrixType>
void
Amesos2Wrapper<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ?
    Teuchos::VERB_LOW : verbLevel;

  // describe() starts, by convention, with a tab before it prints anything.
  OSTab tab0 (out);
  if (vl > Teuchos::VERB_NONE) {
    out << "\"Ifpack2::Amesos2Wrapper\":" << endl;
    OSTab tab1 (out);
    out << "MatrixType: \"" << TypeNameTraits<MatrixType>::name ()
        << "\"" << endl;

    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
    }

    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl;
    out << "Computed: " << (isComputed () ? "true" : "false") << endl;
    out << "Number of initialize calls: " << getNumInitialize () << endl;
    out << "Number of compute calls: " << getNumCompute () << endl;
    out << "Number of apply calls: " << getNumApply () << endl;
    out << "Total time in seconds for initialize: " << getInitializeTime () << endl;
    out << "Total time in seconds for compute: " << getComputeTime () << endl;
    out << "Total time in seconds for apply: " << getApplyTime () << endl;

    if (vl > Teuchos::VERB_LOW) {
      out << "Matrix:" << endl;
      A_local_crs_->describe (out, vl);
    }
  }
}

} // namespace Details
} // namespace Ifpack2

// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_DETAILS_AMESOS2WRAPPER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Details::Amesos2Wrapper< Tpetra::RowMatrix<S, LO, GO, N> >;

#else

#define IFPACK2_DETAILS_AMESOS2WRAPPER_INSTANT(S,LO,GO,N)

#endif // HAVE_IFPACK2_AMESOS2
#endif // IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP
