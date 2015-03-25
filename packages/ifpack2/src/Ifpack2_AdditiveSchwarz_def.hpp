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

#ifndef IFPACK2_ADDITIVESCHWARZ_DEF_HPP
#define IFPACK2_ADDITIVESCHWARZ_DEF_HPP

#include "Ifpack2_AdditiveSchwarz_decl.hpp"

// AdditiveSchwarz uses OneLevelFactory to create a default inner
// preconditioner.
//
// FIXME (mfh 13 Dec 2013) For some inexplicable reason, I have to
// include the _decl and _def headers separately here; including just
// Ifpack2_Details_OneLevelFactory.hpp doesn't work.  It probably has
// something to do with ETI, but I don't fully understand what.
#include "Ifpack2_Details_OneLevelFactory_decl.hpp"
#include "Ifpack2_Details_OneLevelFactory_def.hpp"

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
#include "Xpetra_RowMatrix.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Zoltan2_XpetraRowMatrixAdapter.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#include "Zoltan2_OrderingSolution.hpp"
#endif

#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#include "Ifpack2_SupportGraph_decl.hpp"
#endif

#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_LocalFilter_def.hpp"
#include "Ifpack2_OverlappingRowMatrix_def.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_ReorderFilter_def.hpp"
#include "Ifpack2_SingletonFilter_def.hpp"

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include <locale> // std::toupper

namespace Ifpack2 {

namespace Details {

/// \class OneLevelPreconditionerNamer
/// \brief Map from an Ifpack2::Preconditioner subclass to its string name.
/// \tparam PrecType Specialization of a subclass of Ifpack2::Preconditioner.
///
/// \warning This class is an implementation detail of
///   Ifpack2::AdditiveSchwarz.  Its interface may change or it may go
///   away at any time.
///
/// Ifpack2::AdditiveSchwarz uses this class to map from its
/// compile-time template parameter \c LocalInverseType, to a string
/// name of the inner preconditioner which it can give to
/// Details::OneLevelFactory in order to create the inner
/// preconditioner.  This class will no longer be needed once
/// Ifpack2::AdditiveSchwarz no longer has a LocalInverseType template
/// parameter.
template<class PrecType>
class OneLevelPreconditionerNamer {
public:
  //! Name corresponding to Preconditioner subclass PrecType.
  static std::string name () {
    // The default implementation returns an invalid preconditioner
    // name.  This ensures that AdditiveSchwarz won't try to create a
    // preconditioner it doesn't know how to create.  This is better
    // than providing a valid default that is a different class than
    // the user expects.
    return "INVALID";
  }
};

//
// Partial specialization for Ifpack2::Preconditioner.
// It picks a reasonable default subdomain solver.
//

template<class S, class LO, class GO, class NT>
class OneLevelPreconditionerNamer< ::Ifpack2::Preconditioner<S, LO, GO, NT> > {
public:
  static std::string name () {
    // The default inner preconditioner is "ILUT", for backwards
    // compatibility with the original AdditiveSchwarz implementation.
    return "ILUT";
  }
};

//
// Partial specializations for each single-level preconditioner.
//

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Chebyshev<MatrixType> > {
public:
  static std::string name () {
    return "CHEBYSHEV";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Details::DenseSolver<MatrixType> > {
public:
  static std::string name () {
    return "DENSE";
  }
};

#ifdef HAVE_IFPACK2_AMESOS2
template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Details::Amesos2Wrapper<MatrixType> > {
public:
  static std::string name () {
    return "AMESOS2";
  }
};
#endif // HAVE_IFPACK2_AMESOS2

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Diagonal<MatrixType> > {
public:
  static std::string name () {
    return "DIAGONAL";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::ILUT<MatrixType> > {
public:
  static std::string name () {
    return "ILUT";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Relaxation<MatrixType> > {
public:
  static std::string name () {
    return "RELAXATION";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::RILUK<MatrixType> > {
public:
  static std::string name () {
    return "RILUK";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Krylov<MatrixType> > {
public:
  static std::string name () {
    return "KRYLOV";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::IdentitySolver<MatrixType> > {
public:
  static std::string name () {
    return "IDENTITY";
  }
};

} // namespace Details


template<class MatrixType, class LocalInverseType>
bool
AdditiveSchwarz<MatrixType, LocalInverseType>::hasInnerPrecName () const
{
  const char* options[4] = {
    "inner preconditioner name",
    "subdomain solver name",
    "schwarz: inner preconditioner name",
    "schwarz: subdomain solver name"
  };
  const int numOptions = 4;
  bool match = false;
  for (int k = 0; k < numOptions && ! match; ++k) {
    if (List_.isParameter (options[k])) {
      return true;
    }
  }
  return false;
}


template<class MatrixType, class LocalInverseType>
void
AdditiveSchwarz<MatrixType, LocalInverseType>::removeInnerPrecName ()
{
  const char* options[4] = {
    "inner preconditioner name",
    "subdomain solver name",
    "schwarz: inner preconditioner name",
    "schwarz: subdomain solver name"
  };
  const int numOptions = 4;
  for (int k = 0; k < numOptions; ++k) {
    List_.remove (options[k], false);
  }
}





template<class MatrixType, class LocalInverseType>
std::string
AdditiveSchwarz<MatrixType, LocalInverseType>::innerPrecName () const
{
  const char* options[4] = {
    "inner preconditioner name",
    "subdomain solver name",
    "schwarz: inner preconditioner name",
    "schwarz: subdomain solver name"
  };
  const int numOptions = 4;
  std::string newName;
  bool match = false;

  // As soon as one parameter option matches, ignore all others.
  for (int k = 0; k < numOptions && ! match; ++k) {
    if (List_.isParameter (options[k])) {
      // try-catch block protects against incorrect type errors.
      //
      // FIXME (mfh 04 Jan 2013) We should instead catch and report
      // type errors.
      try {
        newName = List_.get<std::string> (options[k]);
        match = true;
      } catch (...) {}
    }
  }
  return match ? newName : defaultInnerPrecName ();
}


template<class MatrixType, class LocalInverseType>
void
AdditiveSchwarz<MatrixType, LocalInverseType>::removeInnerPrecParams ()
{
  const char* options[4] = {
    "inner preconditioner parameters",
    "subdomain solver parameters",
    "schwarz: inner preconditioner parameters",
    "schwarz: subdomain solver parameters"
  };
  const int numOptions = 4;

  // As soon as one parameter option matches, ignore all others.
  for (int k = 0; k < numOptions; ++k) {
    List_.remove (options[k], false);
  }
}


template<class MatrixType, class LocalInverseType>
std::pair<Teuchos::ParameterList, bool>
AdditiveSchwarz<MatrixType, LocalInverseType>::innerPrecParams () const
{
  const char* options[4] = {
    "inner preconditioner parameters",
    "subdomain solver parameters",
    "schwarz: inner preconditioner parameters",
    "schwarz: subdomain solver parameters"
  };
  const int numOptions = 4;
  Teuchos::ParameterList params;

  // As soon as one parameter option matches, ignore all others.
  bool match = false;
  for (int k = 0; k < numOptions && ! match; ++k) {
    if (List_.isSublist (options[k])) {
      params = List_.sublist (options[k]);
      match = true;
    }
  }
  // Default is an empty list of parameters.
  return std::make_pair (params, match);
}


template<class MatrixType, class LocalInverseType>
std::string
AdditiveSchwarz<MatrixType, LocalInverseType>::defaultInnerPrecName ()
{
  // FIXME (mfh 14 Dec 2013) We want to get rid of the
  // LocalInverseType template parameter.  Soon, we will add an "inner
  // preconditioner" string parameter to the input ParameterList.  For
  // now, we map statically from LocalInverseType to its string name,
  // and use the string name to create the inner preconditioner.
  return Details::OneLevelPreconditionerNamer<LocalInverseType>::name ();
}


template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A) :
  Matrix_ (A),
  IsInitialized_ (false),
  IsComputed_ (false),
  IsOverlapping_ (false),
  OverlapLevel_ (0),
  CombineMode_ (Tpetra::ZERO),
  UseReordering_ (false),
  ReorderingAlgorithm_ ("none"),
  UseSubdomain_ (false),
  FilterSingletons_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0)
{
  Teuchos::ParameterList plist;
  setParameters (plist); // Set parameters to default values
}

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                 const int overlapLevel) :
  Matrix_ (A),
  IsInitialized_ (false),
  IsComputed_ (false),
  IsOverlapping_ (false),
  OverlapLevel_ (overlapLevel),
  CombineMode_ (Tpetra::ZERO),
  UseReordering_ (false),
  ReorderingAlgorithm_ ("none"),
  UseSubdomain_ (false),
  FilterSingletons_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0)
{
  Teuchos::ParameterList plist;
  setParameters (plist); // Set parameters to default values
}


template<class MatrixType,class LocalInverseType>
AdditiveSchwarz<MatrixType,LocalInverseType>::~AdditiveSchwarz () {}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type > >
AdditiveSchwarz<MatrixType,LocalInverseType>::getDomainMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
    "getDomainMap: The matrix to precondition is null.  You must either pass "
    "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
    "input, before you may call this method.");
  return Matrix_->getDomainMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
AdditiveSchwarz<MatrixType,LocalInverseType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
    "getRangeMap: The matrix to precondition is null.  You must either pass "
    "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
    "input, before you may call this method.");
  return Matrix_->getRangeMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > AdditiveSchwarz<MatrixType,LocalInverseType>::getMatrix() const
{
  return Matrix_;
}


template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::string timerName ("Ifpack2::AdditiveSchwarz::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    const scalar_type ZERO = Teuchos::ScalarTraits<scalar_type>::zero ();
    const scalar_type ONE = Teuchos::ScalarTraits<scalar_type>::one ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! IsComputed_, std::runtime_error,
      "Ifpack2::AdditiveSchwarz::apply: "
      "isComputed() must be true before you may call apply().");
    TEUCHOS_TEST_FOR_EXCEPTION(
      Matrix_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::apply: "
      "The input matrix A is null, but the preconditioner says that it has "
      "been computed (isComputed() is true).  This should never happen, since "
      "setMatrix() should always mark the preconditioner as not computed if "
      "its argument is null.  "
      "Please report this bug to the Ifpack2 developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverse_.is_null (), std::runtime_error,
      "Ifpack2::AdditiveSchwarz::apply: The subdomain solver is null.  "
      "This can only happen if you called setInnerPreconditioner() with a null "
      "input, after calling initialize() or compute().  If you choose to call "
      "setInnerPreconditioner() with a null input, you must then call it with "
      "a nonnull input before you may call initialize() or compute().");
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
      "Ifpack2::AdditiveSchwarz::apply: "
      "X and Y must have the same number of columns.  X has "
      << X.getNumVectors() << " columns, but Y has " << Y.getNumVectors() << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      alpha != ONE, std::logic_error,
      "Ifpack2::AdditiveSchwarz::apply: Not implemented for alpha != 1.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      beta != ZERO, std::logic_error,
      "Ifpack2::AdditiveSchwarz::apply: Not implemented for beta != 0.");

    const size_t numVectors = X.getNumVectors ();

    RCP<MV> OverlappingX,OverlappingY;

    if (IsOverlapping_) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        OverlappingMatrix_.is_null (), std::logic_error,
        "Ifpack2::AdditiveSchwarz::apply: The overlapping matrix is null.  "
        "This should never happen if IsOverlapping_ is true.  "
        "Please report this bug to the Ifpack2 developers.");

      // Setup if we're overlapping
      //
      // MV's constructor fills with zeros.
      OverlappingX = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
      OverlappingY = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
      OverlappingMatrix_->importMultiVector (X, *OverlappingX, Tpetra::INSERT);
      // FIXME from Ifpack1: Will not work with non-zero starting solutions.
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        localMap_.is_null (), std::logic_error,
        "Ifpack2::AdditiveSchwarz::apply: localMap_ is null.");

      // MV's constructor fills with zeros.
      //
      // localMap_ has the same number of indices on each process that
      // Matrix_->getRowMap() does on that process.  Thus, we can do
      // the Import step without creating a new MV, just by viewing
      // OverlappingX using Matrix_->getRowMap ().
      OverlappingX = rcp (new MV (localMap_, numVectors));
      OverlappingY = rcp (new MV (localMap_, numVectors));

      RCP<MV> globalOverlappingX =
        OverlappingX->offsetViewNonConst (Matrix_->getRowMap (), 0);

      // Create Import object on demand, if necessary.
      if (DistributedImporter_.is_null ()) {
        // FIXME (mfh 15 Apr 2014) Why can't we just ask the Matrix
        // for its Import object?  Of course a general RowMatrix might
        // not necessarily have one.
        DistributedImporter_ =
          rcp (new import_type (Matrix_->getRowMap (),
                                Matrix_->getDomainMap ()));
      }
      globalOverlappingX->doImport (X, *DistributedImporter_, Tpetra::INSERT);
    }

    if (FilterSingletons_) {
      // process singleton filter
      MV ReducedX (SingletonMatrix_->getRowMap (), numVectors);
      MV ReducedY (SingletonMatrix_->getRowMap (), numVectors);
      SingletonMatrix_->SolveSingletons (*OverlappingX, *OverlappingY);
      SingletonMatrix_->CreateReducedRHS (*OverlappingY, *OverlappingX, ReducedX);

      // process reordering
      if (! UseReordering_) {
        Inverse_->apply (ReducedX, ReducedY);
      }
      else {
        MV ReorderedX (ReducedX, Teuchos::Copy);
        MV ReorderedY (ReducedY, Teuchos::Copy);
        ReorderedLocalizedMatrix_->permuteOriginalToReordered (ReducedX, ReorderedX);
        Inverse_->apply (ReorderedX, ReorderedY);
        ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, ReducedY);
      }

      // finish up with singletons
      SingletonMatrix_->UpdateLHS (ReducedY, *OverlappingY);
    }
    else {

      // process reordering
      if (! UseReordering_) {
        Inverse_->apply (*OverlappingX, *OverlappingY);
      }
      else {
        MV ReorderedX (*OverlappingX, Teuchos::Copy);
        MV ReorderedY (*OverlappingY, Teuchos::Copy);
        ReorderedLocalizedMatrix_->permuteOriginalToReordered (*OverlappingX, ReorderedX);
        Inverse_->apply (ReorderedX, ReorderedY);
        ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, *OverlappingY);
      }
    }

    if (IsOverlapping_) {
      OverlappingMatrix_->exportMultiVector (*OverlappingY, Y, CombineMode_);
    }
    else {
      // mfh 16 Apr 2014: Make a view of Y with the same Map as
      // OverlappingY, so that we can copy OverlappingY into Y.  This
      // replaces code that iterates over all entries of OverlappingY,
      // copying them one at a time into Y.  That code assumed that
      // the rows of Y and the rows of OverlappingY have the same
      // global indices in the same order; see Bug 5992.
      RCP<MV> Y_view = Y.offsetViewNonConst (OverlappingY->getMap (), 0);
      Tpetra::deep_copy (*Y_view, *OverlappingY);
    }
  } // Stop timing here.

  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::
setParameters (const Teuchos::ParameterList& plist)
{
  // mfh 18 Nov 2013: Ifpack2's setParameters() method passes in the
  // input list as const.  This means that we have to copy it before
  // validation or passing into setParameterList().
  List_ = plist;
  this->setParameterList (Teuchos::rcpFromRef (List_));
}



template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::
setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  using Tpetra::CombineMode;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntry;
  using Teuchos::ParameterEntryValidator;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::StringToIntegralParameterEntryValidator;

  if (plist.is_null ()) {
    // Assume that the user meant to set default parameters by passing
    // in an empty list.
    this->setParameterList (Teuchos::parameterList ());
  }
  // At this point, plist should be nonnull.
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameterList: plist is null.  This should never happen, since the "
    "method should have replaced a null input list with a nonnull empty list "
    "by this point.  Please report this bug to the Ifpack2 developers.");

  // try {
  //   List_.validateParameters (* getValidParameters ());
  // }
  // catch (std::exception& e) {
  //   std::cerr << "Ifpack2::AdditiveSchwarz::setParameterList: Validation failed with the following error message: " << e.what () << std::endl;
  //   throw e;
  // }

  // mfh 18 Nov 2013: Supplying the current value as the default value
  // when calling ParameterList::get() ensures "delta" behavior when
  // users pass in new parameters: any unspecified parameters in the
  // new list retain their values in the old list.  This preserves
  // backwards compatiblity with this class' previous behavior.  Note
  // that validateParametersAndSetDefaults() would have different
  // behavior: any parameters not in the new list would get default
  // values, which could be different than their values in the
  // original list.

  bool gotCombineMode = false;
  try {
    CombineMode_ = getIntegralValue<Tpetra::CombineMode> (List_, "schwarz: combine mode");
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterName&) {
    // The caller didn't provide that parameter.  Just keep the
    // existing value of CombineMode_.
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterType&) {
    // The user perhaps supplied it as an Tpetra::CombineMode enum
    // value.  Let's try again (below).  If it doesn't succeed, we
    // know that the type is wrong, so we can let it throw whatever
    // exception it would throw.
  }
  // Try to get the combine mode as an integer.
  if (! gotCombineMode) {
    try {
      CombineMode_ = plist->get ("schwarz: combine mode", CombineMode_);
      gotCombineMode = true;
    }
    catch (Teuchos::Exceptions::InvalidParameterType&) {}
  }
  // Try to get the combine mode as a string.  If this works, use the
  // validator to convert to int.  This is painful, but necessary in
  // order to do validation, since the input list doesn't come with a
  // validator.
  if (! gotCombineMode) {
    const ParameterEntry& validEntry =
      getValidParameters ()->getEntry ("schwarz: combine mode");
    RCP<const ParameterEntryValidator> v = validEntry.validator ();
    typedef StringToIntegralParameterEntryValidator<CombineMode> vs2e_type;
    RCP<const vs2e_type> vs2e = rcp_dynamic_cast<const vs2e_type> (v, true);

    const ParameterEntry& inputEntry = plist->getEntry ("schwarz: combine mode");
    CombineMode_ = vs2e->getIntegralValue (inputEntry, "schwarz: combine mode");
    gotCombineMode = true;
  }
  (void) gotCombineMode; // forestall "set but not used" compiler warning

  OverlapLevel_ = plist->get ("schwarz: overlap level", OverlapLevel_);

  // We set IsOverlapping_ in initialize(), once we know that Matrix_ is nonnull.

  // Will we be doing reordering?  Unlike Ifpack, we'll use a
  // "schwarz: reordering list" to give to Zoltan2.
  UseReordering_ = plist->get ("schwarz: use reordering", UseReordering_);

#if !defined(HAVE_IFPACK2_XPETRA) || !defined(HAVE_IFPACK2_ZOLTAN2)
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseReordering_, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified \"schwarz: use reordering\" = true.  "
    "This is only valid when Trilinos was built with Ifpack2, Xpetra, and "
    "Zoltan2 enabled.  Either Xpetra or Zoltan2 was not enabled in your build "
    "of Trilinos.");
#endif

  // FIXME (mfh 18 Nov 2013) Now would be a good time to validate the
  // "schwarz: reordering list" parameter list.  Currently, that list
  // gets extracted in setup().

  // Subdomain check
  if (plist->isParameter ("schwarz: subdomain id") &&
      plist->get ("schwarz: subdomain id", -1) > 0) {
    UseSubdomain_ = true;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseSubdomain_, std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified the \"schwarz: subdomain id\" parameter, "
    "with a value other than -1.  This parameter is not yet supported.");

  // if true, filter singletons. NOTE: the filtered matrix can still have
  // singletons! A simple example: upper triangular matrix, if I remove
  // the lower node, I still get a matrix with a singleton! However, filter
  // singletons should help for PDE problems with Dirichlet BCs.
  FilterSingletons_ = plist->get ("schwarz: filter singletons", FilterSingletons_);

  // If the inner solver doesn't exist yet, don't create it.
  // initialize() creates it.
  //
  // If the inner solver _does_ exist, there are three cases,
  // depending on what the user put in the input ParameterList.
  //
  //   1. The user did /not/ provide a parameter specifying the inner
  //      solver's type, nor did the user specify a sublist of
  //      parameters for the inner solver
  //   2. The user did /not/ provide a parameter specifying the inner
  //      solver's type, but /did/ specify a sublist of parameters for
  //      the inner solver
  //   3. The user provided a parameter specifying the inner solver's
  //      type (it does not matter in this case whether the user gave
  //      a sublist of parameters for the inner solver)
  //
  // AdditiveSchwarz has "delta" (relative) semantics for setting
  // parameters.  This means that if the user did not specify the
  // inner solver's type, we presume that the type has not changed.
  // Thus, if the inner solver exists, we don't need to recreate it.
  //
  // In Case 3, if the user bothered to specify the inner solver's
  // type, then we must assume it may differ than the current inner
  // solver's type.  Thus, we have to recreate the inner solver.  We
  // achieve this here by assigning null to Inverse_; initialize()
  // will recreate the solver when it is needed.  Our assumption here
  // is necessary because Ifpack2::Preconditioner does not have a
  // method for querying a preconditioner's "type" (i.e., name) as a
  // string.  Remember that the user may have previously set an
  // arbitrary inner solver by calling setInnerPreconditioner().
  //
  // See note at the end of setInnerPreconditioner().

  if (! Inverse_.is_null ()) {
    // "CUSTOM" explicitly indicates that the user called or plans to
    // call setInnerPreconditioner.
    if (hasInnerPrecName () && innerPrecName () != "CUSTOM") {
      // Wipe out the current inner solver.  initialize() will
      // recreate it with the correct type.
      Inverse_ = Teuchos::null;
    }
    else {
      // Extract and apply the sublist of parameters to give to the
      // inner solver, if there is such a sublist of parameters.
      std::pair<Teuchos::ParameterList, bool> result = innerPrecParams ();
      if (result.second) {
        Inverse_->setParameters (result.first);
      }
    }
  }
}



template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Teuchos::ParameterList>
AdditiveSchwarz<MatrixType,LocalInverseType>::
getValidParameters () const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

  if (validParams_.is_null ()) {
    const int overlapLevel = 0;
    const bool useReordering = false;
    const bool filterSingletons = false;
    ParameterList reorderingSublist;
    reorderingSublist.set ("order_method", std::string ("rcm"));

    RCP<ParameterList> plist = parameterList ("Ifpack2::AdditiveSchwarz");

    Tpetra::setCombineModeParameter (*plist, "schwarz: combine mode");
    plist->set ("schwarz: overlap level", overlapLevel);
    plist->set ("schwarz: use reordering", useReordering);
    plist->set ("schwarz: reordering list", reorderingSublist);
    // mfh 24 Mar 2015: We accept this for backwards compatibility
    // ONLY.  It is IGNORED.
    plist->set ("schwarz: compute condest", false);
    plist->set ("schwarz: filter singletons", filterSingletons);

    // FIXME (mfh 18 Nov 2013) Get valid parameters from inner solver.
    //
    // FIXME (mfh 18 Nov 2013) Get valid parameters from Zoltan2, if
    // Zoltan2 was enabled in the build.

    validParams_ = rcp_const_cast<const ParameterList> (plist);
  }
  return validParams_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::initialize ()
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  const std::string timerName ("Ifpack2::AdditiveSchwarz::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
      "initialize: The matrix to precondition is null.  You must either pass "
      "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
      "input, before you may call this method.");

    IsInitialized_ = false;
    IsComputed_ = false;

    RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
    RCP<const map_type> rowMap = Matrix_->getRowMap ();
    RCP<node_type> node = Matrix_->getNode ();
    const global_size_t INVALID =
      Teuchos::OrdinalTraits<global_size_t>::invalid ();

    // If there's only one process in the matrix's communicator,
    // then there's no need to compute overlap.
    if (comm->getSize () == 1) {
      OverlapLevel_ = 0;
      IsOverlapping_ = false;
    } else if (OverlapLevel_ != 0) {
      IsOverlapping_ = true;
    }

    if (OverlapLevel_ == 0) {
      const global_ordinal_type indexBase = rowMap->getIndexBase ();
      RCP<const SerialComm<int> > localComm (new SerialComm<int> ());
      // FIXME (mfh 15 Apr 2014) What if indexBase isn't the least
      // global index in the list of GIDs on this process?
      localMap_ =
        rcp (new map_type (INVALID, rowMap->getNodeNumElements (),
                           indexBase, localComm, node));
    }

    // compute the overlapping matrix if necessary
    if (IsOverlapping_) {
      if (UseSubdomain_) {
        const int sid = List_.get ("subdomain id", -1);
        OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_, sid));
      } else {
        OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_));
      }
    }

    setup (); // This does a lot of the initialization work.

    if (! Inverse_.is_null ()) {
      Inverse_->initialize (); // Initialize subdomain solver.
    }

  } // Stop timing here.

  IsInitialized_ = true;
  ++NumInitialize_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  InitializeTime_ = timer->totalElapsedTime ();
}


template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isInitialized () const
{
  return IsInitialized_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::compute ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  if (! IsInitialized_) {
    initialize ();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isInitialized (), std::logic_error, "Ifpack2::AdditiveSchwarz::compute: "
    "The preconditioner is not yet initialized, "
    "even though initialize() supposedly has been called.  "
    "This should never happen.  "
    "Please report this bug to the Ifpack2 developers.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_.is_null (), std::runtime_error,
    "Ifpack2::AdditiveSchwarz::compute: The subdomain solver is null.  "
    "This can only happen if you called setInnerPreconditioner() with a null "
    "input, after calling initialize() or compute().  If you choose to call "
    "setInnerPreconditioner() with a null input, you must then call it with a "
    "nonnull input before you may call initialize() or compute().");

  const std::string timerName ("Ifpack2::AdditiveSchwarz::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    IsComputed_ = false;
    Inverse_->compute ();
  } // Stop timing here.

  IsComputed_ = true;
  ++NumCompute_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
}

//==============================================================================
// Returns true if the  preconditioner has been successfully computed, false otherwise.
template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isComputed() const
{
  return IsComputed_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumInitialize() const
{
  return NumInitialize_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumCompute() const
{
  return NumCompute_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumApply() const
{
  return NumApply_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getInitializeTime() const
{
  return InitializeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getComputeTime() const
{
  return ComputeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getApplyTime() const
{
  return ApplyTime_;
}


template<class MatrixType,class LocalInverseType>
std::string AdditiveSchwarz<MatrixType,LocalInverseType>::description () const
{
  std::ostringstream out;

  out << "\"Ifpack2::AdditiveSchwarz\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  out << "Initialized: " << (isInitialized () ? "true" : "false")
      << ", Computed: " << (isComputed () ? "true" : "false")
      << ", Overlap level: " << OverlapLevel_
      << ", Subdomain reordering: \"" << ReorderingAlgorithm_ << "\"";
  out << ", Combine mode: \"";
  if (CombineMode_ == Tpetra::INSERT) {
    out << "INSERT";
  } else if (CombineMode_ == Tpetra::ADD) {
    out << "ADD";
  } else if (CombineMode_ == Tpetra::REPLACE) {
    out << "REPLACE";
  } else if (CombineMode_ == Tpetra::ABSMAX) {
    out << "ABSMAX";
  } else if (CombineMode_ == Tpetra::ZERO) {
    out << "ZERO";
  }
  out << "\"";
  if (Matrix_.is_null ()) {
    out << ", Matrix: null";
  }
  else {
    out << ", Global matrix dimensions: ["
        << Matrix_->getGlobalNumRows () << ", "
        << Matrix_->getGlobalNumCols () << "]";
  }
  out << ", Inner solver: ";
  if (!Inverse_.is_null ())
    out << "{" << Inverse_->description() << "}";
  else
    out << "null";

  out << "}";
  return out.str ();
}


template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using std::endl;

  const int myRank = Matrix_->getComm ()->getRank ();
  const int numProcs = Matrix_->getComm ()->getSize ();
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl > Teuchos::VERB_NONE) {
    // describe() starts with a tab, by convention.
    OSTab tab0 (out);
    if (myRank == 0) {
      out << "\"Ifpack2::AdditiveSchwarz\":";
    }
    OSTab tab1 (out);
    if (myRank == 0) {
      out << "MatrixType: " << TypeNameTraits<MatrixType>::name () << endl;
      out << "LocalInverseType: " << TypeNameTraits<LocalInverseType>::name () << endl;
      if (this->getObjectLabel () != "") {
        out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
      }

      out << "Overlap level: " << OverlapLevel_ << endl
          << "Combine mode: \"";
      if (CombineMode_ == Tpetra::INSERT) {
        out << "INSERT";
      } else if (CombineMode_ == Tpetra::ADD) {
        out << "ADD";
      } else if (CombineMode_ == Tpetra::REPLACE) {
        out << "REPLACE";
      } else if (CombineMode_ == Tpetra::ABSMAX) {
        out << "ABSMAX";
      } else if (CombineMode_ == Tpetra::ZERO) {
        out << "ZERO";
      }
      out << "\"" << endl
          << "Subdomain reordering: \"" << ReorderingAlgorithm_ << "\"" << endl;
    }

    if (Matrix_.is_null ()) {
      if (myRank == 0) {
        out << "Matrix: null" << endl;
      }
    }
    else {
      if (myRank == 0) {
        out << "Matrix:" << endl;
        std::flush (out);
      }
      Matrix_->getComm ()->barrier (); // wait for output to finish
      Matrix_->describe (out, Teuchos::VERB_LOW);
    }

    if (myRank == 0) {
      out << "Number of initialize calls: " << getNumInitialize () << endl
          << "Number of compute calls: " << getNumCompute () << endl
          << "Number of apply calls: " << getNumApply () << endl
          << "Total time in seconds for initialize: " << getInitializeTime () << endl
          << "Total time in seconds for compute: " << getComputeTime () << endl
          << "Total time in seconds for apply: " << getApplyTime () << endl;
    }

    if (Inverse_.is_null ()) {
      if (myRank == 0) {
        out << "Subdomain solver: null" << endl;
      }
    }
    else {
      if (vl < Teuchos::VERB_EXTREME) {
        if (myRank == 0) {
          out << "Subdomain solver: not null" << endl;
        }
      }
      else { // vl >= Teuchos::VERB_EXTREME
        for (int p = 0; p < numProcs; ++p) {
          if (p == myRank) {
            out << "Subdomain solver on Process " << myRank << ":";
            if (Inverse_.is_null ()) {
              out << "null" << endl;
            } else {
              out << endl;
              Inverse_->describe (out, vl);
            }
          }
          Matrix_->getComm ()->barrier ();
          Matrix_->getComm ()->barrier ();
          Matrix_->getComm ()->barrier (); // wait for output to finish
        }
      }
    }

    Matrix_->getComm ()->barrier (); // wait for output to finish
  }
}


template<class MatrixType,class LocalInverseType>
std::ostream& AdditiveSchwarz<MatrixType,LocalInverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getOverlapLevel() const
{
  return OverlapLevel_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::setup ()
{
#ifdef HAVE_MPI
  using Teuchos::MpiComm;
#endif // HAVE_MPI
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  typedef Xpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraMatrixType;
  typedef Xpetra::TpetraRowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraTpetraMatrixType;
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
    "initialize: The matrix to precondition is null.  You must either pass "
    "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
    "input, before you may call this method.");

  // Localized version of Matrix_ or OverlappingMatrix_.
  RCP<row_matrix_type> LocalizedMatrix;

  // The "most current local matrix."  At the end of this method, this
  // will be handed off to the inner solver.
  RCP<row_matrix_type> ActiveMatrix;

  // Create localized matrix.
  if (! OverlappingMatrix_.is_null ()) {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
      //
      // FIXME (mfh 18 Nov 2013) btw what's the difference between
      // Ifpack_NodeFilter and Ifpack_LocalFilter?  The former's
      // documentation sounds a lot like what Ifpack2::LocalFilter
      // does.
      //
      //Ifpack2_NodeFilter *tt = new Ifpack2_NodeFilter(OverlappingMatrix_,nodeID); //FIXME
      //LocalizedMatrix = Teuchos::rcp(tt);
    }
    else
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (OverlappingMatrix_));
  }
  else {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
    }
    else {
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (Matrix_));
    }
  }

  // Sanity check; I don't trust the logic above to have created LocalizedMatrix.
  TEUCHOS_TEST_FOR_EXCEPTION(
    LocalizedMatrix.is_null (), std::logic_error,
    "Ifpack2::AdditiveSchwarz::setup: LocalizedMatrix is null, after the code "
    "that claimed to have created it.  This should never be the case.  Please "
    "report this bug to the Ifpack2 developers.");

  // Mark localized matrix as active
  ActiveMatrix = LocalizedMatrix;

  // Singleton Filtering
  if (FilterSingletons_) {
    SingletonMatrix_ = rcp (new SingletonFilter<row_matrix_type> (LocalizedMatrix));
    ActiveMatrix = SingletonMatrix_;
  }

  // Do reordering
  if (UseReordering_) {
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    // Unlike Ifpack, Zoltan2 does all the dirty work here.
    Teuchos::ParameterList zlist = List_.sublist ("schwarz: reordering list");

    // FIXME (mfh 18 Nov 2013) Shouldn't this come from the Zoltan2 sublist?
    ReorderingAlgorithm_ = List_.get<std::string> ("order_method", "rcm");
    XpetraTpetraMatrixType XpetraMatrix (ActiveMatrix);
    typedef Zoltan2::XpetraRowMatrixAdapter<XpetraMatrixType> z2_adapter_type;
    z2_adapter_type Zoltan2Matrix (rcpFromRef (XpetraMatrix));
    typedef Zoltan2::OrderingProblem<z2_adapter_type> ordering_problem_type;
#ifdef HAVE_MPI
    // Grab the MPI Communicator and build the ordering problem with that
    MPI_Comm myRawComm;

    RCP<const MpiComm<int> > mpicomm =
      rcp_dynamic_cast<const MpiComm<int> > (ActiveMatrix->getComm ());
    if (mpicomm == Teuchos::null) {
      myRawComm = MPI_COMM_SELF;
    } else {
      myRawComm = * (mpicomm->getRawMpiComm ());
    }
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist, myRawComm);
#else
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist);
#endif
    MyOrderingProblem.solve ();

    // Now create the reordered matrix & mark it as active
    {
      typedef ReorderFilter<row_matrix_type> reorder_filter_type;
      typedef Zoltan2::OrderingSolution<global_ordinal_type,
        local_ordinal_type> ordering_solution_type;

      ordering_solution_type sol (*MyOrderingProblem.getSolution ());

      // perm[i] gives the where OLD index i shows up in the NEW
      // ordering.  revperm[i] gives the where NEW index i shows
      // up in the OLD ordering.  Note that perm is actually the
      // "inverse permutation," in Zoltan2 terms.
      ArrayRCP<local_ordinal_type> perm = sol.getPermutationRCPConst (true);
      ArrayRCP<local_ordinal_type> revperm = sol.getPermutationRCPConst ();

      ReorderedLocalizedMatrix_ =
        rcp (new reorder_filter_type (ActiveMatrix, perm, revperm));

      ActiveMatrix = ReorderedLocalizedMatrix_;
    }
#else
    // This is a logic_error, not a runtime_error, because
    // setParameters() should have excluded this case already.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Ifpack2::AdditiveSchwarz::setup: "
      "The Zoltan2 and Xpetra packages must be enabled in order "
      "to support reordering.");
#endif
  }

  innerMatrix_ = ActiveMatrix;

  TEUCHOS_TEST_FOR_EXCEPTION(
    innerMatrix_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setup: Inner matrix is null right before constructing inner solver.  "
    "Please report this bug to the Ifpack2 developers.");

  // Construct the inner solver if necessary.
  if (Inverse_.is_null ()) {
    const std::string innerName = innerPrecName ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerName == "INVALID", std::logic_error,
      "Ifpack2::AdditiveSchwarz::initialize: AdditiveSchwarz doesn't "
      "know how to create an instance of your LocalInverseType \""
      << Teuchos::TypeNameTraits<LocalInverseType>::name () << "\".  If "
      "LocalInverseType is a single-level preconditioner (does not take an "
      "inner solver), then you can fix this in one of two ways.  Either (a) "
      "create the LocalInverseType instance yourself and give it to "
      "AdditiveSchwarz by calling setInnerPreconditioner(), before calling "
      "initialize(), or (b) teach Details::OneLevelFactory how to create an "
      "inner preconditioner of that type.  "
      "Please talk to the Ifpack2 developers for details.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      innerName == "CUSTOM", std::runtime_error, "Ifpack2::AdditiveSchwarz::"
      "initialize: If the \"inner preconditioner name\" parameter (or any "
      "alias thereof) has the value \"CUSTOM\", then you must first call "
      "setInnerPreconditioner with a nonnull inner preconditioner input before "
      "you may call initialize().");

    Details::OneLevelFactory<MatrixType> factory;
    RCP<prec_type> innerPrec = factory.create (innerName, innerMatrix_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerPrec.is_null (), std::logic_error,
      "Ifpack2::AdditiveSchwarz::setup: Failed to create inner preconditioner "
      "with name \"" << innerName << "\".");

    // Extract and apply the sublist of parameters to give to the
    // inner solver, if there is such a sublist of parameters.
    std::pair<Teuchos::ParameterList, bool> result = innerPrecParams ();
    if (result.second) {
      innerPrec->setParameters (result.first);
    }
    Inverse_ = innerPrec; // "Commit" the inner solver.
  }
  else if (Inverse_->getMatrix ().getRawPtr () != innerMatrix_.getRawPtr ()) {
    // The new inner matrix is different from the inner
    // preconditioner's current matrix, so give the inner
    // preconditioner the new inner matrix.  First make sure that the
    // inner solver knows how to have its matrix changed.
    typedef Details::CanChangeMatrix<row_matrix_type> can_change_type;
    can_change_type* innerSolver =
      dynamic_cast<can_change_type*> (Inverse_.getRawPtr ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerSolver == NULL, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
      "setup: The current inner preconditioner does not implement the "
      "setMatrix() feature.  Only preconditioners that inherit from "
      "Ifpack2::Details::CanChangeMatrix implement this feature.");

    // Give the new inner matrix to the inner preconditioner.
    innerSolver->setMatrix (innerMatrix_);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setup: Inverse_ is null right after we were supposed to have created it."
    "  Please report this bug to the Ifpack2 developers.");

  // We don't have to call setInnerPreconditioner() here, because we
  // had the inner matrix (innerMatrix_) before creation of the inner
  // preconditioner.  Calling setInnerPreconditioner here would be
  // legal, but it would require an unnecessary reset of the inner
  // preconditioner (i.e., calling initialize() and compute() again).
}


template<class MatrixType, class LocalInverseType>
void AdditiveSchwarz<MatrixType, LocalInverseType>::
setInnerPreconditioner (const Teuchos::RCP<Preconditioner<scalar_type,
                                                          local_ordinal_type,
                                                          global_ordinal_type,
                                                          node_type> >& innerPrec)
{
  if (! innerPrec.is_null ()) {
    // Make sure that the new inner solver knows how to have its matrix changed.
    typedef Details::CanChangeMatrix<row_matrix_type> can_change_type;
    can_change_type* innerSolver = dynamic_cast<can_change_type*> (&*innerPrec);
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerSolver == NULL, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
      "setInnerPreconditioner: The input preconditioner does not implement the "
      "setMatrix() feature.  Only input preconditioners that inherit from "
      "Ifpack2::Details::CanChangeMatrix implement this feature.");

    // If users provide an inner solver, we assume that
    // AdditiveSchwarz's current inner solver parameters no longer
    // apply.  (In fact, we will remove those parameters from
    // AdditiveSchwarz's current list below.)  Thus, we do /not/ apply
    // the current sublist of inner solver parameters to the input
    // inner solver.

    // mfh 03 Jan 2014: Thanks to Paul Tsuji for pointing out that it's
    // perfectly legal for innerMatrix_ to be null here.  This can
    // happen if initialize() has not been called yet.  For example,
    // when Factory creates an AdditiveSchwarz instance, it calls
    // setInnerPreconditioner() without first calling initialize().

    // Give the local matrix to the new inner solver.
    innerSolver->setMatrix (innerMatrix_);

    // If the user previously specified a parameter for the inner
    // preconditioner's type, then clear out that parameter and its
    // associated sublist.  Replace the inner preconditioner's type with
    // "CUSTOM", to make it obvious that AdditiveSchwarz's ParameterList
    // does not necessarily describe the current inner preconditioner.
    // We have to remove all allowed aliases of "inner preconditioner
    // name" before we may set it to "CUSTOM".  Users may also set this
    // parameter to "CUSTOM" themselves, but this is not required.
    removeInnerPrecName ();
    removeInnerPrecParams ();
    List_.set ("inner preconditioner name", "CUSTOM");

    // Bring the new inner solver's current status (initialized or
    // computed) in line with AdditiveSchwarz's current status.
    if (isInitialized ()) {
      innerPrec->initialize ();
    }
    if (isComputed ()) {
      innerPrec->compute ();
    }
  }

  // If the new inner solver is null, we don't change the initialized
  // or computed status of AdditiveSchwarz.  That way, AdditiveSchwarz
  // won't have to recompute innerMatrix_ if the inner solver changes.
  // This does introduce a new error condition in compute() and
  // apply(), but that's OK.

  // Set the new inner solver.
  Inverse_ = innerPrec;
}

template<class MatrixType, class LocalInverseType>
void AdditiveSchwarz<MatrixType, LocalInverseType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Don't set the matrix unless it is different from the current one.
  if (A.getRawPtr () != Matrix_.getRawPtr ()) {
    IsInitialized_ = false;
    IsComputed_ = false;

    // Reset all the state computed in initialize() and compute().
    OverlappingMatrix_ = Teuchos::null;
    ReorderedLocalizedMatrix_ = Teuchos::null;
    innerMatrix_ = Teuchos::null;
    SingletonMatrix_ = Teuchos::null;
    localMap_ = Teuchos::null;
    DistributedImporter_ = Teuchos::null;

    Matrix_ = A;
  }
}

} // namespace Ifpack2

// FIXME (mfh 16 Sep 2014) We should really only use RowMatrix here!
// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.
#define IFPACK2_ADDITIVESCHWARZ_INSTANT(S,LO,GO,N) \
  template class Ifpack2::AdditiveSchwarz< Tpetra::RowMatrix<S, LO, GO, N> >; \
  template class Ifpack2::AdditiveSchwarz< Tpetra::CrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::AdditiveSchwarz< Tpetra::RowMatrix<S, LO, GO, N>, \
                                           Ifpack2::ILUT<Tpetra::RowMatrix< S, LO, GO, N > > >; \
  template class Ifpack2::AdditiveSchwarz< Tpetra::CrsMatrix<S, LO, GO, N>, \
                                           Ifpack2::ILUT<Tpetra::CrsMatrix< S, LO, GO, N > > >;


#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
