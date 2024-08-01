// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_AdditiveSchwarz_def.hpp
/// \brief Definition of Ifpack2::AdditiveSchwarz, which implements
///   additive Schwarz preconditioning with an arbitrary subdomain
///   solver.  For the declaration and class documentation, please see
///   Ifpack2_AdditiveSchwarz_decl.hpp in this directory.
///
/// If you want to use Ifpack2::AdditiveSchwarz directly in your
/// application, please include the automatically generated header
/// file Ifpack2_AdditiveSchwarz.hpp.  In general, you should never
/// need to include Ifpack2_AdditiveSchwarz_def.hpp in your
/// application directly.

#ifndef IFPACK2_ADDITIVESCHWARZ_DEF_HPP
#define IFPACK2_ADDITIVESCHWARZ_DEF_HPP

#include "Ifpack2_AdditiveSchwarz_decl.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
// We need Ifpack2's implementation of LinearSolver, because we use it
// to wrap the user-provided Ifpack2::Preconditioner in
// Ifpack2::AdditiveSchwarz::setInnerPreconditioner.
#include "Ifpack2_Details_LinearSolver.hpp"
#include "Ifpack2_Details_getParamTryingTypes.hpp"

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
#include "Zoltan2_TpetraRowGraphAdapter.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#include "Zoltan2_OrderingSolution.hpp"
#endif

#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_ReorderFilter.hpp"
#include "Ifpack2_SingletonFilter.hpp"
#include "Ifpack2_Details_AdditiveSchwarzFilter.hpp"

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <locale> // std::toupper

#include <Tpetra_BlockMultiVector.hpp>


// FIXME (mfh 25 Aug 2015) Work-around for Bug 6392.  This doesn't
// need to be a weak symbol because it only refers to a function in
// the Ifpack2 package.
namespace Ifpack2 {
namespace Details {
  extern void registerLinearSolverFactory ();
} // namespace Details
} // namespace Ifpack2

#ifdef HAVE_IFPACK2_DEBUG

namespace { // (anonymous)

  template<class MV>
  bool
  anyBad (const MV& X)
  {
    using STS = Teuchos::ScalarTraits<typename MV::scalar_type>;
    using magnitude_type = typename STS::magnitudeType;
    using STM = Teuchos::ScalarTraits<magnitude_type>;

    Teuchos::Array<magnitude_type> norms (X.getNumVectors ());
    X.norm2 (norms ());
    bool good = true;
    for (size_t j = 0; j < X.getNumVectors (); ++j) {
      if (STM::isnaninf (norms[j])) {
        good = false;
        break;
      }
    }
    return ! good;
  }

} // namespace (anonymous)

#endif // HAVE_IFPACK2_DEBUG

namespace Ifpack2 {

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
    const Teuchos::ParameterEntry* paramEnt =
      List_.getEntryPtr (options[k]);
    if (paramEnt != nullptr && paramEnt->isType<std::string> ()) {
      newName = Teuchos::getValue<std::string> (*paramEnt);
      match = true;
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
  // The default inner preconditioner is "ILUT", for backwards
  // compatibility with the original AdditiveSchwarz implementation.
  return "ILUT";
}

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A) :
  Matrix_ (A)
{}

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                 const int overlapLevel) :
  Matrix_ (A),
  OverlapLevel_ (overlapLevel)
{}

template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type > >
AdditiveSchwarz<MatrixType,LocalInverseType>::
getDomainMap () const
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


namespace
{

template<class MatrixType, class map_type>
Teuchos::RCP<const map_type>
pointMapFromMeshMap(const Teuchos::RCP<const map_type> & meshMap,const typename MatrixType::local_ordinal_type blockSize)
{
  using BMV = Tpetra::BlockMultiVector<
      typename MatrixType::scalar_type,
      typename MatrixType::local_ordinal_type,
      typename MatrixType::global_ordinal_type,
      typename MatrixType::node_type>;

  if (blockSize == 1) return meshMap;
  
  return Teuchos::RCP<const map_type>(new map_type(BMV::makePointMap (*meshMap,blockSize)));
}

template <typename MV, typename Map>
void resetMultiVecIfNeeded(std::unique_ptr<MV> &mv_ptr, const Map &map, const size_t numVectors, bool initialize)
{
  if(!mv_ptr || mv_ptr->getNumVectors() != numVectors) {
    mv_ptr.reset(new MV(map, numVectors, initialize));
  }
}

} // namespace

template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &B,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  const char prefix[] = "Ifpack2::AdditiveSchwarz::apply: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (! IsComputed_, std::runtime_error,
     prefix << "isComputed() must be true before you may call apply().");
  TEUCHOS_TEST_FOR_EXCEPTION
    (Matrix_.is_null (), std::logic_error, prefix <<
     "The input matrix A is null, but the preconditioner says that it has "
     "been computed (isComputed() is true).  This should never happen, since "
     "setMatrix() should always mark the preconditioner as not computed if "
     "its argument is null.  "
     "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (Inverse_.is_null (), std::runtime_error,
     prefix << "The subdomain solver is null.  "
     "This can only happen if you called setInnerPreconditioner() with a null "
     "input, after calling initialize() or compute().  If you choose to call "
     "setInnerPreconditioner() with a null input, you must then call it with "
     "a nonnull input before you may call initialize() or compute().");
  TEUCHOS_TEST_FOR_EXCEPTION
    (B.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
     prefix << "B and Y must have the same number of columns.  B has " <<
     B.getNumVectors () << " columns, but Y has " << Y.getNumVectors() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (IsOverlapping_ && OverlappingMatrix_.is_null (), std::logic_error,
     prefix << "The overlapping matrix is null.  "
     "This should never happen if IsOverlapping_ is true.  "
     "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! IsOverlapping_ && localMap_.is_null (), std::logic_error,
     prefix << "localMap_ is null.  "
     "This should never happen if IsOverlapping_ is false.  "
     "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (alpha != STS::one (), std::logic_error,
     prefix << "Not implemented for alpha != 1.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (beta != STS::zero (), std::logic_error,
     prefix << "Not implemented for beta != 0.");

#ifdef HAVE_IFPACK2_DEBUG
  {
    const bool bad = anyBad (B);
    TEUCHOS_TEST_FOR_EXCEPTION
      (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
       "The 2-norm of the input B is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

#ifdef HAVE_IFPACK2_DEBUG
  if (! ZeroStartingSolution_) {
    const bool bad = anyBad (Y);
    TEUCHOS_TEST_FOR_EXCEPTION
      (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
       "On input, the initial guess Y has 2-norm NaN or Inf "
       "(ZeroStartingSolution_ is false).");
  }
#endif // HAVE_IFPACK2_DEBUG

  const std::string timerName ("Ifpack2::AdditiveSchwarz::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }
  double startTime = timer->wallTime();

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    const scalar_type ZERO = Teuchos::ScalarTraits<scalar_type>::zero ();
    const size_t numVectors = B.getNumVectors ();

    // mfh 25 Apr 2015: Fix for currently failing
    // Ifpack2_AdditiveSchwarz_RILUK test.
    if (ZeroStartingSolution_) {
      Y.putScalar (ZERO);
    }

    // set up for overlap communication
    MV* OverlappingB = nullptr;
    MV* OverlappingY = nullptr;
    {
      RCP<const map_type> B_and_Y_map = pointMapFromMeshMap<MatrixType>(IsOverlapping_ ?
        OverlappingMatrix_->getRowMap () : localMap_ , Matrix_->getBlockSize());
      resetMultiVecIfNeeded(overlapping_B_, B_and_Y_map, numVectors, false);
      resetMultiVecIfNeeded(overlapping_Y_, B_and_Y_map, numVectors, false);
      OverlappingB = overlapping_B_.get ();
      OverlappingY = overlapping_Y_.get ();
      // FIXME (mfh 25 Jun 2019) It's not clear whether we really need
      // to fill with zeros here, but that's what was happening before.
      OverlappingB->putScalar (ZERO);
      OverlappingY->putScalar (ZERO);
    }

    RCP<MV> globalOverlappingB;
    if (! IsOverlapping_) {
      auto matrixPointRowMap = pointMapFromMeshMap<MatrixType>(Matrix_->getRowMap (),Matrix_->getBlockSize ());

      globalOverlappingB =
        OverlappingB->offsetViewNonConst (matrixPointRowMap, 0);

      // Create Import object on demand, if necessary.
      if (DistributedImporter_.is_null ()) {
        // FIXME (mfh 15 Apr 2014) Why can't we just ask the Matrix
        // for its Import object?  Of course a general RowMatrix might
        // not necessarily have one.
        DistributedImporter_ =
          rcp (new import_type (matrixPointRowMap,
                                Matrix_->getDomainMap ()));
      }
    }

    resetMultiVecIfNeeded(R_, B.getMap(), numVectors, false);
    resetMultiVecIfNeeded(C_, Y.getMap(), numVectors, false);
    // If taking averages in overlap region, we need to compute
    // the number of procs who have a copy of each overlap dof
    Teuchos::ArrayRCP<scalar_type>  dataNumOverlapCopies;
    if (IsOverlapping_ && AvgOverlap_) {
      if (num_overlap_copies_.get()  == nullptr) {
        num_overlap_copies_.reset (new MV (Y.getMap (), 1, false));
        RCP<MV> onesVec( new MV(OverlappingMatrix_->getRowMap(), 1, false) );
        onesVec->putScalar(Teuchos::ScalarTraits<scalar_type>::one());
        rcp_dynamic_cast<OverlappingRowMatrix<row_matrix_type>> (OverlappingMatrix_)->exportMultiVector (*onesVec, *(num_overlap_copies_.get ()), CombineMode_);
      }
      dataNumOverlapCopies = num_overlap_copies_.get ()->getDataNonConst(0);
    }

    MV* R = R_.get ();
    MV* C = C_.get ();

    // FIXME (mfh 25 Jun 2019) It was never clear whether C had to be
    // initialized to zero.  R definitely should not need this.
    C->putScalar (ZERO);

    for (int ni=0; ni<NumIterations_; ++ni)
    {
#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (Y);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At top of iteration " << ni << ", the 2-norm of Y is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

      Tpetra::deep_copy(*R, B);

      // if (ZeroStartingSolution_ && ni == 0) {
      //   Y.putScalar (STS::zero ());
      // }
      if (!ZeroStartingSolution_ || ni > 0) {
        //calculate residual
        Matrix_->apply (Y, *R, mode, -STS::one(), STS::one());

#ifdef HAVE_IFPACK2_DEBUG
        {
          const bool bad = anyBad (*R);
          TEUCHOS_TEST_FOR_EXCEPTION
            (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
             "At iteration " << ni << ", the 2-norm of R (result of computing "
             "residual with Y) is NaN or Inf.");
        }
#endif // HAVE_IFPACK2_DEBUG
      }

      // do communication if necessary
      if (IsOverlapping_) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (OverlappingMatrix_.is_null (), std::logic_error, prefix <<
           "IsOverlapping_ is true, but OverlappingMatrix_, while nonnull, is "
           "not an OverlappingRowMatrix<row_matrix_type>.  Please report this "
           "bug to the Ifpack2 developers.");
        OverlappingMatrix_->importMultiVector (*R, *OverlappingB, Tpetra::INSERT);

        //JJH We don't need to import the solution Y we are always solving AY=R with initial guess zero
        //if (ZeroStartingSolution_ == false)
        //  overlapMatrix->importMultiVector (Y, *OverlappingY, Tpetra::INSERT);
        /*
          FIXME from Ifpack1: Will not work with non-zero starting solutions.
          TODO  JJH 3/20/15   I don't know whether this comment is still valid.

          Here is the log for the associated commit 720b2fa4 to Ifpack1:

          "Added a note to recall that the nonzero starting solution will not
          work properly if reordering, filtering or wider overlaps are used. This only
          applied to methods like Jacobi, Gauss-Seidel, and SGS (in both point and block
          version), and not to ILU-type preconditioners."
        */

#ifdef HAVE_IFPACK2_DEBUG
        {
          const bool bad = anyBad (*OverlappingB);
          TEUCHOS_TEST_FOR_EXCEPTION
            (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
             "At iteration " << ni << ", result of importMultiVector from R "
             "to OverlappingB, has 2-norm NaN or Inf.");
        }
#endif // HAVE_IFPACK2_DEBUG
      } else {
        globalOverlappingB->doImport (*R, *DistributedImporter_, Tpetra::INSERT);

#ifdef HAVE_IFPACK2_DEBUG
        {
          const bool bad = anyBad (*globalOverlappingB);
          TEUCHOS_TEST_FOR_EXCEPTION
            (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
             "At iteration " << ni << ", result of doImport from R, has 2-norm "
             "NaN or Inf.");
        }
#endif // HAVE_IFPACK2_DEBUG
      }

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (*OverlappingB);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At iteration " << ni << ", right before localApply, the 2-norm of "
           "OverlappingB is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

      // local solve
      localApply(*OverlappingB, *OverlappingY);

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (*OverlappingY);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At iteration " << ni << ", after localApply and before export / "
           "copy, the 2-norm of OverlappingY is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (*C);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At iteration " << ni << ", before export / copy, the 2-norm of C "
           "is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

      // do communication if necessary
      if (IsOverlapping_) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (OverlappingMatrix_.is_null (), std::logic_error, prefix
           << "OverlappingMatrix_ is null when it shouldn't be.  "
           "Please report this bug to the Ifpack2 developers.");
        OverlappingMatrix_->exportMultiVector (*OverlappingY, *C, CombineMode_);

        // average solution in overlap regions if requested via "schwarz: combine mode" "AVG"
        if (AvgOverlap_) {
          Teuchos::ArrayRCP<scalar_type>  dataC = C->getDataNonConst(0);
          for (int i = 0; i < (int) C->getMap()->getLocalNumElements(); i++) {
            dataC[i] = dataC[i]/dataNumOverlapCopies[i];
          }
        }
      }
      else {
        // mfh 16 Apr 2014: Make a view of Y with the same Map as
        // OverlappingY, so that we can copy OverlappingY into Y.  This
        // replaces code that iterates over all entries of OverlappingY,
        // copying them one at a time into Y.  That code assumed that
        // the rows of Y and the rows of OverlappingY have the same
        // global indices in the same order; see Bug 5992.
        RCP<MV> C_view = C->offsetViewNonConst (OverlappingY->getMap (), 0);
        Tpetra::deep_copy (*C_view, *OverlappingY);
      }

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (*C);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At iteration " << ni << ", before Y := C + Y, the 2-norm of C "
           "is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (Y);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "Before Y := C + Y, at iteration " << ni << ", the 2-norm of Y "
           "is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG

      Y.update(UpdateDamping_, *C, STS::one());

#ifdef HAVE_IFPACK2_DEBUG
      {
        const bool bad = anyBad (Y);
        TEUCHOS_TEST_FOR_EXCEPTION
          (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
           "At iteration " << ni << ", after Y := C + Y, the 2-norm of Y "
           "is NaN or Inf.");
      }
#endif // HAVE_IFPACK2_DEBUG
    } // for each iteration

  } // Stop timing here

#ifdef HAVE_IFPACK2_DEBUG
  {
    const bool bad = anyBad (Y);
    TEUCHOS_TEST_FOR_EXCEPTION
      (bad, std::runtime_error, "Ifpack2::AdditiveSchwarz::apply: "
       "The 2-norm of the output Y is NaN or Inf.");
  }
#endif // HAVE_IFPACK2_DEBUG

  ++NumApply_;

  ApplyTime_ += (timer->wallTime() - startTime);
}

template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
localApply (MV& OverlappingB, MV& OverlappingY) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  const size_t numVectors = OverlappingB.getNumVectors ();

  auto additiveSchwarzFilter = rcp_dynamic_cast<Details::AdditiveSchwarzFilter<MatrixType>>(innerMatrix_);
  if(additiveSchwarzFilter)
  {
    //Create the reduced system innerMatrix_ * ReducedY = ReducedB.
    //This effectively fuses 3 tasks:
    //  -SingletonFilter::SolveSingletons (solve entries of OverlappingY corresponding to singletons)
    //  -SingletonFilter::CreateReducedRHS (fill ReducedReorderedB from OverlappingB, with entries in singleton columns eliminated)
    //  -ReorderFilter::permuteOriginalToReordered (apply permutation to ReducedReorderedB)
    resetMultiVecIfNeeded(reduced_reordered_B_, additiveSchwarzFilter->getRowMap(), numVectors, true);
    resetMultiVecIfNeeded(reduced_reordered_Y_, additiveSchwarzFilter->getRowMap(), numVectors, true);
    additiveSchwarzFilter->CreateReducedProblem(OverlappingB, OverlappingY, *reduced_reordered_B_);
    //Apply inner solver
    Inverse_->solve (*reduced_reordered_Y_, *reduced_reordered_B_);
    //Scatter ReducedY back to non-singleton rows of OverlappingY, according to the reordering.
    additiveSchwarzFilter->UpdateLHS(*reduced_reordered_Y_, OverlappingY);
  }
  else
  {
    if (FilterSingletons_) {
      // process singleton filter
      resetMultiVecIfNeeded(reduced_B_, SingletonMatrix_->getRowMap(), numVectors, true);
      resetMultiVecIfNeeded(reduced_Y_, SingletonMatrix_->getRowMap(), numVectors, true);

      RCP<SingletonFilter<row_matrix_type> > singletonFilter =
        rcp_dynamic_cast<SingletonFilter<row_matrix_type> > (SingletonMatrix_);
      TEUCHOS_TEST_FOR_EXCEPTION
        (! SingletonMatrix_.is_null () && singletonFilter.is_null (),
         std::logic_error, "Ifpack2::AdditiveSchwarz::localApply: "
         "SingletonFilter_ is nonnull but is not a SingletonFilter"
         "<row_matrix_type>.  This should never happen.  Please report this bug "
         "to the Ifpack2 developers.");
      singletonFilter->SolveSingletons (OverlappingB, OverlappingY);
      singletonFilter->CreateReducedRHS (OverlappingY, OverlappingB, *reduced_B_);

      // process reordering
      if (! UseReordering_) {
        Inverse_->solve (*reduced_Y_, *reduced_B_);
      }
      else {
        RCP<ReorderFilter<row_matrix_type> > rf =
          rcp_dynamic_cast<ReorderFilter<row_matrix_type> > (ReorderedLocalizedMatrix_);
        TEUCHOS_TEST_FOR_EXCEPTION
          (! ReorderedLocalizedMatrix_.is_null () && rf.is_null (), std::logic_error,
           "Ifpack2::AdditiveSchwarz::localApply: ReorderedLocalizedMatrix_ is "
           "nonnull but is not a ReorderFilter<row_matrix_type>.  This should "
           "never happen.  Please report this bug to the Ifpack2 developers.");
        resetMultiVecIfNeeded(reordered_B_, reduced_B_->getMap(), numVectors, false);
        resetMultiVecIfNeeded(reordered_Y_, reduced_Y_->getMap(), numVectors, false);
        rf->permuteOriginalToReordered (*reduced_B_, *reordered_B_);
        Inverse_->solve (*reordered_Y_, *reordered_B_);
        rf->permuteReorderedToOriginal (*reordered_Y_, *reduced_Y_);
      }

      // finish up with singletons
      singletonFilter->UpdateLHS (*reduced_Y_, OverlappingY);
    }
    else {
      // process reordering
      if (! UseReordering_) {
        Inverse_->solve (OverlappingY, OverlappingB);
      }
      else {
        resetMultiVecIfNeeded(reordered_B_, OverlappingB.getMap(), numVectors, false);
        resetMultiVecIfNeeded(reordered_Y_, OverlappingY.getMap(), numVectors, false);

        RCP<ReorderFilter<row_matrix_type> > rf =
          rcp_dynamic_cast<ReorderFilter<row_matrix_type> > (ReorderedLocalizedMatrix_);
        TEUCHOS_TEST_FOR_EXCEPTION
          (! ReorderedLocalizedMatrix_.is_null () && rf.is_null (), std::logic_error,
           "Ifpack2::AdditiveSchwarz::localApply: ReorderedLocalizedMatrix_ is "
           "nonnull but is not a ReorderFilter<row_matrix_type>.  This should "
           "never happen.  Please report this bug to the Ifpack2 developers.");
        rf->permuteOriginalToReordered (OverlappingB, *reordered_B_);
        Inverse_->solve (*reordered_Y_, *reordered_B_);
        rf->permuteReorderedToOriginal (*reordered_Y_, OverlappingY);
      }
    }
  }
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
  using Teuchos::ParameterEntry;
  using Teuchos::ParameterEntryValidator;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Details::getParamTryingTypes;
  const char prefix[] = "Ifpack2::AdditiveSchwarz: ";

  if (plist.is_null ()) {
    // Assume that the user meant to set default parameters by passing
    // in an empty list.
    this->setParameterList (rcp (new ParameterList ()));
  }
  // FIXME (mfh 26 Aug 2015) It's not necessarily true that plist is
  // nonnull at this point.

  // At this point, plist should be nonnull.
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameterList: plist is null.  This should never happen, since the "
    "method should have replaced a null input list with a nonnull empty list "
    "by this point.  Please report this bug to the Ifpack2 developers.");

  // TODO JJH 24March2015  The list needs to be validated.  Not sure why this is commented out.
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

  const std::string cmParamName ("schwarz: combine mode");
  const ParameterEntry* cmEnt = plist->getEntryPtr (cmParamName);
  if (cmEnt != nullptr) {
    if (cmEnt->isType<CombineMode> ()) {
      CombineMode_ = Teuchos::getValue<CombineMode> (*cmEnt);
    }
    else if (cmEnt->isType<int> ()) {
      const int cm = Teuchos::getValue<int> (*cmEnt);
      CombineMode_ = static_cast<CombineMode> (cm);
    }
    else if (cmEnt->isType<std::string> ()) {
      // Try to get the combine mode as a string.  If this works, use
      // the validator to convert to int.  This is painful, but
      // necessary in order to do validation, since the input list may
      // not necessarily come with a validator.
      const ParameterEntry& validEntry =
        getValidParameters ()->getEntry (cmParamName);
      RCP<const ParameterEntryValidator> v = validEntry.validator ();
      using vs2e_type = StringToIntegralParameterEntryValidator<CombineMode>;
      RCP<const vs2e_type> vs2e = rcp_dynamic_cast<const vs2e_type> (v, true);

      ParameterEntry& inputEntry = plist->getEntry (cmParamName);
      // As AVG is only a Schwarz option and does not exist in Tpetra's
      // version of CombineMode, we use a separate boolean local to
      // Schwarz in conjunction with CombineMode_ == ADD to handle
      // averaging. Here, we change input entry to ADD and set the boolean.
      if (strncmp(Teuchos::getValue<std::string>(inputEntry).c_str(),"AVG",3) == 0) {
        inputEntry.template setValue<std::string>("ADD");
        AvgOverlap_ = true;
      }
      CombineMode_ = vs2e->getIntegralValue (inputEntry, cmParamName);
    }
  }
  // If doing user partitioning with Block Jacobi relaxation and overlapping blocks, we might
  // later need to know whether or not the overlapping Schwarz scheme is "ADD" or "ZERO" (which
  // is really RAS Schwarz. If it is "ADD", communication will be necessary when computing the
  // proper weights needed to combine solution values in overlap regions
  if (plist->isParameter("subdomain solver name")) {
    if (plist->get<std::string>("subdomain solver name") == "BLOCK_RELAXATION") {
      if (plist->isSublist("subdomain solver parameters")) {
        if (plist->sublist("subdomain solver parameters").isParameter("relaxation: type")) {
          if (plist->sublist("subdomain solver parameters").get<std::string>("relaxation: type") == "Jacobi" ) {
            if (plist->sublist("subdomain solver parameters").isParameter("partitioner: type")) {
              if (plist->sublist("subdomain solver parameters").get<std::string>("partitioner: type") == "user") { 
                 if (CombineMode_ == Tpetra::ADD)  plist->sublist("subdomain solver parameters").set("partitioner: combine mode","ADD");
                 if (CombineMode_ == Tpetra::ZERO) plist->sublist("subdomain solver parameters").set("partitioner: combine mode","ZERO");
                 AvgOverlap_ = false;     // averaging already taken care of by  the partitioner: nonsymmetric overlap combine option
              }
            }
          }   
        } 
      }
    }
  }

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

  // if true, filter singletons. NOTE: the filtered matrix can still have
  // singletons! A simple example: upper triangular matrix, if I remove
  // the lower node, I still get a matrix with a singleton! However, filter
  // singletons should help for PDE problems with Dirichlet BCs.
  FilterSingletons_ = plist->get ("schwarz: filter singletons", FilterSingletons_);

  // Allow for damped Schwarz updates
  getParamTryingTypes<scalar_type, scalar_type, double>
    (UpdateDamping_, *plist, "schwarz: update damping", prefix);

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
      std::pair<ParameterList, bool> result = innerPrecParams ();
      if (result.second) {
        // FIXME (mfh 26 Aug 2015) Rewrite innerPrecParams() so this
        // isn't another deep copy.
        Inverse_->setParameters (rcp (new ParameterList (result.first)));
      }
    }
  }

  NumIterations_ = plist->get ("schwarz: num iterations", NumIterations_);
  ZeroStartingSolution_ =
    plist->get ("schwarz: zero starting solution", ZeroStartingSolution_);
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
    const int  overlapLevel         = 0;
    const bool useReordering        = false;
    const bool filterSingletons     = false;
    const int  numIterations        = 1;
    const bool zeroStartingSolution = true;
    const scalar_type updateDamping = Teuchos::ScalarTraits<scalar_type>::one ();
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
    plist->set ("schwarz: num iterations", numIterations);
    plist->set ("schwarz: zero starting solution", zeroStartingSolution);
    plist->set ("schwarz: update damping", updateDamping);

    // FIXME (mfh 18 Nov 2013) Get valid parameters from inner solver.
    //        JJH The inner solver should handle its own validation.
    //
    // FIXME (mfh 18 Nov 2013) Get valid parameters from Zoltan2, if
    // Zoltan2 was enabled in the build.
    //        JJH Zoltan2 should handle its own validation.
    //

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
  double startTime = timer->wallTime();

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
      "initialize: The matrix to precondition is null.  You must either pass "
      "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
      "input, before you may call this method.");

    IsInitialized_ = false;
    IsComputed_ = false;
    overlapping_B_.reset (nullptr);
    overlapping_Y_.reset (nullptr);
    R_.reset (nullptr);
    C_.reset (nullptr);
    reduced_reordered_B_.reset (nullptr);
    reduced_reordered_Y_.reset (nullptr);
    reduced_B_.reset (nullptr);
    reduced_Y_.reset (nullptr);
    reordered_B_.reset (nullptr);
    reordered_Y_.reset (nullptr);

    RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
    RCP<const map_type> rowMap = Matrix_->getRowMap ();
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
        rcp (new map_type (INVALID, rowMap->getLocalNumElements (),
                           indexBase, localComm));
    }

    // compute the overlapping matrix if necessary
    if (IsOverlapping_) {
      Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("OverlappingRowMatrix construction"));
      OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_));
    }

    setup (); // This does a lot of the initialization work.

    if (! Inverse_.is_null ()) {
      Inverse_->symbolic (); // Initialize subdomain solver.
    }

  } // Stop timing here.

  IsInitialized_ = true;
  ++NumInitialize_;

  InitializeTime_ += (timer->wallTime() - startTime);
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
  TimeMonitor timeMon (*timer);
  double startTime = timer->wallTime();

  // compute () assumes that the values of Matrix_ (aka A) have changed.
  // If this has overlap, do an import from the input matrix to the halo.
  if (IsOverlapping_) {
    Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Halo Import"));
    OverlappingMatrix_->doExtImport();
  }
  // At this point, either Matrix_ or OverlappingMatrix_ (depending on whether this is overlapping)
  // has new values and unchanged structure. If we are using AdditiveSchwarzFilter, update the local matrix.
  //
  if(auto asf = Teuchos::rcp_dynamic_cast<Details::AdditiveSchwarzFilter<MatrixType>>(innerMatrix_))
  {
    Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Fill Local Matrix"));
    //NOTE: if this compute() call comes right after the initialize() with no intervening matrix changes, this call is redundant.
    //initialize() already filled the local matrix. However, we have no way to tell if this is the case.
    asf->updateMatrixValues();
  }
  //Now, whether the Inverse_'s matrix is the AdditiveSchwarzFilter's local matrix or simply Matrix_/OverlappingMatrix_,
  //it will be able to see the new values and update itself accordingly.

  { // Start timing here.

    IsComputed_ = false;
    Inverse_->numeric ();
  } // Stop timing here.

  IsComputed_ = true;
  ++NumCompute_;

  ComputeTime_ += (timer->wallTime() - startTime);
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
      << ", Iterations: " << NumIterations_
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
  if (! Inverse_.is_null ()) {
    Teuchos::RCP<Teuchos::Describable> inv =
      Teuchos::rcp_dynamic_cast<Teuchos::Describable> (Inverse_);
    if (! inv.is_null ()) {
      out << "{" << inv->description () << "}";
    } else {
      out << "{" << "Some inner solver" << "}";
    }
  } else {
    out << "null";
  }

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
          auto ifpack2_inverse = Teuchos::rcp_dynamic_cast< Ifpack2::Details::LinearSolver<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > (Inverse_);
          if(ifpack2_inverse.is_null()) 
            out << "Subdomain solver: not null" << endl;
          else {
            out << "Subdomain solver: "; ifpack2_inverse->describe (out, Teuchos::VERB_LOW);
          }
        }
      }
      else { // vl >= Teuchos::VERB_EXTREME
        for (int p = 0; p < numProcs; ++p) {
          if (p == myRank) {
            out << "Subdomain solver on Process " << myRank << ":";
            if (Inverse_.is_null ()) {
              out << "null" << endl;
            } else {
              Teuchos::RCP<Teuchos::Describable> inv =
                Teuchos::rcp_dynamic_cast<Teuchos::Describable> (Inverse_);
              if (! inv.is_null ()) {
                out << endl;
                inv->describe (out, vl);
              } else {
                out << "null" << endl;
              }
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
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::"
    "initialize: The matrix to precondition is null.  You must either pass "
    "a nonnull matrix to the constructor, or call setMatrix() with a nonnull "
    "input, before you may call this method.");
  
  //If the matrix is a CrsMatrix or OverlappingRowMatrix, use the high-performance 
  //AdditiveSchwarzFilter. Otherwise, use composition of Reordered/Singleton/LocalFilter.
  auto matrixCrs = rcp_dynamic_cast<const crs_matrix_type>(Matrix_);
  if(!OverlappingMatrix_.is_null() || !matrixCrs.is_null())
  {
    ArrayRCP<local_ordinal_type> perm;
    ArrayRCP<local_ordinal_type> revperm;
    if (UseReordering_) {
      Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Reordering"));
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
      // Unlike Ifpack, Zoltan2 does all the dirty work here.
      Teuchos::ParameterList zlist = List_.sublist ("schwarz: reordering list");
      ReorderingAlgorithm_ = zlist.get<std::string> ("order_method", "rcm");

      if(ReorderingAlgorithm_ == "user") {
        // User-provided reordering
        perm    = zlist.get<Teuchos::ArrayRCP<local_ordinal_type> >("user ordering");
        revperm = zlist.get<Teuchos::ArrayRCP<local_ordinal_type> >("user reverse ordering");
      }
      else {
        // Zoltan2 reordering
        typedef Tpetra::RowGraph
          <local_ordinal_type, global_ordinal_type, node_type> row_graph_type;
        typedef Zoltan2::TpetraRowGraphAdapter<row_graph_type> z2_adapter_type;
        auto constActiveGraph = Teuchos::rcp_const_cast<const row_graph_type>(
            IsOverlapping_ ? OverlappingMatrix_->getGraph() : Matrix_->getGraph());
        z2_adapter_type Zoltan2Graph (constActiveGraph);

        typedef Zoltan2::OrderingProblem<z2_adapter_type> ordering_problem_type;
#ifdef HAVE_MPI
        // Grab the MPI Communicator and build the ordering problem with that
        MPI_Comm myRawComm;

        RCP<const MpiComm<int> > mpicomm =
          rcp_dynamic_cast<const MpiComm<int> > (Matrix_->getComm ());
        if (mpicomm == Teuchos::null) {
          myRawComm = MPI_COMM_SELF;
        } else {
          myRawComm = * (mpicomm->getRawMpiComm ());
        }
        ordering_problem_type MyOrderingProblem (&Zoltan2Graph, &zlist, myRawComm);
#else
        ordering_problem_type MyOrderingProblem (&Zoltan2Graph, &zlist);
#endif
        MyOrderingProblem.solve ();

        {
          typedef Zoltan2::LocalOrderingSolution<local_ordinal_type>
            ordering_solution_type;

          ordering_solution_type sol (*MyOrderingProblem.getLocalOrderingSolution());

          // perm[i] gives the where OLD index i shows up in the NEW
          // ordering.  revperm[i] gives the where NEW index i shows
          // up in the OLD ordering.  Note that perm is actually the
          // "inverse permutation," in Zoltan2 terms.
          perm = sol.getPermutationRCPConst (true);
          revperm = sol.getPermutationRCPConst ();
        }
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
    else
    {
      local_ordinal_type numLocalRows = OverlappingMatrix_.is_null() ? matrixCrs->getLocalNumRows() : OverlappingMatrix_->getLocalNumRows();
      //Use an identity ordering.
      //TODO: create a non-permuted code path in AdditiveSchwarzFilter, in the case that neither
      //reordering nor singleton filtering are enabled. In this situation it's like LocalFilter.
      perm = ArrayRCP<local_ordinal_type>(numLocalRows);
      revperm = ArrayRCP<local_ordinal_type>(numLocalRows);
      for(local_ordinal_type i = 0; i < numLocalRows; i++)
      {
        perm[i] = i;
        revperm[i] = i;
      }
    }
    //Now, construct the filter
    {
      Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Filter construction"));
      RCP<Details::AdditiveSchwarzFilter<MatrixType>> asf;
      if(OverlappingMatrix_.is_null())
        asf = rcp(new Details::AdditiveSchwarzFilter<MatrixType>(matrixCrs, perm, revperm, FilterSingletons_));
      else
        asf = rcp(new Details::AdditiveSchwarzFilter<MatrixType>(OverlappingMatrix_, perm, revperm, FilterSingletons_));
      innerMatrix_ = asf;
    }
  }
  else
  {
    // Localized version of Matrix_ or OverlappingMatrix_.
    RCP<row_matrix_type> LocalizedMatrix;

    // The "most current local matrix."  At the end of this method, this
    // will be handed off to the inner solver.
    RCP<row_matrix_type> ActiveMatrix;

    // Create localized matrix.
    if (! OverlappingMatrix_.is_null ()) {
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (OverlappingMatrix_));
    }
    else {
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (Matrix_));
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
      typedef ReorderFilter<row_matrix_type> reorder_filter_type;
      Teuchos::ParameterList zlist = List_.sublist ("schwarz: reordering list");
      ReorderingAlgorithm_ = zlist.get<std::string> ("order_method", "rcm");

      ArrayRCP<local_ordinal_type> perm;
      ArrayRCP<local_ordinal_type> revperm;

      if(ReorderingAlgorithm_ == "user") {
        // User-provided reordering
        perm    = zlist.get<Teuchos::ArrayRCP<local_ordinal_type> >("user ordering");
        revperm = zlist.get<Teuchos::ArrayRCP<local_ordinal_type> >("user reverse ordering");
      }
      else {
        // Zoltan2 reordering
        typedef Tpetra::RowGraph
          <local_ordinal_type, global_ordinal_type, node_type> row_graph_type;
        typedef Zoltan2::TpetraRowGraphAdapter<row_graph_type> z2_adapter_type;
        RCP<const row_graph_type> constActiveGraph =
          Teuchos::rcp_const_cast<const row_graph_type>(ActiveMatrix->getGraph());
        z2_adapter_type Zoltan2Graph (constActiveGraph);

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
        ordering_problem_type MyOrderingProblem (&Zoltan2Graph, &zlist, myRawComm);
#else
        ordering_problem_type MyOrderingProblem (&Zoltan2Graph, &zlist);
#endif
        MyOrderingProblem.solve ();

        {
          typedef Zoltan2::LocalOrderingSolution<local_ordinal_type>
            ordering_solution_type;

          ordering_solution_type sol (*MyOrderingProblem.getLocalOrderingSolution());

          // perm[i] gives the where OLD index i shows up in the NEW
          // ordering.  revperm[i] gives the where NEW index i shows
          // up in the OLD ordering.  Note that perm is actually the
          // "inverse permutation," in Zoltan2 terms.
          perm = sol.getPermutationRCPConst (true);
          revperm = sol.getPermutationRCPConst ();
        }
      }
      // All reorderings here...
      ReorderedLocalizedMatrix_ = rcp (new reorder_filter_type (ActiveMatrix, perm, revperm));


      ActiveMatrix = ReorderedLocalizedMatrix_;
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
  }

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
      << Teuchos::TypeNameTraits<LocalInverseType>::name () << "\".  "
      "Please talk to the Ifpack2 developers for details.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      innerName == "CUSTOM", std::runtime_error, "Ifpack2::AdditiveSchwarz::"
      "initialize: If the \"inner preconditioner name\" parameter (or any "
      "alias thereof) has the value \"CUSTOM\", then you must first call "
      "setInnerPreconditioner with a nonnull inner preconditioner input before "
      "you may call initialize().");

    // FIXME (mfh 26 Aug 2015) Once we fix Bug 6392, the following
    // three lines of code can and SHOULD go away.
    if (! Trilinos::Details::Impl::registeredSomeLinearSolverFactory ("Ifpack2")) {
      Ifpack2::Details::registerLinearSolverFactory ();
    }

    // FIXME (mfh 26 Aug 2015) Provide the capability to get inner
    // solvers from packages other than Ifpack2.
    typedef typename MV::mag_type MT;
    RCP<inner_solver_type> innerPrec =
      Trilinos::Details::getLinearSolver<MV, OP, MT> ("Ifpack2", innerName);
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerPrec.is_null (), std::logic_error,
      "Ifpack2::AdditiveSchwarz::setup: Failed to create inner preconditioner "
      "with name \"" << innerName << "\".");
    innerPrec->setMatrix (innerMatrix_);

    // Extract and apply the sublist of parameters to give to the
    // inner solver, if there is such a sublist of parameters.
    std::pair<Teuchos::ParameterList, bool> result = innerPrecParams ();
    if (result.second) {
      // FIXME (mfh 26 Aug 2015) We don't really want to use yet
      // another deep copy of the ParameterList here.
      innerPrec->setParameters (rcp (new ParameterList (result.first)));
    }
    Inverse_ = innerPrec; // "Commit" the inner solver.
  }
  else if (Inverse_->getMatrix ().getRawPtr () != innerMatrix_.getRawPtr ()) {
    // The new inner matrix is different from the inner
    // preconditioner's current matrix, so give the inner
    // preconditioner the new inner matrix.
    Inverse_->setMatrix (innerMatrix_);
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

    // mfh 03 Jan 2014: Thanks to Paul Tsuji for pointing out that
    // it's perfectly legal for innerMatrix_ to be null here.  This
    // can happen if initialize() has not been called yet.  For
    // example, when Ifpack2::Factory creates an AdditiveSchwarz
    // instance, it calls setInnerPreconditioner() without first
    // calling initialize().

    // Give the local matrix to the new inner solver.
    if(auto asf = Teuchos::rcp_dynamic_cast<Details::AdditiveSchwarzFilter<MatrixType>>(innerMatrix_))
      innerSolver->setMatrix (asf->getFilteredMatrix());
    else
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
  typedef Ifpack2::Details::LinearSolver<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> inner_solver_impl_type;
  Inverse_ = Teuchos::rcp (new inner_solver_impl_type (innerPrec, "CUSTOM"));
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
    overlapping_B_.reset (nullptr);
    overlapping_Y_.reset (nullptr);
    R_.reset (nullptr);
    C_.reset (nullptr);
    DistributedImporter_ = Teuchos::null;

    Matrix_ = A;
  }
}

} // namespace Ifpack2

// NOTE (mfh 26 Aug 2015) There's no need to instantiate for CrsMatrix
// too.  All Ifpack2 preconditioners can and should do dynamic casts
// internally, if they need a type more specific than RowMatrix.
#define IFPACK2_ADDITIVESCHWARZ_INSTANT(S,LO,GO,N) \
  template class Ifpack2::AdditiveSchwarz< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP

