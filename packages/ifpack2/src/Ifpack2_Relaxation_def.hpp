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

#ifndef IFPACK2_RELAXATION_DEF_HPP
#define IFPACK2_RELAXATION_DEF_HPP

#include "Ifpack2_Relaxation_decl.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

// mfh 28 Mar 2013: Uncomment out these three lines to compute
// statistics on diagonal entries in compute().
// #ifndef IFPACK2_RELAXATION_COMPUTE_DIAGONAL_STATS
// #  define IFPACK2_RELAXATION_COMPUTE_DIAGONAL_STATS 1
// #endif // IFPACK2_RELAXATION_COMPUTE_DIAGONAL_STATS

namespace {
  // Validate that a given int is nonnegative.
  class NonnegativeIntValidator : public Teuchos::ParameterEntryValidator {
  public:
    // Constructor (does nothing).
    NonnegativeIntValidator () {}

    // ParameterEntryValidator wants this method.
    Teuchos::ParameterEntryValidator::ValidStringsList validStringValues () const {
      return Teuchos::null;
    }

    // Actually validate the parameter's value.
    void
    validate (const Teuchos::ParameterEntry& entry,
              const std::string& paramName,
              const std::string& sublistName) const
    {
      using std::endl;
      Teuchos::any anyVal = entry.getAny (true);
      const std::string entryName = entry.getAny (false).typeName ();

      TEUCHOS_TEST_FOR_EXCEPTION(
        anyVal.type () != typeid (int),
        Teuchos::Exceptions::InvalidParameterType,
        "Parameter \"" << paramName << "\" in sublist \"" << sublistName
        << "\" has the wrong type." << endl << "Parameter: " << paramName
        << endl << "Type specified: " << entryName << endl
        << "Type required: int" << endl);

      const int val = Teuchos::any_cast<int> (anyVal);
      TEUCHOS_TEST_FOR_EXCEPTION(
        val < 0, Teuchos::Exceptions::InvalidParameterValue,
        "Parameter \"" << paramName << "\" in sublist \"" << sublistName
        << "\" is negative." << endl << "Parameter: " << paramName
        << endl << "Value specified: " << val << endl
        << "Required range: [0, INT_MAX]" << endl);
    }

    // ParameterEntryValidator wants this method.
    const std::string getXMLTypeName () const {
      return "NonnegativeIntValidator";
    }

    // ParameterEntryValidator wants this method.
    void
    printDoc (const std::string& docString,
              std::ostream &out) const
    {
      Teuchos::StrUtils::printLines (out, "# ", docString);
      out << "#\tValidator Used: " << std::endl;
      out << "#\t\tNonnegativeIntValidator" << std::endl;
    }
  };

  // A way to get a small positive number (eps() for floating-point
  // types, or 1 for integer types) when Teuchos::ScalarTraits doesn't
  // define it (for example, for integer values).
  template<class Scalar, const bool isOrdinal=Teuchos::ScalarTraits<Scalar>::isOrdinal>
  class SmallTraits {
  public:
    // Return eps if Scalar is a floating-point type, else return 1.
    static const Scalar eps ();
  };

  // Partial specialization for when Scalar is not a floating-point type.
  template<class Scalar>
  class SmallTraits<Scalar, true> {
  public:
    static const Scalar eps () {
      return Teuchos::ScalarTraits<Scalar>::one ();
    }
  };

  // Partial specialization for when Scalar is a floating-point type.
  template<class Scalar>
  class SmallTraits<Scalar, false> {
  public:
    static const Scalar eps () {
      return Teuchos::ScalarTraits<Scalar>::eps ();
    }
  };
} // namespace (anonymous)

namespace Ifpack2 {

template<class MatrixType>
void Relaxation<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != A_.getRawPtr ()) { // it's a different matrix
    Importer_ = Teuchos::null;
    Diagonal_ = Teuchos::null; // ??? what if this comes from the user???
    isInitialized_ = false;
    IsComputed_ = false;
    diagOffsets_ = Teuchos::null;
    savedDiagOffsets_ = false;
    hasBlockCrsMatrix_ = false;
    if (! A.is_null ()) {
      IsParallel_ = (A->getRowMap ()->getComm ()->getSize () > 1);
    }
    A_ = A;
  }
}


template<class MatrixType>
Relaxation<MatrixType>::
Relaxation (const Teuchos::RCP<const row_matrix_type>& A)
: A_ (A),
  Time_ (Teuchos::rcp (new Teuchos::Time ("Ifpack2::Relaxation"))),
  NumSweeps_ (1),
  PrecType_ (Ifpack2::Details::JACOBI),
  DampingFactor_ (STS::one ()),
  IsParallel_ (A.is_null () ? false : A->getRowMap ()->getComm ()->getSize () > 1), // pick a reasonable default if the matrix is null.  It doesn't really matter, because the methods that care can't use a null matrix anyway.
  ZeroStartingSolution_ (true),
  DoBackwardGS_ (false),
  DoL1Method_ (false),
  L1Eta_ (Teuchos::as<magnitude_type> (1.5)),
  MinDiagonalValue_ (STS::zero ()),
  fixTinyDiagEntries_ (false),
  checkDiagEntries_ (false),
  Condest_ (-STM::one ()),
  isInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0), // Times are double anyway, so no need for ScalarTraits.
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  ComputeFlops_ (0.0),
  ApplyFlops_ (0.0),
  globalMinMagDiagEntryMag_ (STM::zero ()),
  globalMaxMagDiagEntryMag_ (STM::zero ()),
  globalNumSmallDiagEntries_ (0),
  globalNumZeroDiagEntries_ (0),
  globalNumNegDiagEntries_ (0),
  globalDiagNormDiff_(Teuchos::ScalarTraits<magnitude_type>::zero()),
  savedDiagOffsets_ (false),
  hasBlockCrsMatrix_ (false)
{
  this->setObjectLabel ("Ifpack2::Relaxation");
}

//==========================================================================
template<class MatrixType>
Relaxation<MatrixType>::~Relaxation() {
}

template<class MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
Relaxation<MatrixType>::getValidParameters () const
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::setStringToIntegralParameter;
  typedef Teuchos::ParameterEntryValidator PEV;

  if (validParams_.is_null ()) {
    RCP<ParameterList> pl = parameterList ("Ifpack2::Relaxation");

    // Set a validator that automatically converts from the valid
    // string options to their enum values.
    Array<std::string> precTypes (3);
    precTypes[0] = "Jacobi";
    precTypes[1] = "Gauss-Seidel";
    precTypes[2] = "Symmetric Gauss-Seidel";
    Array<Details::RelaxationType> precTypeEnums (3);
    precTypeEnums[0] = Details::JACOBI;
    precTypeEnums[1] = Details::GS;
    precTypeEnums[2] = Details::SGS;
    const std::string defaultPrecType ("Jacobi");
    setStringToIntegralParameter<Details::RelaxationType> ("relaxation: type",
      defaultPrecType, "Relaxation method", precTypes (), precTypeEnums (),
      pl.getRawPtr ());

    const int numSweeps = 1;
    RCP<PEV> numSweepsValidator =
      rcp_implicit_cast<PEV> (rcp (new NonnegativeIntValidator));
    pl->set ("relaxation: sweeps", numSweeps, "Number of relaxation sweeps",
             rcp_const_cast<const PEV> (numSweepsValidator));

    const scalar_type dampingFactor = STS::one ();
    pl->set ("relaxation: damping factor", dampingFactor);

    const bool zeroStartingSolution = true;
    pl->set ("relaxation: zero starting solution", zeroStartingSolution);

    const bool doBackwardGS = false;
    pl->set ("relaxation: backward mode", doBackwardGS);

    const bool doL1Method = false;
    pl->set ("relaxation: use l1", doL1Method);

    const magnitude_type l1eta = (STM::one() + STM::one() + STM::one()) /
      (STM::one() + STM::one()); // 1.5
    pl->set ("relaxation: l1 eta", l1eta);

    const scalar_type minDiagonalValue = STS::zero ();
    pl->set ("relaxation: min diagonal value", minDiagonalValue);

    const bool fixTinyDiagEntries = false;
    pl->set ("relaxation: fix tiny diagonal entries", fixTinyDiagEntries);

    const bool checkDiagEntries = false;
    pl->set ("relaxation: check diagonal entries", checkDiagEntries);

    Teuchos::ArrayRCP<local_ordinal_type> localSmoothingIndices = Teuchos::null;
    pl->set("relaxation: local smoothing indices", localSmoothingIndices);

    validParams_ = rcp_const_cast<const ParameterList> (pl);
  }
  return validParams_;
}


template<class MatrixType>
void Relaxation<MatrixType>::setParametersImpl (Teuchos::ParameterList& pl)
{
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  typedef scalar_type ST; // just to make code below shorter

  pl.validateParametersAndSetDefaults (* getValidParameters ());

  const Details::RelaxationType precType =
    getIntegralValue<Details::RelaxationType> (pl, "relaxation: type");
  const int numSweeps = pl.get<int> ("relaxation: sweeps");
  const ST dampingFactor = pl.get<ST> ("relaxation: damping factor");
  const bool zeroStartSol = pl.get<bool> ("relaxation: zero starting solution");
  const bool doBackwardGS = pl.get<bool> ("relaxation: backward mode");
  const bool doL1Method = pl.get<bool> ("relaxation: use l1");
  const magnitude_type l1Eta = pl.get<magnitude_type> ("relaxation: l1 eta");
  const ST minDiagonalValue = pl.get<ST> ("relaxation: min diagonal value");
  const bool fixTinyDiagEntries = pl.get<bool> ("relaxation: fix tiny diagonal entries");
  const bool checkDiagEntries = pl.get<bool> ("relaxation: check diagonal entries");
  Teuchos::ArrayRCP<local_ordinal_type> localSmoothingIndices = pl.get<Teuchos::ArrayRCP<local_ordinal_type> >("relaxation: local smoothing indices");


  // "Commit" the changes, now that we've validated everything.
  PrecType_              = precType;
  NumSweeps_             = numSweeps;
  DampingFactor_         = dampingFactor;
  ZeroStartingSolution_  = zeroStartSol;
  DoBackwardGS_          = doBackwardGS;
  DoL1Method_            = doL1Method;
  L1Eta_                 = l1Eta;
  MinDiagonalValue_      = minDiagonalValue;
  fixTinyDiagEntries_    = fixTinyDiagEntries;
  checkDiagEntries_      = checkDiagEntries;
  localSmoothingIndices_ = localSmoothingIndices;
}


template<class MatrixType>
void Relaxation<MatrixType>::setParameters (const Teuchos::ParameterList& pl)
{
  // FIXME (aprokop 18 Oct 2013) Casting away const is bad here.
  // but otherwise, we will get [unused] in pl
  this->setParametersImpl(const_cast<Teuchos::ParameterList&>(pl));
}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Relaxation<MatrixType>::getComm() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::getComm: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getRowMap ()->getComm ();
}


template<class MatrixType>
Teuchos::RCP<const typename Relaxation<MatrixType>::row_matrix_type>
Relaxation<MatrixType>::getMatrix() const {
  return A_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
Relaxation<MatrixType>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::getDomainMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getDomainMap ();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
Relaxation<MatrixType>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::getRangeMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getRangeMap ();
}

//==========================================================================
template<class MatrixType>
bool Relaxation<MatrixType>::hasTransposeApply () const {
  return true;
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Relaxation<MatrixType>::getCondEst () const
{
  return Condest_;
}


template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Relaxation<MatrixType>::
computeCondEst (CondestType CT,
                typename MatrixType::local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) { // cannot compute right now
    return -Teuchos::ScalarTraits<magnitude_type>::one ();
  }
  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  return Condest_;
}


template<class MatrixType>
void
Relaxation<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::apply: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (),
    std::runtime_error,
    "Ifpack2::Relaxation::apply: You must call compute() on this Ifpack2 "
    "preconditioner instance before you may call apply().  You may call "
    "isComputed() to find out if compute() has been called already.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(),
    std::runtime_error,
    "Ifpack2::Relaxation::apply: X and Y have different numbers of columns.  "
    "X has " << X.getNumVectors() << " columns, but Y has "
    << Y.getNumVectors() << " columns.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    beta != STS::zero (), std::logic_error,
    "Ifpack2::Relaxation::apply: beta = " << beta << " != 0 case not "
    "implemented.");
  {
    // Reset the timer each time, since Relaxation uses the same Time
    // object to track times for different methods.
    Teuchos::TimeMonitor timeMon (*Time_, true);

    // Special case: alpha == 0.
    if (alpha == STS::zero ()) {
      // No floating-point operations, so no need to update a count.
      Y.putScalar (STS::zero ());
    }
    else {
      // If X and Y are pointing to the same memory location,
      // we need to create an auxiliary vector, Xcopy
      RCP<const MV> Xcopy;
      if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
        Xcopy = rcp (new MV(Tpetra::createCopy(X)));
      }
      else {
        Xcopy = rcpFromRef (X);
      }

      // Each of the following methods updates the flop count itself.
      // All implementations handle zeroing out the starting solution
      // (if necessary) themselves.
      switch (PrecType_) {
      case Ifpack2::Details::JACOBI:
        ApplyInverseJacobi(*Xcopy,Y);
        break;
      case Ifpack2::Details::GS:
        ApplyInverseGS(*Xcopy,Y);
        break;
      case Ifpack2::Details::SGS:
        ApplyInverseSGS(*Xcopy,Y);
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Ifpack2::Relaxation::apply: Invalid preconditioner type enum value "
          << PrecType_ << ".  Please report this bug to the Ifpack2 developers.");
      }
      if (alpha != STS::one ()) {
        Y.scale (alpha);
        const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
        const double numVectors = as<double> (Y.getNumVectors ());
        ApplyFlops_ += numGlobalRows * numVectors;
      }
    }
  }
  ApplyTime_ += Time_->totalElapsedTime ();
  ++NumApply_;
}


template<class MatrixType>
void
Relaxation<MatrixType>::
applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
          Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
          Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::applyMat: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error, "Ifpack2::Relaxation::applyMat: "
    "isComputed() must be true before you may call applyMat().  "
    "Please call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::Relaxation::applyMat: X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");
  A_->apply (X, Y, mode);
}


template<class MatrixType>
void Relaxation<MatrixType>::initialize () {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::initialize: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // Initialization for Relaxation is trivial, so we say it takes zero time.
  //InitializeTime_ += Time_->totalElapsedTime ();

  ++NumInitialize_;
  isInitialized_ = true;

  const block_crs_matrix_type* blockCrsA = dynamic_cast<const block_crs_matrix_type*>(A_.getRawPtr());
  if (blockCrsA != NULL)
    hasBlockCrsMatrix_ = true;

}

template<class MatrixType>
void Relaxation<MatrixType>::computeBlockCrs ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MAX;
  using Teuchos::REDUCE_MIN;
  using Teuchos::REDUCE_SUM;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::reduceAll;

  {
    // Reset the timer each time, since Relaxation uses the same Time
    // object to track times for different methods.
    Teuchos::TimeMonitor timeMon (*Time_, true);

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::"
      "computeBlockCrs: The input matrix A is null.  Please call setMatrix() "
      "with a nonnull input matrix, then call initialize(), before calling "
      "this method.");
    const block_crs_matrix_type* blockCrsA =
      dynamic_cast<const block_crs_matrix_type*> (A_.getRawPtr ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockCrsA == NULL, std::logic_error, "Ifpack2::Relaxation::"
      "computeBlockCrs: A_ is not a BlockCrsMatrix, but it should be if we "
      "got this far.  Please report this bug to the Ifpack2 developers.");

    const scalar_type one = STS::one ();

    // Reset state.
    IsComputed_ = false;
    Condest_ = -STM::one ();

    const local_ordinal_type blockSize = blockCrsA->getBlockSize ();

    if (! savedDiagOffsets_) {
      BlockDiagonal_ = Teuchos::null;
      BlockDiagonal_ =
        rcp (new block_crs_matrix_type (* (blockCrsA->getDiagonalGraph ()),
                                        blockSize));
      blockCrsA->getLocalDiagOffsets (diagOffsets_);
      savedDiagOffsets_ = true;
    }
    blockCrsA->getLocalDiagCopy (*BlockDiagonal_, diagOffsets_ ());

    const size_t numMyRows = A_->getNodeNumRows ();

    if (DoL1Method_ && IsParallel_) {
      const scalar_type two = one + one;
      const size_t maxLength = A_->getNodeMaxNumRowEntries ();
      Array<local_ordinal_type> indices (maxLength);
      Array<scalar_type> values (maxLength * blockSize * blockSize);
      size_t numEntries = 0;

      for (size_t i = 0; i < numMyRows; ++i) {
        blockCrsA->getLocalRowCopy (i, indices (), values (), numEntries);
        scalar_type* diagBlock = BlockDiagonal_->getLocalBlock (i,i).getRawPtr ();
        for (local_ordinal_type subRow = 0; subRow < blockSize; ++subRow) {
          magnitude_type diagonal_boost = STM::zero ();
          for (size_t k = 0 ; k < numEntries ; ++k) {
            if (static_cast<size_t> (indices[k]) > numMyRows) {
              const size_t offset = blockSize*blockSize*k + subRow*blockSize;
              for (local_ordinal_type subCol = 0; subCol < blockSize; ++subCol) {
                diagonal_boost += STS::magnitude (values[offset+subCol] / two);
              }
            }
          }
          if (STS::magnitude (diagBlock[subRow*(blockSize+1)]) < L1Eta_ * diagonal_boost) {
            diagBlock[subRow*(blockSize+1)] += diagonal_boost;
          }
        }
      }
    }

    blockDiagonalFactorizationPivots.resize (numMyRows * blockSize);
    int info;
    for (size_t i = 0 ; i < numMyRows; ++i) {
      typename block_crs_matrix_type::little_block_type diagBlock =
        BlockDiagonal_->getLocalBlock (i, i);
      diagBlock.factorize (&blockDiagonalFactorizationPivots[i*blockSize], info);
    }

    Importer_ = A_->getGraph ()->getImporter ();
  } // end TimeMonitor scope

  ComputeTime_ += Time_->totalElapsedTime ();
  ++NumCompute_;
  IsComputed_ = true;
}

template<class MatrixType>
void Relaxation<MatrixType>::compute ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MAX;
  using Teuchos::REDUCE_MIN;
  using Teuchos::REDUCE_SUM;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::reduceAll;
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vector_type;
  const scalar_type zero = STS::zero ();
  const scalar_type one = STS::one ();

  // We don't count initialization in compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  if (hasBlockCrsMatrix_)
  {
    computeBlockCrs();
    return;
  }

  {
    // Reset the timer each time, since Relaxation uses the same Time
    // object to track times for different methods.
    Teuchos::TimeMonitor timeMon (*Time_, true);

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Relaxation::compute: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix, then call initialize(), before calling this method.");


    // Reset state.
    IsComputed_ = false;
    Condest_ = -STM::one ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      NumSweeps_ < 0, std::logic_error,
      "Ifpack2::Relaxation::compute: NumSweeps_ = " << NumSweeps_ << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");

    Diagonal_ = rcp (new vector_type (A_->getRowMap ()));

    TEUCHOS_TEST_FOR_EXCEPTION(
      Diagonal_.is_null (), std::logic_error,
      "Ifpack2::Relaxation::compute: Vector of diagonal entries has not been "
      "created yet.  Please report this bug to the Ifpack2 developers.");

    // Extract the diagonal entries.  The CrsMatrix static graph
    // version is faster for subsequent calls to compute(), since it
    // caches the diagonal offsets.
    //
    // isStaticGraph() == true means that the matrix was created with
    // a const graph.  The only requirement is that the structure of
    // the matrix never changes, so isStaticGraph() == true is a bit
    // more conservative than we need.  However, Tpetra doesn't (as of
    // 05 Apr 2013) have a way to tell if the graph hasn't changed
    // since the last time we used it.
    {
      // NOTE (mfh 07 Jul 2013): We must cast here to CrsMatrix
      // instead of MatrixType, because isStaticGraph is a CrsMatrix
      // method (not inherited from RowMatrix's interface).  It's
      // perfectly valid to do relaxation on a RowMatrix which is not
      // a CrsMatrix.
      const crs_matrix_type* crsMat =
        dynamic_cast<const crs_matrix_type*> (A_.getRawPtr ());
      if (crsMat == NULL || ! crsMat->isStaticGraph ()) {
        A_->getLocalDiagCopy (*Diagonal_); // slow path
      } else {
        if (! savedDiagOffsets_) { // we haven't precomputed offsets
          crsMat->getLocalDiagOffsets (diagOffsets_);
          savedDiagOffsets_ = true;
        }
        crsMat->getLocalDiagCopy (*Diagonal_, diagOffsets_ ());
#ifdef HAVE_TPETRA_DEBUG
        // Validate the fast-path diagonal against the slow-path diagonal.
        vector_type D_copy (A_->getRowMap ());
        A_->getLocalDiagCopy (D_copy);
        D_copy.update (STS::one (), *Diagonal_, -STS::one ());
        const magnitude_type err = D_copy.normInf ();
        // The two diagonals should be exactly the same, so their
        // difference should be exactly zero.
        TEUCHOS_TEST_FOR_EXCEPTION(
          err != STM::zero(), std::logic_error, "Ifpack2::Relaxation::compute: "
          "\"fast-path\" diagonal computation failed.  \\|D1 - D2\\|_inf = "
          << err << ".");
#endif // HAVE_TPETRA_DEBUG
      }
    }

    // If we're checking the computed inverse diagonal, then keep a
    // copy of the original diagonal entries for later comparison.
    RCP<vector_type> origDiag;
    if (checkDiagEntries_) {
      origDiag = rcp (new vector_type (A_->getRowMap ()));
      *origDiag = *Diagonal_;
    }

    // "Host view" means that if the Node type is a GPU Node, the
    // ArrayRCP points to host memory, not device memory.  It will
    // write back to device memory (Diagonal_) at end of scope.
    ArrayRCP<scalar_type> diagHostView = Diagonal_->get1dViewNonConst ();

    // The view below is only valid as long as diagHostView is within
    // scope.  We extract a raw pointer in release mode because we
    // don't trust older versions of compilers (like GCC 4.4.x or
    // Intel < 13) to optimize away ArrayView::operator[].
#ifdef HAVE_IFPACK2_DEBUG
    ArrayView<scalar_type> diag = diagHostView ();
#else
    scalar_type* KOKKOSCLASSIC_RESTRICT const diag = diagHostView.getRawPtr ();
#endif // HAVE_IFPACK2_DEBUG

    const size_t numMyRows = A_->getNodeNumRows ();

    // Setup for L1 Methods.
    // Here we add half the value of the off-processor entries in the row,
    // but only if diagonal isn't sufficiently large.
    //
    // This follows from Equation (6.5) in: Baker, Falgout, Kolev and
    // Yang.  "Multigrid Smoothers for Ultraparallel Computing."  SIAM
    // J. Sci. Comput., Vol. 33, No. 5. (2011), pp. 2864-2887.
    if (DoL1Method_ && IsParallel_) {
      const scalar_type two = one + one;
      const size_t maxLength = A_->getNodeMaxNumRowEntries ();
      Array<local_ordinal_type> indices (maxLength);
      Array<scalar_type> values (maxLength);
      size_t numEntries;

      for (size_t i = 0; i < numMyRows; ++i) {
        A_->getLocalRowCopy (i, indices (), values (), numEntries);
        magnitude_type diagonal_boost = STM::zero ();
        for (size_t k = 0 ; k < numEntries ; ++k) {
          if (static_cast<size_t> (indices[k]) > numMyRows) {
            diagonal_boost += STS::magnitude (values[k] / two);
          }
        }
        if (STS::magnitude (diag[i]) < L1Eta_ * diagonal_boost) {
          diag[i] += diagonal_boost;
        }
      }

    }

    //
    // Invert the diagonal entries of the matrix (not in place).
    //

    // Precompute some quantities for "fixing" zero or tiny diagonal
    // entries.  We'll only use them if this "fixing" is enabled.
    //
    // SmallTraits covers for the lack of eps() in
    // Teuchos::ScalarTraits when its template parameter is not a
    // floating-point type.  (Ifpack2 sometimes gets instantiated for
    // integer Scalar types.)
    const scalar_type oneOverMinDiagVal = (MinDiagonalValue_ == zero) ?
      one / SmallTraits<scalar_type>::eps () :
      one / MinDiagonalValue_;
    // It's helpful not to have to recompute this magnitude each time.
    const magnitude_type minDiagValMag = STS::magnitude (MinDiagonalValue_);

    if (checkDiagEntries_) {
      //
      // Check diagonal elements, replace zero elements with the minimum
      // diagonal value, and store their inverses.  Count the number of
      // "small" and zero diagonal entries while we're at it.
      //
      size_t numSmallDiagEntries = 0; // "small" includes zero
      size_t numZeroDiagEntries = 0; // # zero diagonal entries
      size_t numNegDiagEntries = 0; // # negative (real parts of) diagonal entries

      // As we go, keep track of the diagonal entries with the least and
      // greatest magnitude.  We could use the trick of starting the min
      // with +Inf and the max with -Inf, but that doesn't work if
      // scalar_type is a built-in integer type.  Thus, we have to start
      // by reading the first diagonal entry redundantly.
      // scalar_type minMagDiagEntry = zero;
      // scalar_type maxMagDiagEntry = zero;
      magnitude_type minMagDiagEntryMag = STM::zero ();
      magnitude_type maxMagDiagEntryMag = STM::zero ();
      if (numMyRows > 0) {
        const scalar_type d_0 = diag[0];
        const magnitude_type d_0_mag = STS::magnitude (d_0);
        // minMagDiagEntry = d_0;
        // maxMagDiagEntry = d_0;
        minMagDiagEntryMag = d_0_mag;
        maxMagDiagEntryMag = d_0_mag;
      }

      // Go through all the diagonal entries.  Compute counts of
      // small-magnitude, zero, and negative-real-part entries.  Invert
      // the diagonal entries that aren't too small.  For those that are
      // too small in magnitude, replace them with 1/MinDiagonalValue_
      // (or 1/eps if MinDiagonalValue_ happens to be zero).
      for (size_t i = 0 ; i < numMyRows; ++i) {
        const scalar_type d_i = diag[i];
        const magnitude_type d_i_mag = STS::magnitude (d_i);
        const magnitude_type d_i_real = STS::real (d_i);

        // We can't compare complex numbers, but we can compare their
        // real parts.
        if (d_i_real < STM::zero ()) {
          ++numNegDiagEntries;
        }
        if (d_i_mag < minMagDiagEntryMag) {
          // minMagDiagEntry = d_i;
          minMagDiagEntryMag = d_i_mag;
        }
        if (d_i_mag > maxMagDiagEntryMag) {
          // maxMagDiagEntry = d_i;
          maxMagDiagEntryMag = d_i_mag;
        }

        if (fixTinyDiagEntries_) {
          // <= not <, in case minDiagValMag is zero.
          if (d_i_mag <= minDiagValMag) {
            ++numSmallDiagEntries;
            if (d_i_mag == STM::zero ()) {
              ++numZeroDiagEntries;
            }
            diag[i] = oneOverMinDiagVal;
          } else {
            diag[i] = one / d_i;
          }
        }
        else { // Don't fix zero or tiny diagonal entries.
          // <= not <, in case minDiagValMag is zero.
          if (d_i_mag <= minDiagValMag) {
            ++numSmallDiagEntries;
            if (d_i_mag == STM::zero ()) {
              ++numZeroDiagEntries;
            }
          }
          diag[i] = one / d_i;
        }
      }

      // We're done computing the inverse diagonal, so invalidate the view.
      // This ensures that the operations below, that use methods on Vector,
      // produce correct results even on a GPU Node.
      diagHostView = Teuchos::null;

      // Count floating-point operations of computing the inverse diagonal.
      //
      // FIXME (mfh 30 Mar 2013) Shouldn't counts be global, not local?
      if (STS::isComplex) { // magnitude: at least 3 flops per diagonal entry
        ComputeFlops_ += 4.0 * numMyRows;
      } else {
        ComputeFlops_ += numMyRows;
      }

      // Collect global data about the diagonal entries.  Only do this
      // (which involves O(1) all-reduces) if the user asked us to do
      // the extra work.
      //
      // FIXME (mfh 28 Mar 2013) This is wrong if some processes have
      // zero rows.  Fixing this requires one of two options:
      //
      // 1. Use a custom reduction operation that ignores processes
      //    with zero diagonal entries.
      // 2. Split the communicator, compute all-reduces using the
      //    subcommunicator over processes that have a nonzero number
      //    of diagonal entries, and then broadcast from one of those
      //    processes (if there is one) to the processes in the other
      //    subcommunicator.
      const Comm<int>& comm = * (A_->getRowMap ()->getComm ());

      // Compute global min and max magnitude of entries.
      Array<magnitude_type> localVals (2);
      localVals[0] = minMagDiagEntryMag;
      // (- (min (- x))) is the same as (max x).
      localVals[1] = -maxMagDiagEntryMag;
      Array<magnitude_type> globalVals (2);
      reduceAll<int, magnitude_type> (comm, REDUCE_MIN, 2,
                                      localVals.getRawPtr (),
                                      globalVals.getRawPtr ());
      globalMinMagDiagEntryMag_ = globalVals[0];
      globalMaxMagDiagEntryMag_ = -globalVals[1];

      Array<size_t> localCounts (3);
      localCounts[0] = numSmallDiagEntries;
      localCounts[1] = numZeroDiagEntries;
      localCounts[2] = numNegDiagEntries;
      Array<size_t> globalCounts (3);
      reduceAll<int, size_t> (comm, REDUCE_SUM, 3,
                              localCounts.getRawPtr (),
                              globalCounts.getRawPtr ());
      globalNumSmallDiagEntries_ = globalCounts[0];
      globalNumZeroDiagEntries_ = globalCounts[1];
      globalNumNegDiagEntries_ = globalCounts[2];

      // Forestall "set but not used" compiler warnings.
      // (void) minMagDiagEntry;
      // (void) maxMagDiagEntry;

      // Compute and save the difference between the computed inverse
      // diagonal, and the original diagonal's inverse.
      RCP<vector_type> diff = rcp (new vector_type (A_->getRowMap ()));
      diff->reciprocal (*origDiag);
      diff->update (-one, *Diagonal_, one);
      globalDiagNormDiff_ = diff->norm2 ();
    }
    else { // don't check diagonal elements
      if (fixTinyDiagEntries_) {
        // Go through all the diagonal entries.  Invert those that
        // aren't too small in magnitude.  For those that are too
        // small in magnitude, replace them with oneOverMinDiagVal.
        for (size_t i = 0 ; i < numMyRows; ++i) {
          const scalar_type d_i = diag[i];
          const magnitude_type d_i_mag = STS::magnitude (d_i);

          // <= not <, in case minDiagValMag is zero.
          if (d_i_mag <= minDiagValMag) {
            diag[i] = oneOverMinDiagVal;
          } else {
            diag[i] = one / d_i;
          }
        }
      }
      else { // don't fix tiny or zero diagonal entries
        for (size_t i = 0 ; i < numMyRows; ++i) {
          diag[i] = one / diag[i];
        }
      }
      if (STS::isComplex) { // magnitude: at least 3 flops per diagonal entry
        ComputeFlops_ += 4.0 * numMyRows;
      } else {
        ComputeFlops_ += numMyRows;
      }
    }

    if (IsParallel_ && (PrecType_ == Ifpack2::Details::GS ||
                        PrecType_ == Ifpack2::Details::SGS)) {
      Importer_ = A_->getGraph ()->getImporter ();
      // mfh 21 Mar 2013: The Import object may be null, but in that
      // case, the domain and column Maps are the same and we don't
      // need to Import anyway.
    }
  } // end TimeMonitor scope

  ComputeTime_ += Time_->totalElapsedTime ();
  ++NumCompute_;
  IsComputed_ = true;
}


template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseJacobi (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::as;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> MV;

  if (hasBlockCrsMatrix_) {
    ApplyInverseJacobi_BlockCrsMatrix (X, Y);
    return;
  }

  const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
  const double numVectors = as<double> (X.getNumVectors ());
  if (ZeroStartingSolution_) {
    // For the first Jacobi sweep, if we are allowed to assume that
    // the initial guess is zero, then Jacobi is just diagonal
    // scaling.  (A_ij * x_j = 0 for i != j, since x_j = 0.)
    // Compute it as Y(i,j) = DampingFactor_ * X(i,j) * D(i).
    Y.elementWiseMultiply (DampingFactor_, *Diagonal_, X, STS::zero ());

    // Count (global) floating-point operations.  Ifpack2 represents
    // this as a floating-point number rather than an integer, so that
    // overflow (for a very large number of calls, or a very large
    // problem) is approximate instead of catastrophic.
    double flopUpdate = 0.0;
    if (DampingFactor_ == STS::one ()) {
      // Y(i,j) = X(i,j) * D(i): one multiply for each entry of Y.
      flopUpdate = numGlobalRows * numVectors;
    } else {
      // Y(i,j) = DampingFactor_ * X(i,j) * D(i):
      // Two multiplies per entry of Y.
      flopUpdate = 2.0 * numGlobalRows * numVectors;
    }
    ApplyFlops_ += flopUpdate;
    if (NumSweeps_ == 1) {
      return;
    }
  }
  // If we were allowed to assume that the starting guess was zero,
  // then we have already done the first sweep above.
  const int startSweep = ZeroStartingSolution_ ? 1 : 0;
  // We don't need to initialize the result MV, since the sparse
  // mat-vec will clobber its contents anyway.
  MV A_times_Y (Y.getMap (), as<size_t>(numVectors), false);
  for (int j = startSweep; j < NumSweeps_; ++j) {
    // Each iteration: Y = Y + \omega D^{-1} (X - A*Y)
    applyMat (Y, A_times_Y);
    A_times_Y.update (STS::one (), X, -STS::one ());
    Y.elementWiseMultiply (DampingFactor_, *Diagonal_, A_times_Y, STS::one ());
  }

  // For each column of output, for each pass over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const double numGlobalNonzeros = as<double> (A_->getGlobalNumEntries ());
  const double dampingFlops = (DampingFactor_ == STS::one ()) ? 0.0 : 1.0;
  ApplyFlops_ += as<double> (NumSweeps_ - startSweep) * numVectors *
    (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseJacobi_BlockCrsMatrix (const Tpetra::MultiVector<scalar_type, local_ordinal_type,
                                                             global_ordinal_type, node_type>& X,
                                   Tpetra::MultiVector<scalar_type, local_ordinal_type,
                                                       global_ordinal_type,node_type>& Y) const
{
  typedef Tpetra::Experimental::BlockMultiVector<scalar_type,
    local_ordinal_type, global_ordinal_type,node_type> BMV;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> MV;
  typedef typename block_crs_matrix_type::little_block_type little_block_type;
  typedef typename block_crs_matrix_type::little_vec_type little_vec_type;

  if (ZeroStartingSolution_) {
    Y.putScalar (STS::zero ());
  }

  const int numVectors = X.getNumVectors ();

  const block_crs_matrix_type* blockMat =
    dynamic_cast<const block_crs_matrix_type*> (A_.getRawPtr ());

  const local_ordinal_type blockSize = blockMat->getBlockSize();
  BMV A_times_Y_Block (*blockMat->getRowMap(), *Y.getMap (), blockSize, numVectors);
  MV A_times_Y = A_times_Y_Block.getMultiVectorView();
  BMV yBlock(Y, *blockMat->getRowMap(), blockSize);
  for (int j = 0; j < NumSweeps_; ++j) {
    blockMat->apply (Y, A_times_Y);
    A_times_Y.update (STS::one (), X, -STS::one ());
    const size_t numRows = blockMat->getNodeNumRows ();
    for (int i = 0; i < numVectors; ++i) {
      for (size_t k = 0; k < numRows; ++k) {
        // Each iteration: Y = Y + \omega D^{-1} (X - A*Y)
        little_block_type factorizedBlockDiag = BlockDiagonal_->getLocalBlock (k, k);
        little_vec_type xloc = A_times_Y_Block.getLocalBlock (k, i);
        little_vec_type yloc = yBlock.getLocalBlock (k, i);
        factorizedBlockDiag.solve (xloc, &blockDiagonalFactorizationPivots[i*blockSize]);
        yloc.update (DampingFactor_, xloc);
      }
    }
  }
}

template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseGS (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  typedef Relaxation<MatrixType> this_type;
  // The CrsMatrix version is faster, because it can access the sparse
  // matrix data directly, rather than by copying out each row's data
  // in turn.  Thus, we check whether the RowMatrix is really a
  // CrsMatrix.
  //
  // FIXME (mfh 07 Jul 2013) See note on crs_matrix_type typedef
  // declaration in Ifpack2_Relaxation_decl.hpp header file.  The code
  // will still be correct if the cast fails, but it will use an
  // unoptimized kernel.
  const block_crs_matrix_type* blockCrsMat =
    dynamic_cast<const block_crs_matrix_type*> (A_.getRawPtr ());
  const crs_matrix_type* crsMat =
    dynamic_cast<const crs_matrix_type*> (A_.getRawPtr ());
  if (blockCrsMat != NULL)  {
    const_cast<this_type*> (this)->ApplyInverseGS_BlockCrsMatrix (*blockCrsMat, X, Y);
  } else if (crsMat != NULL) {
    ApplyInverseGS_CrsMatrix (*crsMat, X, Y);
  } else {
    ApplyInverseGS_RowMatrix (X, Y);
  }
}


template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;



  // Tpetra's GS implementation for CrsMatrix handles zeroing out the
  // starting multivector itself.  The generic RowMatrix version here
  // does not, so we have to zero out Y here.
  if (ZeroStartingSolution_) {
    Y.putScalar (STS::zero ());
  }

  const size_t NumVectors = X.getNumVectors ();
  const size_t maxLength = A_->getNodeMaxNumRowEntries ();
  Array<local_ordinal_type> Indices (maxLength);
  Array<scalar_type> Values (maxLength);

  // Local smoothing stuff
  const size_t numMyRows             = A_->getNodeNumRows();
  const local_ordinal_type * rowInd  = 0;
  size_t numActive                   = numMyRows;
  bool do_local = !localSmoothingIndices_.is_null();
  if(do_local) {
    rowInd    = localSmoothingIndices_.getRawPtr();
    numActive = localSmoothingIndices_.size();
  }

  RCP<MV> Y2;
  if (IsParallel_) {
    if (Importer_.is_null ()) { // domain and column Maps are the same.
      // We will copy Y into Y2 below, so no need to fill with zeros here.
      Y2 = rcp (new MV (Y.getMap (), NumVectors, false));
    } else {
      // FIXME (mfh 21 Mar 2013) We probably don't need to fill with
      // zeros here, since we are doing an Import into Y2 below
      // anyway.  However, it doesn't hurt correctness.
      Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
    }
  }
  else {
    Y2 = rcpFromRef (Y);
  }

  // Diagonal
  ArrayRCP<const scalar_type>  d_rcp = Diagonal_->get1dView();
  ArrayView<const scalar_type> d_ptr = d_rcp();

  // Constant stride check
  bool constant_stride = X.isConstantStride() && Y2->isConstantStride();

  if(constant_stride) {
    // extract 1D RCPs
    size_t                    x_stride = X.getStride();
    size_t                   y2_stride = Y2->getStride();
    ArrayRCP<scalar_type>       y2_rcp = Y2->get1dViewNonConst();
    ArrayRCP<const scalar_type>  x_rcp = X.get1dView();
    ArrayView<scalar_type>      y2_ptr = y2_rcp();
    ArrayView<const scalar_type> x_ptr = x_rcp();
    Array<scalar_type> dtemp(NumVectors,STS::zero());

    for (int j = 0; j < NumSweeps_; j++) {
      // data exchange is here, once per sweep
      if (IsParallel_) {
        if (Importer_.is_null ()) {
          *Y2 = Y; // just copy, since domain and column Maps are the same
        } else {
          Y2->doImport (Y, *Importer_, Tpetra::INSERT);
        }
      }

      if (! DoBackwardGS_) { // Forward sweep
        for (size_t ii = 0; ii < numActive; ++ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);
          dtemp.assign(NumVectors,STS::zero());

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            for (size_t m = 0; m < NumVectors; ++m)
              dtemp[m] += Values[k] * y2_ptr[col + y2_stride*m];
          }

          for (size_t m = 0; m < NumVectors; ++m)
            y2_ptr[i + y2_stride*m] += DampingFactor_ * d_ptr[i] * (x_ptr[i + x_stride*m] - dtemp[m]);
        }
      }
      else { // Backward sweep
        // ptrdiff_t is the same size as size_t, but is signed.  Being
        // signed is important so that i >= 0 is not trivially true.
        for (ptrdiff_t ii = as<ptrdiff_t> (numActive) - 1; ii >= 0; --ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);
          dtemp.assign(NumVectors,STS::zero());

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            for (size_t m = 0; m < NumVectors; ++m)
              dtemp[m] += Values[k] * y2_ptr[col + y2_stride*m];
          }

          for (size_t m = 0; m < NumVectors; ++m)
            y2_ptr[i + y2_stride*m] += DampingFactor_ * d_ptr[i] * (x_ptr[i + x_stride*m] - dtemp[m]);
        }
      }
      // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
      if (IsParallel_) {
        Y = *Y2;
      }
    }
  }
  else {
    // extract 2D RCPS
    ArrayRCP<ArrayRCP<scalar_type> > y2_ptr      = Y2->get2dViewNonConst();
    ArrayRCP<ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();

    for (int j = 0; j < NumSweeps_; j++) {
      // data exchange is here, once per sweep
      if (IsParallel_) {
        if (Importer_.is_null ()) {
          *Y2 = Y; // just copy, since domain and column Maps are the same
        } else {
          Y2->doImport (Y, *Importer_, Tpetra::INSERT);
        }
      }

      if (! DoBackwardGS_) { // Forward sweep
        for (size_t ii = 0; ii < numActive; ++ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);

          for (size_t m = 0; m < NumVectors; ++m) {
            scalar_type dtemp = STS::zero ();
            ArrayView<const scalar_type> x_local = (x_ptr())[m]();
            ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();

            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              dtemp += Values[k] * y2_local[col];
            }
            y2_local[i] += DampingFactor_ * d_ptr[i] * (x_local[i] - dtemp);
          }
        }
      }
      else { // Backward sweep
        // ptrdiff_t is the same size as size_t, but is signed.  Being
        // signed is important so that i >= 0 is not trivially true.
        for (ptrdiff_t ii = as<ptrdiff_t> (numActive) - 1; ii >= 0; --ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);

          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);

          for (size_t m = 0; m < NumVectors; ++m) {
            scalar_type dtemp = STS::zero ();
            ArrayView<const scalar_type> x_local = (x_ptr())[m]();
            ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();

            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              dtemp += Values[k] * y2_local[col];
            }
            y2_local[i] += DampingFactor_ * d_ptr[i] * (x_local[i] - dtemp);
          }
        }
      }

      // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
      if (IsParallel_) {
        Y = *Y2;
      }
    }
  }


  // See flop count discussion in implementation of ApplyInverseGS_CrsMatrix().
  const double dampingFlops = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors = as<double> (X.getNumVectors ());
  const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
  const double numGlobalNonzeros = as<double> (A_->getGlobalNumEntries ());
  ApplyFlops_ += 2.0 * NumSweeps_ * numVectors *
    (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}


template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseGS_CrsMatrix (const crs_matrix_type& A,
                          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::as;
  const Tpetra::ESweepDirection direction =
    DoBackwardGS_ ? Tpetra::Backward : Tpetra::Forward;
  if(localSmoothingIndices_.is_null())
    A.gaussSeidelCopy (Y, X, *Diagonal_, DampingFactor_, direction,
                       NumSweeps_, ZeroStartingSolution_);
  else
    A.reorderedGaussSeidelCopy (Y, X, *Diagonal_, localSmoothingIndices_(), DampingFactor_, direction,
                                NumSweeps_, ZeroStartingSolution_);

  // For each column of output, for each sweep over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const double dampingFlops = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors = as<double> (X.getNumVectors ());
  const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
  const double numGlobalNonzeros = as<double> (A_->getGlobalNumEntries ());
  ApplyFlops_ += NumSweeps_ * numVectors *
    (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseGS_BlockCrsMatrix (const block_crs_matrix_type& A,
                          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y)
{

  typedef Tpetra::Experimental::BlockMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> BMV;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  //FIXME: (tcf) 8/21/2014 -- may be problematic for multiple right hand sides
  BMV yBlock(Y, *A.getGraph()->getDomainMap(), A.getBlockSize());
  const BMV xBlock(X, *A.getColMap(), A.getBlockSize());

  bool performImport = false;
  Teuchos::RCP<BMV> yBlockCol;
  if (Importer_.is_null())
    yBlockCol = Teuchos::rcpFromRef(yBlock);
  else
  {
    if (yBlockColumnPointMap_.is_null () ||
        yBlockColumnPointMap_->getNumVectors () != yBlock.getNumVectors () ||
        yBlockColumnPointMap_->getBlockSize () != yBlock.getBlockSize ())
    {
      yBlockColumnPointMap_ = Teuchos::rcp (new BMV (*A.getColMap(), A.getBlockSize (),
          static_cast<local_ordinal_type> (yBlock.getNumVectors ())));
    }
    yBlockCol = yBlockColumnPointMap_;
    performImport = true;
  }

  if (ZeroStartingSolution_) {
    yBlockCol->putScalar (STS::zero ());
  }
  else
  {
    yBlockCol->doImport(yBlock, *Importer_, Tpetra::INSERT);
  }

  const Tpetra::ESweepDirection direction =
    DoBackwardGS_ ? Tpetra::Backward : Tpetra::Forward;

  for (int sweep = 0; sweep < NumSweeps_; ++sweep)
  {
    if (performImport && sweep > 0) yBlockCol->doImport(yBlock, *Importer_, Tpetra::INSERT);
    A.localGaussSeidel(xBlock, *yBlockCol, *BlockDiagonal_, &blockDiagonalFactorizationPivots[0], DampingFactor_, direction);
    if (performImport)
    {
      Teuchos::RCP<const MV> yBlockColPointDomain = yBlockCol->getMultiVectorView().offsetView(A.getDomainMap(), 0);
      Tpetra::deep_copy(Y, *yBlockColPointDomain);
    }
  }


}


template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseSGS (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                 Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  typedef Relaxation<MatrixType> this_type;
  // The CrsMatrix version is faster, because it can access the sparse
  // matrix data directly, rather than by copying out each row's data
  // in turn.  Thus, we check whether the RowMatrix is really a
  // CrsMatrix.
  //
  // FIXME (mfh 07 Jul 2013) See note on crs_matrix_type typedef
  // declaration in Ifpack2_Relaxation_decl.hpp header file.  The code
  // will still be correct if the cast fails, but it will use an
  // unoptimized kernel.
  const block_crs_matrix_type* blockCrsMat = dynamic_cast<const block_crs_matrix_type*> (A_.getRawPtr());
  const crs_matrix_type* crsMat = dynamic_cast<const crs_matrix_type*> (&(*A_));
  if (blockCrsMat != NULL)  {
    const_cast<this_type*> (this)->ApplyInverseSGS_BlockCrsMatrix(*blockCrsMat, X, Y);
  }
  else if (crsMat != NULL) {
    ApplyInverseSGS_CrsMatrix (*crsMat, X, Y);
  } else {
    ApplyInverseSGS_RowMatrix (X, Y);
  }
}


template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseSGS_RowMatrix (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;

  // Tpetra's GS implementation for CrsMatrix handles zeroing out the
  // starting multivector itself.  The generic RowMatrix version here
  // does not, so we have to zero out Y here.
  if (ZeroStartingSolution_) {
    Y.putScalar (STS::zero ());
  }

  const size_t NumVectors = X.getNumVectors ();
  const size_t maxLength = A_->getNodeMaxNumRowEntries ();
  Array<local_ordinal_type> Indices (maxLength);
  Array<scalar_type> Values (maxLength);

  // Local smoothing stuff
  const size_t numMyRows             = A_->getNodeNumRows();
  const local_ordinal_type * rowInd  = 0;
  size_t numActive                   = numMyRows;
  bool do_local = localSmoothingIndices_.is_null();
  if(do_local) {
    rowInd    = localSmoothingIndices_.getRawPtr();
    numActive = localSmoothingIndices_.size();
  }


  RCP<MV> Y2;
  if (IsParallel_) {
    if (Importer_.is_null ()) { // domain and column Maps are the same.
      // We will copy Y into Y2 below, so no need to fill with zeros here.
      Y2 = rcp (new MV (Y.getMap (), NumVectors, false));
    } else {
      // FIXME (mfh 21 Mar 2013) We probably don't need to fill with
      // zeros here, since we are doing an Import into Y2 below
      // anyway.  However, it doesn't hurt correctness.
      Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
    }
  }
  else {
    Y2 = rcpFromRef (Y);
  }

  // Diagonal
  ArrayRCP<const scalar_type>  d_rcp = Diagonal_->get1dView();
  ArrayView<const scalar_type> d_ptr = d_rcp();

  // Constant stride check
  bool constant_stride = X.isConstantStride() && Y2->isConstantStride();

  if(constant_stride) {
    // extract 1D RCPs
    size_t                    x_stride = X.getStride();
    size_t                   y2_stride = Y2->getStride();
    ArrayRCP<scalar_type>       y2_rcp = Y2->get1dViewNonConst();
    ArrayRCP<const scalar_type>  x_rcp = X.get1dView();
    ArrayView<scalar_type>      y2_ptr = y2_rcp();
    ArrayView<const scalar_type> x_ptr = x_rcp();
    Array<scalar_type> dtemp(NumVectors,STS::zero());
    for (int iter = 0; iter < NumSweeps_; ++iter) {
      // only one data exchange per sweep
      if (IsParallel_) {
        if (Importer_.is_null ()) {
          *Y2 = Y; // just copy, since domain and column Maps are the same
        } else {
          Y2->doImport (Y, *Importer_, Tpetra::INSERT);
        }
      }

      for (int j = 0; j < NumSweeps_; j++) {
        // data exchange is here, once per sweep
        if (IsParallel_) {
          if (Importer_.is_null ()) {
            *Y2 = Y; // just copy, since domain and column Maps are the same
          } else {
            Y2->doImport (Y, *Importer_, Tpetra::INSERT);
          }
        }

        for (size_t ii = 0; ii < numActive; ++ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);
          dtemp.assign(NumVectors,STS::zero());

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            for (size_t m = 0; m < NumVectors; ++m)
              dtemp[m] += Values[k] * y2_ptr[col + y2_stride*m];
          }

          for (size_t m = 0; m < NumVectors; ++m)
            y2_ptr[i + y2_stride*m] += DampingFactor_ * d_ptr[i] * (x_ptr[i + x_stride*m] - dtemp[m]);
        }

        // ptrdiff_t is the same size as size_t, but is signed.  Being
        // signed is important so that i >= 0 is not trivially true.
        for (ptrdiff_t ii = as<ptrdiff_t> (numActive) - 1; ii >= 0; --ii) {
          local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
          size_t NumEntries;
          A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);
          dtemp.assign(NumVectors,STS::zero());

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            for (size_t m = 0; m < NumVectors; ++m)
              dtemp[m] += Values[k] * y2_ptr[col + y2_stride*m];
          }

          for (size_t m = 0; m < NumVectors; ++m)
            y2_ptr[i + y2_stride*m] += DampingFactor_ * d_ptr[i] * (x_ptr[i + x_stride*m] - dtemp[m]);
        }
      }

      // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
      if (IsParallel_) {
        Y = *Y2;
      }
    }
  }
  else {
    // extract 2D RCPs
    ArrayRCP<ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst ();
    ArrayRCP<ArrayRCP<const scalar_type> > x_ptr =  X.get2dView ();

    for (int iter = 0; iter < NumSweeps_; ++iter) {
      // only one data exchange per sweep
      if (IsParallel_) {
        if (Importer_.is_null ()) {
          *Y2 = Y; // just copy, since domain and column Maps are the same
        } else {
          Y2->doImport (Y, *Importer_, Tpetra::INSERT);
        }
      }

      for (size_t ii = 0; ii < numActive; ++ii) {
        local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
        const scalar_type diag = d_ptr[i];
        size_t NumEntries;
        A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

        for (size_t m = 0; m < NumVectors; ++m) {
          scalar_type dtemp = STS::zero ();
          ArrayView<const scalar_type> x_local = (x_ptr())[m]();
          ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_local[col];
          }
          y2_local[i] += DampingFactor_ * (x_local[i] - dtemp) * diag;
        }
      }

      // ptrdiff_t is the same size as size_t, but is signed.  Being
      // signed is important so that i >= 0 is not trivially true.
      for (ptrdiff_t ii = as<ptrdiff_t> (numActive) - 1; ii >= 0; --ii) {
        local_ordinal_type i = as<local_ordinal_type>(do_local ? rowInd[ii] : ii);
        const scalar_type diag = d_ptr[i];
        size_t NumEntries;
        A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

        for (size_t m = 0; m < NumVectors; ++m) {
          scalar_type dtemp = STS::zero ();
          ArrayView<const scalar_type> x_local = (x_ptr())[m]();
          ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();

          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_local[col];
          }
          y2_local[i] += DampingFactor_ * (x_local[i] - dtemp) * diag;
        }
      }

      // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
      if (IsParallel_) {
        Y = *Y2;
      }
    }
  }

  // See flop count discussion in implementation of ApplyInverseSGS_CrsMatrix().
  const double dampingFlops = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors = as<double> (X.getNumVectors ());
  const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
  const double numGlobalNonzeros = as<double> (A_->getGlobalNumEntries ());
  ApplyFlops_ += 2.0 * NumSweeps_ * numVectors *
    (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}


template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseSGS_CrsMatrix (const crs_matrix_type& A,
                           const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::as;
  const Tpetra::ESweepDirection direction = Tpetra::Symmetric;
  if(localSmoothingIndices_.is_null())
    A.gaussSeidelCopy (Y, X, *Diagonal_, DampingFactor_, direction,
                       NumSweeps_, ZeroStartingSolution_);
  else
    A.reorderedGaussSeidelCopy (Y, X, *Diagonal_, localSmoothingIndices_(), DampingFactor_, direction,
                                NumSweeps_, ZeroStartingSolution_);

  // For each column of output, for each sweep over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)
  //
  // Each sweep of symmetric Gauss-Seidel / SOR counts as two sweeps,
  // one forward and one backward.

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const double dampingFlops = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors = as<double> (X.getNumVectors ());
  const double numGlobalRows = as<double> (A_->getGlobalNumRows ());
  const double numGlobalNonzeros = as<double> (A_->getGlobalNumEntries ());
  ApplyFlops_ += 2.0 * NumSweeps_ * numVectors *
    (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseSGS_BlockCrsMatrix (const block_crs_matrix_type& A,
                           const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y)
{
  typedef Tpetra::Experimental::BlockMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> BMV;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  //FIXME: (tcf) 8/21/2014 -- may be problematic for multiple right hand sides
  BMV yBlock(Y, *A.getGraph()->getDomainMap(), A.getBlockSize());
  const BMV xBlock(X, *A.getColMap(), A.getBlockSize());

  bool performImport = false;
  Teuchos::RCP<BMV> yBlockCol;
  if (Importer_.is_null())
    yBlockCol = Teuchos::rcpFromRef(yBlock);
  else
  {
    if (yBlockColumnPointMap_.is_null () ||
        yBlockColumnPointMap_->getNumVectors () != yBlock.getNumVectors () ||
        yBlockColumnPointMap_->getBlockSize () != yBlock.getBlockSize ())
    {
      yBlockColumnPointMap_ = Teuchos::rcp (new BMV (*A.getColMap(), A.getBlockSize (),
          static_cast<local_ordinal_type> (yBlock.getNumVectors ())));
    }
    yBlockCol = yBlockColumnPointMap_;
    performImport = true;
  }

  if (ZeroStartingSolution_) {
    yBlockCol->putScalar (STS::zero ());
  }
  else
  {
    yBlockCol->doImport(yBlock, *Importer_, Tpetra::INSERT);
  }


  const Tpetra::ESweepDirection direction = Tpetra::Symmetric;

  for (int sweep = 0; sweep < NumSweeps_; ++sweep)
  {
    if (performImport && sweep > 0) yBlockCol->doImport(yBlock, *Importer_, Tpetra::INSERT);
    A.localGaussSeidel(xBlock, *yBlockCol, *BlockDiagonal_, &blockDiagonalFactorizationPivots[0], DampingFactor_, direction);
    if (performImport)
    {
      Teuchos::RCP<const MV> yBlockColPointDomain = yBlockCol->getMultiVectorView().offsetView(A.getDomainMap(), 0);
      MV yBlockView = yBlock.getMultiVectorView();
      Tpetra::deep_copy(yBlockView, *yBlockColPointDomain);
    }
  }


}


template<class MatrixType>
std::string Relaxation<MatrixType>::description () const {
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Relaxation\": {";

  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  // It's useful to print this instance's relaxation method (Jacobi,
  // Gauss-Seidel, or symmetric Gauss-Seidel).  If you want more info
  // than that, call describe() instead.
  os << "Type: ";
  if (PrecType_ == Ifpack2::Details::JACOBI) {
    os << "Jacobi";
  } else if (PrecType_ == Ifpack2::Details::GS) {
    os << "Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::SGS) {
    os << "Symmetric Gauss-Seidel";
  } else {
    os << "INVALID";
  }

  os  << ", " << "sweeps: " << NumSweeps_ << ", "
      << "damping factor: " << DampingFactor_ << ", ";
  if (DoL1Method_) {
    os << "use l1: " << DoL1Method_ << ", "
       << "l1 eta: " << L1Eta_ << ", ";
  }

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


template<class MatrixType>
void
Relaxation<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  const int myRank = this->getComm ()->getRank ();

  //    none: print nothing
  //     low: print O(1) info from Proc 0
  //  medium:
  //    high:
  // extreme:

  if (vl != VERB_NONE && myRank == 0) {
    // Describable interface asks each implementation to start with a tab.
    OSTab tab1 (out);
    // Output is valid YAML; hence the quotes, to protect the colons.
    out << "\"Ifpack2::Relaxation\":" << endl;
    OSTab tab2 (out);
    out << "MatrixType: \"" << TypeNameTraits<MatrixType>::name () << "\""
        << endl;
    if (this->getObjectLabel () != "") {
      out << "Label: " << this->getObjectLabel () << endl;
    }
    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl
        << "Computed: " << (isComputed () ? "true" : "false") << endl
        << "Parameters: " << endl;
    {
      OSTab tab3 (out);
      out << "\"relaxation: type\": ";
      if (PrecType_ == Ifpack2::Details::JACOBI) {
        out << "Jacobi";
      } else if (PrecType_ == Ifpack2::Details::GS) {
        out << "Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::SGS) {
        out << "Symmetric Gauss-Seidel";
      } else {
        out << "INVALID";
      }
      // We quote these parameter names because they contain colons.
      // YAML uses the colon to distinguish key from value.
      out << endl
          << "\"relaxation: sweeps\": " << NumSweeps_ << endl
          << "\"relaxation: damping factor\": " << DampingFactor_ << endl
          << "\"relaxation: min diagonal value\": " << MinDiagonalValue_ << endl
          << "\"relaxation: zero starting solution\": " << ZeroStartingSolution_ << endl
          << "\"relaxation: backward mode\": " << DoBackwardGS_ << endl
          << "\"relaxation: use l1\": " << DoL1Method_ << endl
          << "\"relaxation: l1 eta\": " << L1Eta_ << endl;
    }
    out << "Computed quantities:" << endl;
    {
      OSTab tab3 (out);
      out << "Condition number estimate: " << Condest_ << endl
          << "Global number of rows: " << A_->getGlobalNumRows () << endl
          << "Global number of columns: " << A_->getGlobalNumCols () << endl;
    }
    if (checkDiagEntries_ && isComputed ()) {
      out << "Properties of input diagonal entries:" << endl;
      {
        OSTab tab3 (out);
        out << "Magnitude of minimum-magnitude diagonal entry: "
            << globalMinMagDiagEntryMag_ << endl
            << "Magnitude of maximum-magnitude diagonal entry: "
            << globalMaxMagDiagEntryMag_ << endl
            << "Number of diagonal entries with small magnitude: "
            << globalNumSmallDiagEntries_ << endl
            << "Number of zero diagonal entries: "
            << globalNumZeroDiagEntries_ << endl
            << "Number of diagonal entries with negative real part: "
            << globalNumNegDiagEntries_ << endl
            << "Abs 2-norm diff between computed and actual inverse "
            << "diagonal: " << globalDiagNormDiff_ << endl;
      }
    }
    if (isComputed ()) {
      out << "Saved diagonal offsets: "
          << (savedDiagOffsets_ ? "true" : "false") << endl;
    }
    out << "Call counts and total times (in seconds): " << endl;
    {
      OSTab tab3 (out);
      out << "initialize: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumInitialize_ << endl;
        out << "Total time: " << InitializeTime_ << endl;
      }
      out << "compute: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumCompute_ << endl;
        out << "Total time: " << ComputeTime_ << endl;
      }
      out << "apply: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumApply_ << endl;
        out << "Total time: " << ApplyTime_ << endl;
      }
    }
  }
}

} // namespace Ifpack2

#define IFPACK2_RELAXATION_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::Relaxation< Tpetra::CrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::Relaxation< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_RELAXATION_DEF_HPP
