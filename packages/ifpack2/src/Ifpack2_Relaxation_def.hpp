// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_RELAXATION_DEF_HPP
#define IFPACK2_RELAXATION_DEF_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockView.hpp"
#include "Ifpack2_Utilities.hpp"
#include "Ifpack2_Details_getCrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Details_residual.hpp"
#include <cstdlib>
#include <memory>
#include <sstream>
#include "KokkosSparse_gauss_seidel.hpp"

namespace {
// Validate that a given int is nonnegative.
class NonnegativeIntValidator : public Teuchos::ParameterEntryValidator {
 public:
  // Constructor (does nothing).
  NonnegativeIntValidator() {}

  // ParameterEntryValidator wants this method.
  Teuchos::ParameterEntryValidator::ValidStringsList validStringValues() const {
    return Teuchos::null;
  }

  // Actually validate the parameter's value.
  void
  validate(const Teuchos::ParameterEntry& entry,
           const std::string& paramName,
           const std::string& sublistName) const {
    using std::endl;
    Teuchos::any anyVal         = entry.getAny(true);
    const std::string entryName = entry.getAny(false).typeName();

    TEUCHOS_TEST_FOR_EXCEPTION(
        anyVal.type() != typeid(int),
        Teuchos::Exceptions::InvalidParameterType,
        "Parameter \"" << paramName << "\" in sublist \"" << sublistName
                       << "\" has the wrong type." << endl
                       << "Parameter: " << paramName
                       << endl
                       << "Type specified: " << entryName << endl
                       << "Type required: int" << endl);

    const int val = Teuchos::any_cast<int>(anyVal);
    TEUCHOS_TEST_FOR_EXCEPTION(
        val < 0, Teuchos::Exceptions::InvalidParameterValue,
        "Parameter \"" << paramName << "\" in sublist \"" << sublistName
                       << "\" is negative." << endl
                       << "Parameter: " << paramName
                       << endl
                       << "Value specified: " << val << endl
                       << "Required range: [0, INT_MAX]" << endl);
  }

  // ParameterEntryValidator wants this method.
  const std::string getXMLTypeName() const {
    return "NonnegativeIntValidator";
  }

  // ParameterEntryValidator wants this method.
  void
  printDoc(const std::string& docString,
           std::ostream& out) const {
    Teuchos::StrUtils::printLines(out, "# ", docString);
    out << "#\tValidator Used: " << std::endl;
    out << "#\t\tNonnegativeIntValidator" << std::endl;
  }
};

// A way to get a small positive number (eps() for floating-point
// types, or 1 for integer types) when Teuchos::ScalarTraits doesn't
// define it (for example, for integer values).
template <class Scalar, const bool isOrdinal = Teuchos::ScalarTraits<Scalar>::isOrdinal>
class SmallTraits {
 public:
  // Return eps if Scalar is a floating-point type, else return 1.
  static const Scalar eps();
};

// Partial specialization for when Scalar is not a floating-point type.
template <class Scalar>
class SmallTraits<Scalar, true> {
 public:
  static const Scalar eps() {
    return Teuchos::ScalarTraits<Scalar>::one();
  }
};

// Partial specialization for when Scalar is a floating-point type.
template <class Scalar>
class SmallTraits<Scalar, false> {
 public:
  static const Scalar eps() {
    return Teuchos::ScalarTraits<Scalar>::eps();
  }
};

// Work-around for GitHub Issue #5269.
template <class Scalar,
          const bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex>
struct RealTraits {};

template <class Scalar>
struct RealTraits<Scalar, false> {
  using val_type = Scalar;
  using mag_type = Scalar;
  static KOKKOS_INLINE_FUNCTION mag_type real(const val_type& z) {
    return z;
  }
};

template <class Scalar>
struct RealTraits<Scalar, true> {
  using val_type = Scalar;
  using mag_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  static KOKKOS_INLINE_FUNCTION mag_type real(const val_type& z) {
#if KOKKOS_VERSION >= 40799
    return KokkosKernels::ArithTraits<val_type>::real(z);
#else
    return Kokkos::ArithTraits<val_type>::real(z);
#endif
  }
};

template <class Scalar>
KOKKOS_INLINE_FUNCTION typename RealTraits<Scalar>::mag_type
getRealValue(const Scalar& z) {
  return RealTraits<Scalar>::real(z);
}

}  // namespace

namespace Ifpack2 {

template <class MatrixType>
void Relaxation<MatrixType>::
    updateCachedMultiVector(const Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>>& map,
                            size_t numVecs) const {
  // Allocate a multivector if the cached one isn't perfect.  Checking
  // for map pointer equality is much cheaper than Map::isSameAs.
  if (cachedMV_.is_null() ||
      map.get() != cachedMV_->getMap().get() ||
      cachedMV_->getNumVectors() != numVecs) {
    using MV  = Tpetra::MultiVector<scalar_type, local_ordinal_type,
                                   global_ordinal_type, node_type>;
    cachedMV_ = Teuchos::rcp(new MV(map, numVecs, false));
  }
}

template <class MatrixType>
void Relaxation<MatrixType>::
    setMatrix(const Teuchos::RCP<const row_matrix_type>& A) {
  if (A.getRawPtr() != A_.getRawPtr()) {  // it's a different matrix
    Importer_          = Teuchos::null;
    pointImporter_     = Teuchos::null;
    Diagonal_          = Teuchos::null;  // ??? what if this comes from the user???
    isInitialized_     = false;
    IsComputed_        = false;
    diagOffsets_       = Kokkos::View<size_t*, device_type>();
    savedDiagOffsets_  = false;
    hasBlockCrsMatrix_ = false;
    serialGaussSeidel_ = Teuchos::null;
    if (!A.is_null()) {
      IsParallel_ = (A->getRowMap()->getComm()->getSize() > 1);
    }
    A_ = A;
  }
}

template <class MatrixType>
Relaxation<MatrixType>::
    Relaxation(const Teuchos::RCP<const row_matrix_type>& A)
  : A_(A)
  , IsParallel_((A.is_null() || A->getRowMap().is_null() || A->getRowMap()->getComm().is_null()) ? false :  // a reasonable default if there's no communicator
                    A->getRowMap()->getComm()->getSize() > 1) {
  this->setObjectLabel("Ifpack2::Relaxation");
}

template <class MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
Relaxation<MatrixType>::getValidParameters() const {
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::setStringToIntegralParameter;
  typedef Teuchos::ParameterEntryValidator PEV;

  if (validParams_.is_null()) {
    RCP<ParameterList> pl = parameterList("Ifpack2::Relaxation");

    // Set a validator that automatically converts from the valid
    // string options to their enum values.
    Array<std::string> precTypes(8);
    precTypes[0] = "Jacobi";
    precTypes[1] = "Gauss-Seidel";
    precTypes[2] = "Symmetric Gauss-Seidel";
    precTypes[3] = "MT Gauss-Seidel";
    precTypes[4] = "MT Symmetric Gauss-Seidel";
    precTypes[5] = "Richardson";
    precTypes[6] = "Two-stage Gauss-Seidel";
    precTypes[7] = "Two-stage Symmetric Gauss-Seidel";
    Array<Details::RelaxationType> precTypeEnums(8);
    precTypeEnums[0] = Details::JACOBI;
    precTypeEnums[1] = Details::GS;
    precTypeEnums[2] = Details::SGS;
    precTypeEnums[3] = Details::MTGS;
    precTypeEnums[4] = Details::MTSGS;
    precTypeEnums[5] = Details::RICHARDSON;
    precTypeEnums[6] = Details::GS2;
    precTypeEnums[7] = Details::SGS2;
    const std::string defaultPrecType("Jacobi");
    setStringToIntegralParameter<Details::RelaxationType>("relaxation: type",
                                                          defaultPrecType, "Relaxation method", precTypes(), precTypeEnums(),
                                                          pl.getRawPtr());

    const int numSweeps = 1;
    RCP<PEV> numSweepsValidator =
        rcp_implicit_cast<PEV>(rcp(new NonnegativeIntValidator));
    pl->set("relaxation: sweeps", numSweeps, "Number of relaxation sweeps",
            rcp_const_cast<const PEV>(numSweepsValidator));

    // number of 'local' outer sweeps for two-stage GS
    const int numOuterSweeps = 1;
    RCP<PEV> numOuterSweepsValidator =
        rcp_implicit_cast<PEV>(rcp(new NonnegativeIntValidator));
    pl->set("relaxation: outer sweeps", numOuterSweeps, "Number of outer local relaxation sweeps for two-stage GS",
            rcp_const_cast<const PEV>(numOuterSweepsValidator));
    // number of 'local' inner sweeps for two-stage GS
    const int numInnerSweeps = 1;
    RCP<PEV> numInnerSweepsValidator =
        rcp_implicit_cast<PEV>(rcp(new NonnegativeIntValidator));
    pl->set("relaxation: inner sweeps", numInnerSweeps, "Number of inner local relaxation sweeps for two-stage GS",
            rcp_const_cast<const PEV>(numInnerSweepsValidator));
    // specify damping factor for the inner sweeps of two-stage GS
    const scalar_type innerDampingFactor = STS::one();
    pl->set("relaxation: inner damping factor", innerDampingFactor, "Damping factor for the inner sweep of two-stage GS");
    // specify whether to use sptrsv instead of inner-iterations for two-stage GS
    const bool innerSpTrsv = false;
    pl->set("relaxation: inner sparse-triangular solve", innerSpTrsv, "Specify whether to use sptrsv instead of JR iterations for two-stage GS");
    // specify whether to use compact form of recurrence for two-stage GS
    const bool compactForm = false;
    pl->set("relaxation: compact form", compactForm, "Specify whether to use compact form of recurrence for two-stage GS");

    const scalar_type dampingFactor = STS::one();
    pl->set("relaxation: damping factor", dampingFactor);

    const bool zeroStartingSolution = true;
    pl->set("relaxation: zero starting solution", zeroStartingSolution);

    const bool doBackwardGS = false;
    pl->set("relaxation: backward mode", doBackwardGS);

    const bool doL1Method = false;
    pl->set("relaxation: use l1", doL1Method);

    const magnitude_type l1eta = (STM::one() + STM::one() + STM::one()) /
                                 (STM::one() + STM::one());  // 1.5
    pl->set("relaxation: l1 eta", l1eta);

    const scalar_type minDiagonalValue = STS::zero();
    pl->set("relaxation: min diagonal value", minDiagonalValue);

    const bool fixTinyDiagEntries = false;
    pl->set("relaxation: fix tiny diagonal entries", fixTinyDiagEntries);

    const bool checkDiagEntries = false;
    pl->set("relaxation: check diagonal entries", checkDiagEntries);

    Teuchos::ArrayRCP<local_ordinal_type> localSmoothingIndices = Teuchos::null;
    pl->set("relaxation: local smoothing indices", localSmoothingIndices);

    const bool is_matrix_structurally_symmetric = false;
    pl->set("relaxation: symmetric matrix structure", is_matrix_structurally_symmetric);

    const bool ifpack2_dump_matrix = false;
    pl->set("relaxation: ifpack2 dump matrix", ifpack2_dump_matrix);

    const int cluster_size = 1;
    pl->set("relaxation: mtgs cluster size", cluster_size);

    pl->set("relaxation: mtgs coloring algorithm", "Default");

    const int long_row_threshold = 0;
    pl->set("relaxation: long row threshold", long_row_threshold);

    const bool timer_for_apply = true;
    pl->set("timer for apply", timer_for_apply);

    validParams_ = rcp_const_cast<const ParameterList>(pl);
  }
  return validParams_;
}

template <class MatrixType>
void Relaxation<MatrixType>::setParametersImpl(Teuchos::ParameterList& pl) {
  using Teuchos::getIntegralValue;
  using Teuchos::getStringValue;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  typedef scalar_type ST;  // just to make code below shorter

  if (pl.isType<double>("relaxation: damping factor")) {
    // Make sure that ST=complex can run with a damping factor that is
    // a double.
    ST df = pl.get<double>("relaxation: damping factor");
    pl.remove("relaxation: damping factor");
    pl.set("relaxation: damping factor", df);
  }

  pl.validateParametersAndSetDefaults(*getValidParameters());

  const Details::RelaxationType precType =
      getIntegralValue<Details::RelaxationType>(pl, "relaxation: type");
  const std::string precTypeStr = getStringValue<Details::RelaxationType>(pl, "relaxation: type");
  // We only access "relaxation: type" using strings in the rest of the code
  pl.set<std::string>("relaxation: type", precTypeStr);
  pl.get<std::string>("relaxation: type");  // We need to mark the parameter as "used"
  const int numSweeps                         = pl.get<int>("relaxation: sweeps");
  const ST dampingFactor                      = pl.get<ST>("relaxation: damping factor");
  const bool zeroStartSol                     = pl.get<bool>("relaxation: zero starting solution");
  const bool doBackwardGS                     = pl.get<bool>("relaxation: backward mode");
  const bool doL1Method                       = pl.get<bool>("relaxation: use l1");
  const magnitude_type l1Eta                  = pl.get<magnitude_type>("relaxation: l1 eta");
  const ST minDiagonalValue                   = pl.get<ST>("relaxation: min diagonal value");
  const bool fixTinyDiagEntries               = pl.get<bool>("relaxation: fix tiny diagonal entries");
  const bool checkDiagEntries                 = pl.get<bool>("relaxation: check diagonal entries");
  const bool is_matrix_structurally_symmetric = pl.get<bool>("relaxation: symmetric matrix structure");
  const bool ifpack2_dump_matrix              = pl.get<bool>("relaxation: ifpack2 dump matrix");
  const bool timer_for_apply                  = pl.get<bool>("timer for apply");
  int cluster_size                            = 1;
  if (pl.isParameter("relaxation: mtgs cluster size"))  // optional parameter
    cluster_size = pl.get<int>("relaxation: mtgs cluster size");
  int long_row_threshold = 0;
  if (pl.isParameter("relaxation: long row threshold"))  // optional parameter
    long_row_threshold = pl.get<int>("relaxation: long row threshold");
  std::string color_algo_name = pl.get<std::string>("relaxation: mtgs coloring algorithm");
  // convert to lowercase
  for (char& c : color_algo_name)
    c = tolower(c);
  if (color_algo_name == "default")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_DEFAULT;
  else if (color_algo_name == "serial")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_SERIAL;
  else if (color_algo_name == "vb")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_VB;
  else if (color_algo_name == "vbbit")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_VBBIT;
  else if (color_algo_name == "vbcs")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_VBCS;
  else if (color_algo_name == "vbd")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_VBD;
  else if (color_algo_name == "vbdbit")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_VBDBIT;
  else if (color_algo_name == "eb")
    this->mtColoringAlgorithm_ = KokkosGraph::COLORING_EB;
  else {
    std::ostringstream msg;
    msg << "Ifpack2::Relaxation: 'relaxation: mtgs coloring algorithm' = '" << color_algo_name << "' is not valid.\n";
    msg << "Choices (not case sensitive) are: Default, Serial, VB, VBBIT, VBCS, VBD, VBDBIT, EB.";
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, msg.str());
  }

  Teuchos::ArrayRCP<local_ordinal_type> localSmoothingIndices = pl.get<Teuchos::ArrayRCP<local_ordinal_type>>("relaxation: local smoothing indices");

  // for Two-stage Gauss-Seidel
  if (!std::is_same<double, ST>::value && pl.isType<double>("relaxation: inner damping factor")) {
    // Make sure that ST=complex can run with a damping factor that is
    // a double.
    ST df = pl.get<double>("relaxation: inner damping factor");
    pl.remove("relaxation: inner damping factor");
    pl.set("relaxation: inner damping factor", df);
  }
  // If long row algorithm was requested, make sure non-cluster (point) multicolor Gauss-Seidel (aka MTGS/MTSGS) will be used.
  if (long_row_threshold > 0) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        cluster_size != 1, std::invalid_argument,
        "Ifpack2::Relaxation: "
        "Requested long row MTGS/MTSGS algorithm and cluster GS/SGS, but those are not compatible.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        precType != Details::RelaxationType::MTGS && precType != Details::RelaxationType::MTSGS,
        std::invalid_argument,
        "Ifpack2::Relaxation: "
        "Requested long row MTGS/MTSGS algorithm, but this is only compatible with preconditioner types "
        "'MT Gauss-Seidel' and 'MT Symmetric Gauss-Seidel'.");
  }

  const ST innerDampingFactor = pl.get<ST>("relaxation: inner damping factor");
  const int numInnerSweeps    = pl.get<int>("relaxation: inner sweeps");
  const int numOuterSweeps    = pl.get<int>("relaxation: outer sweeps");
  const bool innerSpTrsv      = pl.get<bool>("relaxation: inner sparse-triangular solve");
  const bool compactForm      = pl.get<bool>("relaxation: compact form");

  // "Commit" the changes, now that we've validated everything.
  PrecType_                         = precType;
  NumSweeps_                        = numSweeps;
  DampingFactor_                    = dampingFactor;
  ZeroStartingSolution_             = zeroStartSol;
  DoBackwardGS_                     = doBackwardGS;
  DoL1Method_                       = doL1Method;
  L1Eta_                            = l1Eta;
  MinDiagonalValue_                 = minDiagonalValue;
  fixTinyDiagEntries_               = fixTinyDiagEntries;
  checkDiagEntries_                 = checkDiagEntries;
  clusterSize_                      = cluster_size;
  longRowThreshold_                 = long_row_threshold;
  is_matrix_structurally_symmetric_ = is_matrix_structurally_symmetric;
  ifpack2_dump_matrix_              = ifpack2_dump_matrix;
  localSmoothingIndices_            = localSmoothingIndices;
  // for Two-stage GS
  NumInnerSweeps_     = numInnerSweeps;
  NumOuterSweeps_     = numOuterSweeps;
  InnerSpTrsv_        = innerSpTrsv;
  InnerDampingFactor_ = innerDampingFactor;
  CompactForm_        = compactForm;
  TimerForApply_      = timer_for_apply;
}

template <class MatrixType>
void Relaxation<MatrixType>::setParameters(const Teuchos::ParameterList& pl) {
  // FIXME (aprokop 18 Oct 2013) Casting away const is bad here.
  // but otherwise, we will get [unused] in pl
  this->setParametersImpl(const_cast<Teuchos::ParameterList&>(pl));
}

template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int>>
Relaxation<MatrixType>::getComm() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::getComm: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix before calling this method.");
  return A_->getRowMap()->getComm();
}

template <class MatrixType>
Teuchos::RCP<const typename Relaxation<MatrixType>::row_matrix_type>
Relaxation<MatrixType>::getMatrix() const {
  return A_;
}

template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type>>
Relaxation<MatrixType>::getDomainMap() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::getDomainMap: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix before calling this method.");
  return A_->getDomainMap();
}

template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type>>
Relaxation<MatrixType>::getRangeMap() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::getRangeMap: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix before calling this method.");
  return A_->getRangeMap();
}

template <class MatrixType>
bool Relaxation<MatrixType>::hasTransposeApply() const {
  return true;
}

template <class MatrixType>
int Relaxation<MatrixType>::getNumInitialize() const {
  return (NumInitialize_);
}

template <class MatrixType>
int Relaxation<MatrixType>::getNumCompute() const {
  return (NumCompute_);
}

template <class MatrixType>
int Relaxation<MatrixType>::getNumApply() const {
  return (NumApply_);
}

template <class MatrixType>
double Relaxation<MatrixType>::getInitializeTime() const {
  return (InitializeTime_);
}

template <class MatrixType>
double Relaxation<MatrixType>::getComputeTime() const {
  return (ComputeTime_);
}

template <class MatrixType>
double Relaxation<MatrixType>::getApplyTime() const {
  return (ApplyTime_);
}

template <class MatrixType>
double Relaxation<MatrixType>::getComputeFlops() const {
  return (ComputeFlops_);
}

template <class MatrixType>
double Relaxation<MatrixType>::getApplyFlops() const {
  return (ApplyFlops_);
}

template <class MatrixType>
size_t Relaxation<MatrixType>::getNodeSmootherComplexity() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::getNodeSmootherComplexity: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix, then call compute(), before calling this method.");
  // Relaxation methods cost roughly one apply + one diagonal inverse per iteration
  return A_->getLocalNumRows() + A_->getLocalNumEntries();
}

template <class MatrixType>
void Relaxation<MatrixType>::
    apply(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
          Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
          Teuchos::ETransp /* mode */,
          scalar_type alpha,
          scalar_type beta) const {
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type>
      MV;
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::apply: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix, then call compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !isComputed(),
      std::runtime_error,
      "Ifpack2::Relaxation::apply: You must call compute() on this Ifpack2 "
      "preconditioner instance before you may call apply().  You may call "
      "isComputed() to find out if compute() has been called already.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(),
      std::runtime_error,
      "Ifpack2::Relaxation::apply: X and Y have different numbers of columns.  "
      "X has "
          << X.getNumVectors() << " columns, but Y has "
          << Y.getNumVectors() << " columns.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      beta != STS::zero(), std::logic_error,
      "Ifpack2::Relaxation::apply: beta = " << beta << " != 0 case not "
                                                       "implemented.");

  Teuchos::RCP<Teuchos::Time> timer;
  const std::string timerName("Ifpack2::Relaxation::apply");
  if (TimerForApply_) {
    timer = Teuchos::TimeMonitor::lookupCounter(timerName);
    if (timer.is_null()) {
      timer = Teuchos::TimeMonitor::getNewCounter(timerName);
    }
  }

  Teuchos::Time time = Teuchos::Time(timerName);
  double startTime   = time.wallTime();

  {
    Teuchos::RCP<Teuchos::TimeMonitor> timeMon;
    if (TimerForApply_)
      timeMon = Teuchos::rcp(new Teuchos::TimeMonitor(*timer));

    // Special case: alpha == 0.
    if (alpha == STS::zero()) {
      // No floating-point operations, so no need to update a count.
      Y.putScalar(STS::zero());
    } else {
      // If X and Y alias one another, then we need to create an
      // auxiliary vector, Xcopy (a deep copy of X).
      RCP<const MV> Xcopy;
      {
        if (X.aliases(Y)) {
          Xcopy = rcp(new MV(X, Teuchos::Copy));
        } else {
          Xcopy = rcpFromRef(X);
        }
      }

      // Each of the following methods updates the flop count itself.
      // All implementations handle zeroing out the starting solution
      // (if necessary) themselves.
      switch (PrecType_) {
        case Ifpack2::Details::JACOBI:
          ApplyInverseJacobi(*Xcopy, Y);
          break;
        case Ifpack2::Details::GS:
          ApplyInverseSerialGS(*Xcopy, Y, DoBackwardGS_ ? Tpetra::Backward : Tpetra::Forward);
          break;
        case Ifpack2::Details::SGS:
          ApplyInverseSerialGS(*Xcopy, Y, Tpetra::Symmetric);
          break;
        case Ifpack2::Details::MTGS:
        case Ifpack2::Details::GS2:
          ApplyInverseMTGS_CrsMatrix(*Xcopy, Y, DoBackwardGS_ ? Tpetra::Backward : Tpetra::Forward);
          break;
        case Ifpack2::Details::MTSGS:
        case Ifpack2::Details::SGS2:
          ApplyInverseMTGS_CrsMatrix(*Xcopy, Y, Tpetra::Symmetric);
          break;
        case Ifpack2::Details::RICHARDSON:
          ApplyInverseRichardson(*Xcopy, Y);
          break;

        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                     "Ifpack2::Relaxation::apply: Invalid preconditioner type enum value "
                                         << PrecType_ << ".  Please report this bug to the Ifpack2 developers.");
      }
      if (alpha != STS::one()) {
        Y.scale(alpha);
        const double numGlobalRows = as<double>(A_->getGlobalNumRows());
        const double numVectors    = as<double>(Y.getNumVectors());
        ApplyFlops_ += numGlobalRows * numVectors;
      }
    }
  }
  ApplyTime_ += (time.wallTime() - startTime);
  ++NumApply_;
}

template <class MatrixType>
void Relaxation<MatrixType>::
    applyMat(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
             Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
             Teuchos::ETransp mode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error,
      "Ifpack2::Relaxation::applyMat: "
      "The input matrix A is null.  Please call setMatrix() with a nonnull "
      "input matrix, then call compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !isComputed(), std::runtime_error,
      "Ifpack2::Relaxation::applyMat: "
      "isComputed() must be true before you may call applyMat().  "
      "Please call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
      "Ifpack2::Relaxation::applyMat: X.getNumVectors() = " << X.getNumVectors()
                                                            << " != Y.getNumVectors() = " << Y.getNumVectors() << ".");
  A_->apply(X, Y, mode);
}

template <class MatrixType>
void Relaxation<MatrixType>::initialize() {
  const char methodName[] = "Ifpack2::Relaxation::initialize";

  TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null(), std::runtime_error, methodName << ": The "
                                                                             "input matrix A is null.  Please call setMatrix() with "
                                                                             "a nonnull input matrix before calling this method.");

  Teuchos::RCP<Teuchos::Time> timer =
      Teuchos::TimeMonitor::getNewCounter(methodName);

  double startTime = timer->wallTime();

  {  // Timing of initialize starts here
    Teuchos::TimeMonitor timeMon(*timer);
    isInitialized_ = false;

    {
      auto rowMap = A_->getRowMap();
      auto comm   = rowMap.is_null() ? Teuchos::null : rowMap->getComm();
      IsParallel_ = !comm.is_null() && comm->getSize() > 1;
    }

    // mfh 21 Mar 2013, 07 May 2019: The Import object may be null,
    // but in that case, the domain and column Maps are the same and
    // we don't need to Import anyway.
    //
    // mfh 07 May 2019: A_->getGraph() might be an
    // OverlappingRowGraph, which doesn't have an Import object.
    // However, in that case, the comm should have size 1.
    Importer_ = IsParallel_ ? A_->getGraph()->getImporter() : Teuchos::null;

    {
      Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
          Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(A_);
      hasBlockCrsMatrix_ = !A_bcrs.is_null();
    }

    serialGaussSeidel_ = Teuchos::null;

    if (PrecType_ == Details::MTGS || PrecType_ == Details::MTSGS ||
        PrecType_ == Details::GS2 || PrecType_ == Details::SGS2) {
      auto crsMat = Details::getCrsMatrix(A_);
      TEUCHOS_TEST_FOR_EXCEPTION(crsMat.is_null(), std::logic_error, methodName << ": "
                                                                                   "Multithreaded Gauss-Seidel methods currently only work "
                                                                                   "when the input matrix is a Tpetra::CrsMatrix.");

      // FIXME (mfh 27 May 2019) Dumping the matrix belongs in
      // compute, not initialize, since users may change the matrix's
      // values at any time before calling compute.
      if (ifpack2_dump_matrix_) {
        static int sequence_number  = 0;
        const std::string file_name = "Ifpack2_MT_GS_" +
                                      std::to_string(sequence_number++) + ".mtx";
        using writer_type = Tpetra::MatrixMarket::Writer<crs_matrix_type>;
        writer_type::writeSparseFile(file_name, crsMat);
      }

      this->mtKernelHandle_ = Teuchos::rcp(new mt_kernel_handle_type());
      if (mtKernelHandle_->get_gs_handle() == nullptr) {
        if (PrecType_ == Details::GS2 || PrecType_ == Details::SGS2)
          mtKernelHandle_->create_gs_handle(KokkosSparse::GS_TWOSTAGE);
        else if (this->clusterSize_ == 1) {
          mtKernelHandle_->create_gs_handle(KokkosSparse::GS_DEFAULT, this->mtColoringAlgorithm_);
          mtKernelHandle_->get_point_gs_handle()->set_long_row_threshold(longRowThreshold_);
        } else
          mtKernelHandle_->create_gs_handle(KokkosSparse::CLUSTER_DEFAULT, this->clusterSize_, this->mtColoringAlgorithm_);
      }
      local_matrix_device_type kcsr = crsMat->getLocalMatrixDevice();
      if (PrecType_ == Details::GS2 || PrecType_ == Details::SGS2) {
        // set parameters for two-stage GS
        mtKernelHandle_->set_gs_set_num_inner_sweeps(NumInnerSweeps_);
        mtKernelHandle_->set_gs_set_num_outer_sweeps(NumOuterSweeps_);
        mtKernelHandle_->set_gs_set_inner_damp_factor(InnerDampingFactor_);
        mtKernelHandle_->set_gs_twostage(!InnerSpTrsv_, A_->getLocalNumRows());
        mtKernelHandle_->set_gs_twostage_compact_form(CompactForm_);
      }

      KokkosSparse::gauss_seidel_symbolic(
          mtKernelHandle_.getRawPtr(),
          A_->getLocalNumRows(),
          A_->getLocalNumCols(),
          kcsr.graph.row_map,
          kcsr.graph.entries,
          is_matrix_structurally_symmetric_);
    }
  }  // timing of initialize stops here

  InitializeTime_ += (timer->wallTime() - startTime);
  ++NumInitialize_;
  isInitialized_ = true;
}

namespace Impl {
template <typename BlockDiagView>
struct InvertDiagBlocks {
  typedef typename BlockDiagView::size_type Size;

 private:
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> Unmanaged;
  template <typename View>
  using UnmanagedView = Kokkos::View<typename View::data_type, typename View::array_layout,
                                     typename View::device_type, Unmanaged>;

  typedef typename BlockDiagView::non_const_value_type Scalar;
  typedef typename BlockDiagView::device_type Device;
  typedef Kokkos::View<Scalar**, Kokkos::LayoutRight, Device> RWrk;
  typedef Kokkos::View<int**, Kokkos::LayoutRight, Device> IWrk;

  UnmanagedView<BlockDiagView> block_diag_;
  // TODO Use thread team and scratch memory space. In this first
  // pass, provide workspace for each block.
  RWrk rwrk_buf_;
  UnmanagedView<RWrk> rwrk_;
  IWrk iwrk_buf_;
  UnmanagedView<IWrk> iwrk_;

 public:
  InvertDiagBlocks(BlockDiagView& block_diag)
    : block_diag_(block_diag) {
    const auto blksz = block_diag.extent(1);
    Kokkos::resize(rwrk_buf_, block_diag_.extent(0), blksz);
    rwrk_ = rwrk_buf_;
    Kokkos::resize(iwrk_buf_, block_diag_.extent(0), blksz);
    iwrk_ = iwrk_buf_;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const Size i, int& jinfo) const {
    auto D_cur = Kokkos::subview(block_diag_, i, Kokkos::ALL(), Kokkos::ALL());
    auto ipiv  = Kokkos::subview(iwrk_, i, Kokkos::ALL());
    auto work  = Kokkos::subview(rwrk_, i, Kokkos::ALL());
    int info   = 0;
    Tpetra::GETF2(D_cur, ipiv, info);
    if (info) {
      ++jinfo;
      return;
    }
    Tpetra::GETRI(D_cur, ipiv, work, info);
    if (info) ++jinfo;
  }
};
}  // namespace Impl

template <class MatrixType>
void Relaxation<MatrixType>::computeBlockCrs() {
  using Kokkos::ALL;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::REDUCE_MAX;
  using Teuchos::REDUCE_MIN;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef local_ordinal_type LO;

  const std::string timerName("Ifpack2::Relaxation::computeBlockCrs");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter(timerName);
  if (timer.is_null()) {
    timer = Teuchos::TimeMonitor::getNewCounter(timerName);
  }
  double startTime = timer->wallTime();
  {
    Teuchos::TimeMonitor timeMon(*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
        A_.is_null(), std::runtime_error,
        "Ifpack2::Relaxation::"
        "computeBlockCrs: The input matrix A is null.  Please call setMatrix() "
        "with a nonnull input matrix, then call initialize(), before calling "
        "this method.");
    auto blockCrsA = rcp_dynamic_cast<const block_crs_matrix_type>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(
        blockCrsA.is_null(), std::logic_error,
        "Ifpack2::Relaxation::"
        "computeBlockCrs: A_ is not a BlockCrsMatrix, but it should be if we "
        "got this far.  Please report this bug to the Ifpack2 developers.");

    const scalar_type one = STS::one();

    // Reset state.
    IsComputed_ = false;

    const LO lclNumMeshRows =
        blockCrsA->getCrsGraph().getLocalNumRows();
    const LO blockSize = blockCrsA->getBlockSize();

    if (!savedDiagOffsets_) {
      blockDiag_ = block_diag_type();  // clear it before reallocating
      blockDiag_ = block_diag_type("Ifpack2::Relaxation::blockDiag_",
                                   lclNumMeshRows, blockSize, blockSize);
      if (Teuchos::as<LO>(diagOffsets_.extent(0)) < lclNumMeshRows) {
        // Clear diagOffsets_ first (by assigning an empty View to it)
        // to save memory, before reallocating.
        diagOffsets_ = Kokkos::View<size_t*, device_type>();
        diagOffsets_ = Kokkos::View<size_t*, device_type>("offsets", lclNumMeshRows);
      }
      blockCrsA->getCrsGraph().getLocalDiagOffsets(diagOffsets_);
      TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(diagOffsets_.extent(0)) !=
                                     static_cast<size_t>(blockDiag_.extent(0)),
                                 std::logic_error, "diagOffsets_.extent(0) = " << diagOffsets_.extent(0) << " != blockDiag_.extent(0) = " << blockDiag_.extent(0) << ".  Please report this bug to the Ifpack2 developers.");
      savedDiagOffsets_ = true;
    }
    blockCrsA->getLocalDiagCopy(blockDiag_, diagOffsets_);

    // Use an unmanaged View in this method, so that when we take
    // subviews of it (to get each diagonal block), we don't have to
    // touch the reference count.  Reference count updates are a
    // thread scalability bottleneck and have a performance cost even
    // without using threads.
    unmanaged_block_diag_type blockDiag = blockDiag_;

    if (DoL1Method_ && IsParallel_) {
      const scalar_type two  = one + one;
      const size_t maxLength = A_->getLocalMaxNumRowEntries();
      nonconst_local_inds_host_view_type indices("indices", maxLength);
      nonconst_values_host_view_type values_("values", maxLength * blockSize * blockSize);
      size_t numEntries = 0;

      for (LO i = 0; i < lclNumMeshRows; ++i) {
        // FIXME (mfh 16 Dec 2015) Get views instead of copies.
        blockCrsA->getLocalRowCopy(i, indices, values_, numEntries);
        scalar_type* values = reinterpret_cast<scalar_type*>(values_.data());

        auto diagBlock = Kokkos::subview(blockDiag, i, ALL(), ALL());
        for (LO subRow = 0; subRow < blockSize; ++subRow) {
          magnitude_type diagonal_boost = STM::zero();
          for (size_t k = 0; k < numEntries; ++k) {
            if (indices[k] >= lclNumMeshRows) {
              const size_t offset = blockSize * blockSize * k + subRow * blockSize;
              for (LO subCol = 0; subCol < blockSize; ++subCol) {
                diagonal_boost += STS::magnitude(values[offset + subCol] / two);
              }
            }
          }
          if (STS::magnitude(diagBlock(subRow, subRow)) < L1Eta_ * diagonal_boost) {
            diagBlock(subRow, subRow) += diagonal_boost;
          }
        }
      }
    }

    int info = 0;
    {
      Impl::InvertDiagBlocks<unmanaged_block_diag_type> idb(blockDiag);
      typedef typename unmanaged_block_diag_type::execution_space exec_space;
      typedef Kokkos::RangePolicy<exec_space, LO> range_type;

      Kokkos::parallel_reduce(range_type(0, lclNumMeshRows), idb, info);
      // FIXME (mfh 19 May 2017) Different processes may not
      // necessarily have this error all at the same time, so it would
      // be better just to preserve a local error state and let users
      // check.
      TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::runtime_error, "GETF2 or GETRI failed on = " << info << " diagonal blocks.");
    }

    // In a debug build, do an extra test to make sure that all the
    // factorizations were computed correctly.
#ifdef HAVE_IFPACK2_DEBUG
    const int numResults = 2;
    // Use "max = -min" trick to get min and max in a single all-reduce.
    int lclResults[2], gblResults[2];
    lclResults[0] = info;
    lclResults[1] = -info;
    gblResults[0] = 0;
    gblResults[1] = 0;
    reduceAll<int, int>(*(A_->getGraph()->getComm()), REDUCE_MIN,
                        numResults, lclResults, gblResults);
    TEUCHOS_TEST_FOR_EXCEPTION(gblResults[0] != 0 || gblResults[1] != 0, std::runtime_error,
                               "Ifpack2::Relaxation::compute: When processing the input "
                               "Tpetra::BlockCrsMatrix, one or more diagonal block LU factorizations "
                               "failed on one or more (MPI) processes.");
#endif  // HAVE_IFPACK2_DEBUG
    serialGaussSeidel_ = rcp(new SerialGaussSeidel(blockCrsA, blockDiag_, localSmoothingIndices_, DampingFactor_));
  }  // end TimeMonitor scope

  ComputeTime_ += (timer->wallTime() - startTime);
  ++NumCompute_;
  IsComputed_ = true;
}

template <class MatrixType>
void Relaxation<MatrixType>::compute() {
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
  using Teuchos::reduceAll;
  using LO          = local_ordinal_type;
  using vector_type = Tpetra::Vector<scalar_type, local_ordinal_type,
                                     global_ordinal_type, node_type>;
  using IST         = typename vector_type::impl_scalar_type;
#if KOKKOS_VERSION >= 40799
  using KAT = KokkosKernels::ArithTraits<IST>;
#else
  using KAT = Kokkos::ArithTraits<IST>;
#endif

  const char methodName[] = "Ifpack2::Relaxation::compute";
  const scalar_type zero  = STS::zero();
  const scalar_type one   = STS::one();

  TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null(), std::runtime_error, methodName << ": "
                                                                             "The input matrix A is null.  Please call setMatrix() with a nonnull "
                                                                             "input matrix, then call initialize(), before calling this method.");

  // We don't count initialization in compute() time.
  if (!isInitialized()) {
    initialize();
  }

  if (hasBlockCrsMatrix_) {
    computeBlockCrs();
    return;
  }

  Teuchos::RCP<Teuchos::Time> timer =
      Teuchos::TimeMonitor::getNewCounter(methodName);
  double startTime = timer->wallTime();

  {  // Timing of compute starts here.
    Teuchos::TimeMonitor timeMon(*timer);

    // Precompute some quantities for "fixing" zero or tiny diagonal
    // entries.  We'll only use them if this "fixing" is enabled.
    //
    // SmallTraits covers for the lack of eps() in
    // Teuchos::ScalarTraits when its template parameter is not a
    // floating-point type.  (Ifpack2 sometimes gets instantiated for
    // integer Scalar types.)
    IST oneOverMinDiagVal = KAT::one() / static_cast<IST>(SmallTraits<scalar_type>::eps());
    if (MinDiagonalValue_ != zero)
      oneOverMinDiagVal = KAT::one() / static_cast<IST>(MinDiagonalValue_);

    // It's helpful not to have to recompute this magnitude each time.
    const magnitude_type minDiagValMag = STS::magnitude(MinDiagonalValue_);

    const LO numMyRows = static_cast<LO>(A_->getLocalNumRows());

    TEUCHOS_TEST_FOR_EXCEPTION(NumSweeps_ < 0, std::logic_error, methodName << ": NumSweeps_ = " << NumSweeps_ << " < 0.  "
                                                                                                                  "Please report this bug to the Ifpack2 developers.");
    IsComputed_ = false;

    if (Diagonal_.is_null()) {
      // A_->getLocalDiagCopy fills in all Vector entries, even if the
      // matrix has no stored entries in the corresponding rows.
      Diagonal_ = rcp(new vector_type(A_->getRowMap(), false));
    }

    if (checkDiagEntries_) {
      //
      // Check diagonal elements, replace zero elements with the minimum
      // diagonal value, and store their inverses.  Count the number of
      // "small" and zero diagonal entries while we're at it.
      //
      size_t numSmallDiagEntries        = 0;  // "small" includes zero
      size_t numZeroDiagEntries         = 0;  // # zero diagonal entries
      size_t numNegDiagEntries          = 0;  // # negative (real parts of) diagonal entries
      magnitude_type minMagDiagEntryMag = STM::zero();
      magnitude_type maxMagDiagEntryMag = STM::zero();

      // FIXME: We are extracting the diagonal more than once. That
      //        isn't very efficient, but we shouldn't be checking
      //        diagonal entries at scale anyway.
      A_->getLocalDiagCopy(*Diagonal_);  // slow path for row matrix
      std::unique_ptr<vector_type> origDiag;
      origDiag    = std::unique_ptr<vector_type>(new vector_type(*Diagonal_, Teuchos::Copy));
      auto diag2d = Diagonal_->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto diag   = Kokkos::subview(diag2d, Kokkos::ALL(), 0);
      // As we go, keep track of the diagonal entries with the
      // least and greatest magnitude.  We could use the trick of
      // starting min with +Inf and max with -Inf, but that
      // doesn't work if scalar_type is a built-in integer type.
      // Thus, we have to start by reading the first diagonal
      // entry redundantly.
      if (numMyRows != 0) {
        const magnitude_type d_0_mag = KAT::abs(diag(0));
        minMagDiagEntryMag           = d_0_mag;
        maxMagDiagEntryMag           = d_0_mag;
      }

      // Go through all the diagonal entries.  Compute counts of
      // small-magnitude, zero, and negative-real-part entries.
      // Invert the diagonal entries that aren't too small.  For
      // those too small in magnitude, replace them with
      // 1/MinDiagonalValue_ (or 1/eps if MinDiagonalValue_
      // happens to be zero).
      for (LO i = 0; i < numMyRows; ++i) {
        const IST d_i                = diag(i);
        const magnitude_type d_i_mag = KAT::abs(d_i);
        // Work-around for GitHub Issue #5269.
        // const magnitude_type d_i_real = KAT::real (d_i);
        const auto d_i_real = getRealValue(d_i);

        // We can't compare complex numbers, but we can compare their
        // real parts.
        if (d_i_real < STM::zero()) {
          ++numNegDiagEntries;
        }
        if (d_i_mag < minMagDiagEntryMag) {
          minMagDiagEntryMag = d_i_mag;
        }
        if (d_i_mag > maxMagDiagEntryMag) {
          maxMagDiagEntryMag = d_i_mag;
        }

        if (fixTinyDiagEntries_) {
          // <= not <, in case minDiagValMag is zero.
          if (d_i_mag <= minDiagValMag) {
            ++numSmallDiagEntries;
            if (d_i_mag == STM::zero()) {
              ++numZeroDiagEntries;
            }
            diag(i) = oneOverMinDiagVal;
          } else {
            diag(i) = KAT::one() / d_i;
          }
        } else {  // Don't fix zero or tiny diagonal entries.
          // <= not <, in case minDiagValMag is zero.
          if (d_i_mag <= minDiagValMag) {
            ++numSmallDiagEntries;
            if (d_i_mag == STM::zero()) {
              ++numZeroDiagEntries;
            }
          }
          diag(i) = KAT::one() / d_i;
        }
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
      const Comm<int>& comm = *(A_->getRowMap()->getComm());

      // Compute global min and max magnitude of entries.
      Array<magnitude_type> localVals(2);
      localVals[0] = minMagDiagEntryMag;
      // (- (min (- x))) is the same as (max x).
      localVals[1] = -maxMagDiagEntryMag;
      Array<magnitude_type> globalVals(2);
      reduceAll<int, magnitude_type>(comm, REDUCE_MIN, 2,
                                     localVals.getRawPtr(),
                                     globalVals.getRawPtr());
      globalMinMagDiagEntryMag_ = globalVals[0];
      globalMaxMagDiagEntryMag_ = -globalVals[1];

      Array<size_t> localCounts(3);
      localCounts[0] = numSmallDiagEntries;
      localCounts[1] = numZeroDiagEntries;
      localCounts[2] = numNegDiagEntries;
      Array<size_t> globalCounts(3);
      reduceAll<int, size_t>(comm, REDUCE_SUM, 3,
                             localCounts.getRawPtr(),
                             globalCounts.getRawPtr());
      globalNumSmallDiagEntries_ = globalCounts[0];
      globalNumZeroDiagEntries_  = globalCounts[1];
      globalNumNegDiagEntries_   = globalCounts[2];

      // Compute and save the difference between the computed inverse
      // diagonal, and the original diagonal's inverse.
      vector_type diff(A_->getRowMap());
      diff.reciprocal(*origDiag);
      diff.update(-one, *Diagonal_, one);
      globalDiagNormDiff_ = diff.norm2();
    }

    // Extract the diagonal entries. The CrsMatrix graph version is
    // faster than the row matrix version for subsequent calls to
    // compute(), since it caches the diagonal offsets.

    bool debugAgainstSlowPath = false;

    auto crsMat = Details::getCrsMatrix(A_);

    if (crsMat.get() && crsMat->isFillComplete()) {
      // The invDiagKernel object computes diagonal offsets if
      // necessary. The "compute" call extracts diagonal enties,
      // optionally applies the L1 method and replacement of small
      // entries, and then inverts.
      if (invDiagKernel_.is_null())
        invDiagKernel_ = rcp(new Ifpack2::Details::InverseDiagonalKernel<op_type>(crsMat));
      else
        invDiagKernel_->setMatrix(crsMat);
      invDiagKernel_->compute(*Diagonal_,
                              DoL1Method_ && IsParallel_, L1Eta_,
                              fixTinyDiagEntries_, minDiagValMag);
      savedDiagOffsets_ = true;

      // mfh 27 May 2019: Later on, we should introduce an IFPACK2_DEBUG
      // environment variable to control this behavior at run time.
#ifdef HAVE_IFPACK2_DEBUG
      debugAgainstSlowPath = true;
#endif
    }

    if (crsMat.is_null() || !crsMat->isFillComplete() || debugAgainstSlowPath) {
      // We could not call the CrsMatrix version above, or want to
      // debug by comparing against the slow path.

      // FIXME: The L1 method in this code path runs on host, since we
      //        don't have device access for row matrices.

      RCP<vector_type> Diagonal;
      if (debugAgainstSlowPath)
        Diagonal = rcp(new vector_type(A_->getRowMap()));
      else
        Diagonal = Diagonal_;

      A_->getLocalDiagCopy(*Diagonal);  // slow path for row matrix

      // Setup for L1 Methods.
      // Here we add half the value of the off-processor entries in the row,
      // but only if diagonal isn't sufficiently large.
      //
      // This follows from Equation (6.5) in: Baker, Falgout, Kolev and
      // Yang.  "Multigrid Smoothers for Ultraparallel Computing."  SIAM
      // J. Sci. Comput., Vol. 33, No. 5. (2011), pp. 2864-2887.
      //
      if (DoL1Method_ && IsParallel_) {
        const row_matrix_type& A_row = *A_;
        auto diag                    = Diagonal->getLocalViewHost(Tpetra::Access::ReadWrite);
        const magnitude_type two     = STM::one() + STM::one();
        const size_t maxLength       = A_row.getLocalMaxNumRowEntries();
        nonconst_local_inds_host_view_type indices("indices", maxLength);
        nonconst_values_host_view_type values("values", maxLength);
        size_t numEntries;

        for (LO i = 0; i < numMyRows; ++i) {
          A_row.getLocalRowCopy(i, indices, values, numEntries);
          magnitude_type diagonal_boost = STM::zero();
          for (size_t k = 0; k < numEntries; ++k) {
            if (indices[k] >= numMyRows) {
              diagonal_boost += STS::magnitude(values[k] / two);
            }
          }
          if (KAT::magnitude(diag(i, 0)) < L1Eta_ * diagonal_boost) {
            diag(i, 0) += diagonal_boost;
          }
        }
      }

      //
      // Invert the matrix's diagonal entries (in Diagonal).
      //
      if (fixTinyDiagEntries_) {
        // Go through all the diagonal entries.  Invert those that
        // aren't too small in magnitude.  For those that are too
        // small in magnitude, replace them with oneOverMinDiagVal.
        auto localDiag = Diagonal->getLocalViewDevice(Tpetra::Access::ReadWrite);
        Kokkos::parallel_for(
            Kokkos::RangePolicy<MyExecSpace>(0, localDiag.extent(0)),
            KOKKOS_LAMBDA(local_ordinal_type i) {
              auto d_i                     = localDiag(i, 0);
              const magnitude_type d_i_mag = KAT::magnitude(d_i);
              // <= not <, in case minDiagValMag is zero.
              if (d_i_mag <= minDiagValMag) {
                d_i = oneOverMinDiagVal;
              } else {
                // For Stokhos types, operator/ returns an expression
                // type.  Explicitly convert to IST before returning.
                d_i = IST(KAT::one() / d_i);
              }
              localDiag(i, 0) = d_i;
            });
      } else {  // don't fix tiny or zero diagonal entries
        Diagonal->reciprocal(*Diagonal);
      }

      if (debugAgainstSlowPath) {
        // Validate the fast-path diagonal against the slow-path diagonal.
        Diagonal->update(STS::one(), *Diagonal_, -STS::one());
        const magnitude_type err = Diagonal->normInf();
        // The two diagonals should be exactly the same, so their
        // difference should be exactly zero.
        TEUCHOS_TEST_FOR_EXCEPTION(err > 100 * STM::eps(), std::logic_error, methodName << ": "
                                                                                        << "\"fast-path\" diagonal computation failed.  "
                                                                                           "\\|D1 - D2\\|_inf = "
                                                                                        << err << ".");
      }
    }

    // Count floating-point operations of computing the inverse diagonal.
    //
    // FIXME (mfh 30 Mar 2013) Shouldn't counts be global, not local?
    if (STS::isComplex) {  // magnitude: at least 3 flops per diagonal entry
      ComputeFlops_ += 4.0 * numMyRows;
    } else {
      ComputeFlops_ += numMyRows;
    }

    if (PrecType_ == Ifpack2::Details::MTGS ||
        PrecType_ == Ifpack2::Details::MTSGS ||
        PrecType_ == Ifpack2::Details::GS2 ||
        PrecType_ == Ifpack2::Details::SGS2) {
      // KokkosKernels GaussSeidel Initialization.
      using scalar_view_t = typename local_matrix_device_type::values_type;

      TEUCHOS_TEST_FOR_EXCEPTION(crsMat.is_null(), std::logic_error, methodName << ": "
                                                                                   "Multithreaded Gauss-Seidel methods currently only work "
                                                                                   "when the input matrix is a Tpetra::CrsMatrix.");
      local_matrix_device_type kcsr = crsMat->getLocalMatrixDevice();

      // TODO BMK: This should be ReadOnly, and KokkosKernels should accept a
      // const-valued view for user-provided D^-1. OK for now, Diagonal_ is nonconst.
      auto diagView_2d          = Diagonal_->getLocalViewDevice(Tpetra::Access::ReadWrite);
      scalar_view_t diagView_1d = Kokkos::subview(diagView_2d, Kokkos::ALL(), 0);
      KokkosSparse::gauss_seidel_numeric(
          mtKernelHandle_.getRawPtr(),
          A_->getLocalNumRows(),
          A_->getLocalNumCols(),
          kcsr.graph.row_map,
          kcsr.graph.entries,
          kcsr.values,
          diagView_1d,
          is_matrix_structurally_symmetric_);
    } else if (PrecType_ == Ifpack2::Details::GS || PrecType_ == Ifpack2::Details::SGS) {
      if (crsMat)
        serialGaussSeidel_ = rcp(new SerialGaussSeidel(crsMat, Diagonal_, localSmoothingIndices_, DampingFactor_));
      else
        serialGaussSeidel_ = rcp(new SerialGaussSeidel(A_, Diagonal_, localSmoothingIndices_, DampingFactor_));
    }
  }  // end TimeMonitor scope

  ComputeTime_ += (timer->wallTime() - startTime);
  ++NumCompute_;
  IsComputed_ = true;
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseRichardson(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                           Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y) const {
  using Teuchos::as;
  const double numGlobalRows = as<double>(A_->getGlobalNumRows());
  const double numVectors    = as<double>(X.getNumVectors());
  if (ZeroStartingSolution_) {
    // For the first Richardson sweep, if we are allowed to assume that
    // the initial guess is zero, then Richardson is just alpha times the RHS
    // Compute it as Y(i,j) = DampingFactor_ * X(i,j).
    Y.scale(DampingFactor_, X);

    // Count (global) floating-point operations.  Ifpack2 represents
    // this as a floating-point number rather than an integer, so that
    // overflow (for a very large number of calls, or a very large
    // problem) is approximate instead of catastrophic.
    double flopUpdate = 0.0;
    if (DampingFactor_ == STS::one()) {
      // Y(i,j) = X(i,j): one multiply for each entry of Y.
      flopUpdate = numGlobalRows * numVectors;
    } else {
      // Y(i,j) = DampingFactor_ * X(i,j):
      // One multiplies per entry of Y.
      flopUpdate = numGlobalRows * numVectors;
    }
    ApplyFlops_ += flopUpdate;
    if (NumSweeps_ == 1) {
      return;
    }
  }
  // If we were allowed to assume that the starting guess was zero,
  // then we have already done the first sweep above.
  const int startSweep = ZeroStartingSolution_ ? 1 : 0;

  // Allocate a multivector if the cached one isn't perfect
  updateCachedMultiVector(Y.getMap(), as<size_t>(numVectors));

  for (int j = startSweep; j < NumSweeps_; ++j) {
    // Each iteration: Y = Y + \omega D^{-1} (X - A*Y)
    Tpetra::Details::residual(*A_, Y, X, *cachedMV_);
    Y.update(DampingFactor_, *cachedMV_, STS::one());
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
  const double numGlobalNonzeros = as<double>(A_->getGlobalNumEntries());
  const double dampingFlops      = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  ApplyFlops_ += as<double>(NumSweeps_ - startSweep) * numVectors *
                 (2.0 * numGlobalNonzeros + dampingFlops);
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseJacobi(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y) const {
  using Teuchos::as;
  if (hasBlockCrsMatrix_) {
    ApplyInverseJacobi_BlockCrsMatrix(X, Y);
    return;
  }

  const double numGlobalRows = as<double>(A_->getGlobalNumRows());
  const double numVectors    = as<double>(X.getNumVectors());
  if (ZeroStartingSolution_) {
    // For the first Jacobi sweep, if we are allowed to assume that
    // the initial guess is zero, then Jacobi is just diagonal
    // scaling.  (A_ij * x_j = 0 for i != j, since x_j = 0.)
    // Compute it as Y(i,j) = DampingFactor_ * X(i,j) * D(i).
    Y.elementWiseMultiply(DampingFactor_, *Diagonal_, X, STS::zero());

    // Count (global) floating-point operations.  Ifpack2 represents
    // this as a floating-point number rather than an integer, so that
    // overflow (for a very large number of calls, or a very large
    // problem) is approximate instead of catastrophic.
    double flopUpdate = 0.0;
    if (DampingFactor_ == STS::one()) {
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

  // Allocate a multivector if the cached one isn't perfect
  updateCachedMultiVector(Y.getMap(), as<size_t>(numVectors));

  for (int j = startSweep; j < NumSweeps_; ++j) {
    // Each iteration: Y = Y + \omega D^{-1} (X - A*Y)
    Tpetra::Details::residual(*A_, Y, X, *cachedMV_);
    Y.elementWiseMultiply(DampingFactor_, *Diagonal_, *cachedMV_, STS::one());
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
  const double numGlobalNonzeros = as<double>(A_->getGlobalNumEntries());
  const double dampingFlops      = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  ApplyFlops_ += as<double>(NumSweeps_ - startSweep) * numVectors *
                 (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseJacobi_BlockCrsMatrix(const Tpetra::MultiVector<scalar_type,
                                                                local_ordinal_type,
                                                                global_ordinal_type,
                                                                node_type>& X,
                                      Tpetra::MultiVector<scalar_type,
                                                          local_ordinal_type,
                                                          global_ordinal_type,
                                                          node_type>& Y) const {
  using Tpetra::BlockMultiVector;
  using BMV = BlockMultiVector<scalar_type, local_ordinal_type,
                               global_ordinal_type, node_type>;

  const block_crs_matrix_type* blockMatConst =
      dynamic_cast<const block_crs_matrix_type*>(A_.getRawPtr());
  TEUCHOS_TEST_FOR_EXCEPTION(blockMatConst == nullptr, std::logic_error,
                             "This method should "
                             "never be called if the matrix A_ is not a BlockCrsMatrix.  "
                             "Please report this bug to the Ifpack2 developers.");
  // mfh 23 Jan 2016: Unfortunately, the const cast is necessary.
  // This is because applyBlock() is nonconst (more accurate), while
  // apply() is const (required by Tpetra::Operator interface, but a
  // lie, because it possibly allocates temporary buffers).
  block_crs_matrix_type* blockMat =
      const_cast<block_crs_matrix_type*>(blockMatConst);

  auto meshRowMap                    = blockMat->getRowMap();
  auto meshColMap                    = blockMat->getColMap();
  const local_ordinal_type blockSize = blockMat->getBlockSize();

  BMV X_blk(X, *meshColMap, blockSize);
  BMV Y_blk(Y, *meshRowMap, blockSize);

  if (ZeroStartingSolution_) {
    // For the first sweep, if we are allowed to assume that the
    // initial guess is zero, then block Jacobi is just block diagonal
    // scaling.  (A_ij * x_j = 0 for i != j, since x_j = 0.)
    Y_blk.blockWiseMultiply(DampingFactor_, blockDiag_, X_blk);
    if (NumSweeps_ == 1) {
      return;
    }
  }

  auto pointRowMap     = Y.getMap();
  const size_t numVecs = X.getNumVectors();

  // We don't need to initialize the result MV, since the sparse
  // mat-vec will clobber its contents anyway.
  BMV A_times_Y(*meshRowMap, *pointRowMap, blockSize, numVecs);

  // If we were allowed to assume that the starting guess was zero,
  // then we have already done the first sweep above.
  const int startSweep = ZeroStartingSolution_ ? 1 : 0;

  for (int j = startSweep; j < NumSweeps_; ++j) {
    blockMat->applyBlock(Y_blk, A_times_Y);
    // Y := Y + \omega D^{-1} (X - A*Y).  Use A_times_Y as scratch.
    Y_blk.blockJacobiUpdate(DampingFactor_, blockDiag_,
                            X_blk, A_times_Y, STS::one());
  }
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseSerialGS(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                         Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
                         Tpetra::ESweepDirection direction) const {
  using this_type = Relaxation<MatrixType>;
  // The CrsMatrix version is faster, because it can access the sparse
  // matrix data directly, rather than by copying out each row's data
  // in turn.  Thus, we check whether the RowMatrix is really a
  // CrsMatrix.
  //
  // FIXME (mfh 07 Jul 2013) See note on crs_matrix_type typedef
  // declaration in Ifpack2_Relaxation_decl.hpp header file.  The code
  // will still be correct if the cast fails, but it will use an
  // unoptimized kernel.
  auto blockCrsMat = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(A_);
  auto crsMat      = Details::getCrsMatrix(A_);
  if (blockCrsMat.get()) {
    const_cast<this_type&>(*this).ApplyInverseSerialGS_BlockCrsMatrix(*blockCrsMat, X, Y, direction);
  } else if (crsMat.get()) {
    ApplyInverseSerialGS_CrsMatrix(*crsMat, X, Y, direction);
  } else {
    ApplyInverseSerialGS_RowMatrix(X, Y, direction);
  }
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseSerialGS_RowMatrix(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                                   Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
                                   Tpetra::ESweepDirection direction) const {
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  // Tpetra's GS implementation for CrsMatrix handles zeroing out the
  // starting multivector itself.  The generic RowMatrix version here
  // does not, so we have to zero out Y here.
  if (ZeroStartingSolution_) {
    Y.putScalar(STS::zero());
  }

  size_t NumVectors = X.getNumVectors();

  RCP<MV> Y2;
  if (IsParallel_) {
    if (Importer_.is_null()) {  // domain and column Maps are the same.
      updateCachedMultiVector(Y.getMap(), NumVectors);
    } else {
      updateCachedMultiVector(Importer_->getTargetMap(), NumVectors);
    }
    Y2 = cachedMV_;
  } else {
    Y2 = rcpFromRef(Y);
  }

  for (int j = 0; j < NumSweeps_; ++j) {
    // data exchange is here, once per sweep
    if (IsParallel_) {
      if (Importer_.is_null()) {
        // FIXME (mfh 27 May 2019) This doesn't deep copy -- not
        // clear if this is correct.  Reevaluate at some point.

        *Y2 = Y;  // just copy, since domain and column Maps are the same
      } else {
        Y2->doImport(Y, *Importer_, Tpetra::INSERT);
      }
    }
    serialGaussSeidel_->apply(*Y2, X, direction);

    // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
    if (IsParallel_) {
      Tpetra::deep_copy(Y, *Y2);
    }
  }

  // See flop count discussion in implementation of ApplyInverseGS_CrsMatrix().
  const double dampingFlops      = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors        = as<double>(X.getNumVectors());
  const double numGlobalRows     = as<double>(A_->getGlobalNumRows());
  const double numGlobalNonzeros = as<double>(A_->getGlobalNumEntries());
  ApplyFlops_ += 2.0 * NumSweeps_ * numVectors *
                 (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseSerialGS_CrsMatrix(const crs_matrix_type& A,
                                   const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& B,
                                   Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                                   Tpetra::ESweepDirection direction) const {
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;
  typedef scalar_type Scalar;
  const char prefix[]    = "Ifpack2::Relaxation::SerialGS: ";
  const scalar_type ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  TEUCHOS_TEST_FOR_EXCEPTION(
      !A.isFillComplete(), std::runtime_error,
      prefix << "The matrix is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      NumSweeps_ < 0, std::invalid_argument,
      prefix << "The number of sweeps must be nonnegative, "
                "but you provided numSweeps = "
             << NumSweeps_ << " < 0.");

  if (NumSweeps_ == 0) {
    return;
  }

  RCP<const import_type> importer = A.getGraph()->getImporter();

  RCP<const map_type> domainMap = A.getDomainMap();
  RCP<const map_type> rangeMap  = A.getRangeMap();
  RCP<const map_type> rowMap    = A.getGraph()->getRowMap();
  RCP<const map_type> colMap    = A.getGraph()->getColMap();

#ifdef HAVE_IFPACK2_DEBUG
  {
    // The relation 'isSameAs' is transitive.  It's also a
    // collective, so we don't have to do a "shared" test for
    // exception (i.e., a global reduction on the test value).
    TEUCHOS_TEST_FOR_EXCEPTION(
        !X.getMap()->isSameAs(*domainMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "multivector X be in the domain Map of the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !B.getMap()->isSameAs(*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "B be in the range Map of the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !Diagonal_->getMap()->isSameAs(*rowMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "D be in the row Map of the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !rowMap->isSameAs(*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !domainMap->isSameAs(*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the domain Map and "
        "the range Map of the matrix be the same.");
  }
#endif

  // Fetch a (possibly cached) temporary column Map multivector
  // X_colMap, and a domain Map view X_domainMap of it.  Both have
  // constant stride by construction.  We know that the domain Map
  // must include the column Map, because our Gauss-Seidel kernel
  // requires that the row Map, domain Map, and range Map are all
  // the same, and that each process owns all of its own diagonal
  // entries of the matrix.

  RCP<multivector_type> X_colMap;
  RCP<multivector_type> X_domainMap;
  bool copyBackOutput = false;
  if (importer.is_null()) {
    X_colMap    = Teuchos::rcpFromRef(X);
    X_domainMap = Teuchos::rcpFromRef(X);
    // Column Map and domain Map are the same, so there are no
    // remote entries.  Thus, if we are not setting the initial
    // guess to zero, we don't have to worry about setting remote
    // entries to zero, even though we are not doing an Import in
    // this case.
    if (ZeroStartingSolution_) {
      X_colMap->putScalar(ZERO);
    }
    // No need to copy back to X at end.
  } else {  // Column Map and domain Map are _not_ the same.
    updateCachedMultiVector(colMap, X.getNumVectors());
    X_colMap    = cachedMV_;
    X_domainMap = X_colMap->offsetViewNonConst(domainMap, 0);

    if (ZeroStartingSolution_) {
      // No need for an Import, since we're filling with zeros.
      X_colMap->putScalar(ZERO);
    } else {
      // We could just copy X into X_domainMap.  However, that
      // wastes a copy, because the Import also does a copy (plus
      // communication).  Since the typical use case for
      // Gauss-Seidel is a small number of sweeps (2 is typical), we
      // don't want to waste that copy.  Thus, we do the Import
      // here, and skip the first Import in the first sweep.
      // Importing directly from X effects the copy into X_domainMap
      // (which is a view of X_colMap).
      X_colMap->doImport(X, *importer, Tpetra::INSERT);
    }
    copyBackOutput = true;  // Don't forget to copy back at end.
  }                         // if column and domain Maps are (not) the same

  for (int sweep = 0; sweep < NumSweeps_; ++sweep) {
    if (!importer.is_null() && sweep > 0) {
      // We already did the first Import for the zeroth sweep above,
      // if it was necessary.
      X_colMap->doImport(*X_domainMap, *importer, Tpetra::INSERT);
    }
    // Do local Gauss-Seidel (forward, backward or symmetric)
    serialGaussSeidel_->apply(*X_colMap, B, direction);
  }

  if (copyBackOutput) {
    try {
      deep_copy(X, *X_domainMap);  // Copy result back into X.
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, prefix << "deep_copy(X, *X_domainMap) "
                                              "threw an exception: "
                                           << e.what());
    }
  }

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
  const double dampingFlops      = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors        = X.getNumVectors();
  const double numGlobalRows     = A_->getGlobalNumRows();
  const double numGlobalNonzeros = A_->getGlobalNumEntries();
  ApplyFlops_ += NumSweeps_ * numVectors *
                 (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseSerialGS_BlockCrsMatrix(const block_crs_matrix_type& A,
                                        const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
                                        Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
                                        Tpetra::ESweepDirection direction) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Tpetra::INSERT;

  // FIXME: (tcf) 8/21/2014 -- may be problematic for multiple right hand sides
  //
  //  NOTE (mfh 12 Sep 2014) I don't think it should be a problem for
  //  multiple right-hand sides, unless the input or output MultiVector
  //  does not have constant stride.  We should check for that case
  //  here, in case it doesn't work in localGaussSeidel (which is
  //  entirely possible).
  block_multivector_type yBlock(Y, *(A.getGraph()->getDomainMap()), A.getBlockSize());
  const block_multivector_type xBlock(X, *(A.getColMap()), A.getBlockSize());

  bool performImport = false;
  RCP<block_multivector_type> yBlockCol;
  if (Importer_.is_null()) {
    yBlockCol = rcpFromRef(yBlock);
  } else {
    if (yBlockColumnPointMap_.is_null() ||
        yBlockColumnPointMap_->getNumVectors() != yBlock.getNumVectors() ||
        yBlockColumnPointMap_->getBlockSize() != yBlock.getBlockSize()) {
      yBlockColumnPointMap_ =
          rcp(new block_multivector_type(*(A.getColMap()), A.getBlockSize(),
                                         static_cast<local_ordinal_type>(yBlock.getNumVectors())));
    }
    yBlockCol = yBlockColumnPointMap_;
    if (pointImporter_.is_null()) {
      auto srcMap    = rcp(new map_type(yBlock.getPointMap()));
      auto tgtMap    = rcp(new map_type(yBlockCol->getPointMap()));
      pointImporter_ = rcp(new import_type(srcMap, tgtMap));
    }
    performImport = true;
  }

  multivector_type yBlock_mv;
  multivector_type yBlockCol_mv;
  RCP<const multivector_type> yBlockColPointDomain;
  if (performImport) {  // create views (shallow copies)
    yBlock_mv    = yBlock.getMultiVectorView();
    yBlockCol_mv = yBlockCol->getMultiVectorView();
    yBlockColPointDomain =
        yBlockCol_mv.offsetView(A.getDomainMap(), 0);
  }

  if (ZeroStartingSolution_) {
    yBlockCol->putScalar(STS::zero());
  } else if (performImport) {
    yBlockCol_mv.doImport(yBlock_mv, *pointImporter_, Tpetra::INSERT);
  }

  for (int sweep = 0; sweep < NumSweeps_; ++sweep) {
    if (performImport && sweep > 0) {
      yBlockCol_mv.doImport(yBlock_mv, *pointImporter_, Tpetra::INSERT);
    }
    serialGaussSeidel_->applyBlock(*yBlockCol, xBlock, direction);
    if (performImport) {
      Tpetra::deep_copy(Y, *yBlockColPointDomain);
    }
  }
}

template <class MatrixType>
void Relaxation<MatrixType>::
    ApplyInverseMTGS_CrsMatrix(
        const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& B,
        Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
        Tpetra::ESweepDirection direction) const {
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;

  typedef scalar_type Scalar;

  const char prefix[] = "Ifpack2::Relaxation::(reordered)MTGaussSeidel: ";
  const Scalar ZERO   = Teuchos::ScalarTraits<Scalar>::zero();

  auto crsMat = Details::getCrsMatrix(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(crsMat.is_null(), std::logic_error,
                             "Ifpack2::Relaxation::apply: "
                             "Multithreaded Gauss-Seidel methods currently only work when the "
                             "input matrix is a Tpetra::CrsMatrix.");

  // Teuchos::ArrayView<local_ordinal_type> rowIndices; // unused, as of 04 Jan 2017
  TEUCHOS_TEST_FOR_EXCEPTION(!localSmoothingIndices_.is_null(), std::logic_error,
                             "Ifpack2's implementation of Multithreaded Gauss-Seidel does not "
                             "implement the case where the user supplies an iteration order.  "
                             "This error used to appear as \"MT GaussSeidel ignores the given "
                             "order\".  "
                             "I tried to add more explanation, but I didn't implement \"MT "
                             "GaussSeidel\" [sic].  "
                             "You'll have to ask the person who did.");

  TEUCHOS_TEST_FOR_EXCEPTION(!crsMat->isFillComplete(), std::runtime_error, prefix << "The "
                                                                                      "input CrsMatrix is not fill complete.  Please call fillComplete "
                                                                                      "on the matrix, then do setup again, before calling apply().  "
                                                                                      "\"Do setup\" means that if only the matrix's values have changed "
                                                                                      "since last setup, you need only call compute().  If the matrix's "
                                                                                      "structure may also have changed, you must first call initialize(), "
                                                                                      "then call compute().  If you have not set up this preconditioner "
                                                                                      "for this matrix before, you must first call initialize(), then "
                                                                                      "call compute().");
  TEUCHOS_TEST_FOR_EXCEPTION(NumSweeps_ < 0, std::logic_error, prefix << ": NumSweeps_ = " << NumSweeps_ << " < 0.  Please report this bug to the Ifpack2 "
                                                                                                            "developers.");
  if (NumSweeps_ == 0) {
    return;
  }

  RCP<const import_type> importer = crsMat->getGraph()->getImporter();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !crsMat->getGraph()->getExporter().is_null(), std::runtime_error,
      "This method's implementation currently requires that the matrix's row, "
      "domain, and range Maps be the same.  This cannot be the case, because "
      "the matrix has a nontrivial Export object.");

  RCP<const map_type> domainMap = crsMat->getDomainMap();
  RCP<const map_type> rangeMap  = crsMat->getRangeMap();
  RCP<const map_type> rowMap    = crsMat->getGraph()->getRowMap();
  RCP<const map_type> colMap    = crsMat->getGraph()->getColMap();

#ifdef HAVE_IFPACK2_DEBUG
  {
    // The relation 'isSameAs' is transitive.  It's also a
    // collective, so we don't have to do a "shared" test for
    // exception (i.e., a global reduction on the test value).
    TEUCHOS_TEST_FOR_EXCEPTION(
        !X.getMap()->isSameAs(*domainMap), std::runtime_error,
        "Ifpack2::Relaxation::MTGaussSeidel requires that the input "
        "multivector X be in the domain Map of the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !B.getMap()->isSameAs(*rangeMap), std::runtime_error,
        "Ifpack2::Relaxation::MTGaussSeidel requires that the input "
        "B be in the range Map of the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !rowMap->isSameAs(*rangeMap), std::runtime_error,
        "Ifpack2::Relaxation::MTGaussSeidel requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !domainMap->isSameAs(*rangeMap), std::runtime_error,
        "Ifpack2::Relaxation::MTGaussSeidel requires that the domain Map and "
        "the range Map of the matrix be the same.");
  }
#else
  // Forestall any compiler warnings for unused variables.
  (void)rangeMap;
  (void)rowMap;
#endif  // HAVE_IFPACK2_DEBUG

  // Fetch a (possibly cached) temporary column Map multivector
  // X_colMap, and a domain Map view X_domainMap of it.  Both have
  // constant stride by construction.  We know that the domain Map
  // must include the column Map, because our Gauss-Seidel kernel
  // requires that the row Map, domain Map, and range Map are all
  // the same, and that each process owns all of its own diagonal
  // entries of the matrix.

  RCP<multivector_type> X_colMap;
  RCP<multivector_type> X_domainMap;
  bool copyBackOutput = false;
  if (importer.is_null()) {
    if (X.isConstantStride()) {
      X_colMap    = rcpFromRef(X);
      X_domainMap = rcpFromRef(X);

      // Column Map and domain Map are the same, so there are no
      // remote entries.  Thus, if we are not setting the initial
      // guess to zero, we don't have to worry about setting remote
      // entries to zero, even though we are not doing an Import in
      // this case.
      if (ZeroStartingSolution_) {
        X_colMap->putScalar(ZERO);
      }
      // No need to copy back to X at end.
    } else {
      // We must copy X into a constant stride multivector.
      // Just use the cached column Map multivector for that.
      // force=true means fill with zeros, so no need to fill
      // remote entries (not in domain Map) with zeros.
      updateCachedMultiVector(colMap, as<size_t>(X.getNumVectors()));
      X_colMap = cachedMV_;
      // X_domainMap is always a domain Map view of the column Map
      // multivector.  In this case, the domain and column Maps are
      // the same, so X_domainMap _is_ X_colMap.
      X_domainMap = X_colMap;
      if (!ZeroStartingSolution_) {  // Don't copy if zero initial guess
        try {
          deep_copy(*X_domainMap, X);  // Copy X into constant stride MV
        } catch (std::exception& e) {
          std::ostringstream os;
          os << "Ifpack2::Relaxation::MTGaussSeidel: "
                "deep_copy(*X_domainMap, X) threw an exception: "
             << e.what() << ".";
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what());
        }
      }
      copyBackOutput = true;  // Don't forget to copy back at end.
      /*
      TPETRA_EFFICIENCY_WARNING(
        ! X.isConstantStride (),
        std::runtime_error,
        "MTGaussSeidel: The current implementation of the Gauss-Seidel "
        "kernel requires that X and B both have constant stride.  Since X "
        "does not have constant stride, we had to make a copy.  This is a "
        "limitation of the current implementation and not your fault, but we "
        "still report it as an efficiency warning for your information.");
        */
    }
  } else {  // Column Map and domain Map are _not_ the same.
    updateCachedMultiVector(colMap, as<size_t>(X.getNumVectors()));
    X_colMap = cachedMV_;

    X_domainMap = X_colMap->offsetViewNonConst(domainMap, 0);

    if (ZeroStartingSolution_) {
      // No need for an Import, since we're filling with zeros.
      X_colMap->putScalar(ZERO);
    } else {
      // We could just copy X into X_domainMap.  However, that
      // wastes a copy, because the Import also does a copy (plus
      // communication).  Since the typical use case for
      // Gauss-Seidel is a small number of sweeps (2 is typical), we
      // don't want to waste that copy.  Thus, we do the Import
      // here, and skip the first Import in the first sweep.
      // Importing directly from X effects the copy into X_domainMap
      // (which is a view of X_colMap).
      X_colMap->doImport(X, *importer, Tpetra::CombineMode::INSERT);
    }
    copyBackOutput = true;  // Don't forget to copy back at end.
  }                         // if column and domain Maps are (not) the same

  // The Gauss-Seidel / SOR kernel expects multivectors of constant
  // stride.  X_colMap is by construction, but B might not be.  If
  // it's not, we have to make a copy.
  RCP<const multivector_type> B_in;
  if (B.isConstantStride()) {
    B_in = rcpFromRef(B);
  } else {
    // Range Map and row Map are the same in this case, so we can
    // use the cached row Map multivector to store a constant stride
    // copy of B.
    RCP<multivector_type> B_in_nonconst = rcp(new multivector_type(rowMap, B.getNumVectors()));
    try {
      deep_copy(*B_in_nonconst, B);
    } catch (std::exception& e) {
      std::ostringstream os;
      os << "Ifpack2::Relaxation::MTGaussSeidel: "
            "deep_copy(*B_in_nonconst, B) threw an exception: "
         << e.what() << ".";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what());
    }
    B_in = rcp_const_cast<const multivector_type>(B_in_nonconst);

    /*
    TPETRA_EFFICIENCY_WARNING(
      ! B.isConstantStride (),
      std::runtime_error,
      "MTGaussSeidel: The current implementation requires that B have "
      "constant stride.  Since B does not have constant stride, we had to "
      "copy it into a separate constant-stride multivector.  This is a "
      "limitation of the current implementation and not your fault, but we "
      "still report it as an efficiency warning for your information.");
      */
  }

  local_matrix_device_type kcsr = crsMat->getLocalMatrixDevice();

  bool update_y_vector = true;
  // false as it was done up already, and we dont want to zero it in each sweep.
  bool zero_x_vector = false;

  for (int sweep = 0; sweep < NumSweeps_; ++sweep) {
    if (!importer.is_null() && sweep > 0) {
      // We already did the first Import for the zeroth sweep above,
      // if it was necessary.
      X_colMap->doImport(*X_domainMap, *importer, Tpetra::CombineMode::INSERT);
    }

    if (direction == Tpetra::Symmetric) {
      KokkosSparse::symmetric_gauss_seidel_apply(mtKernelHandle_.getRawPtr(), A_->getLocalNumRows(), A_->getLocalNumCols(),
                                                 kcsr.graph.row_map, kcsr.graph.entries, kcsr.values,
                                                 X_colMap->getLocalViewDevice(Tpetra::Access::ReadWrite),
                                                 B_in->getLocalViewDevice(Tpetra::Access::ReadOnly),
                                                 zero_x_vector, update_y_vector, DampingFactor_, 1);
    } else if (direction == Tpetra::Forward) {
      KokkosSparse::forward_sweep_gauss_seidel_apply(mtKernelHandle_.getRawPtr(), A_->getLocalNumRows(), A_->getLocalNumCols(),
                                                     kcsr.graph.row_map, kcsr.graph.entries, kcsr.values,
                                                     X_colMap->getLocalViewDevice(Tpetra::Access::ReadWrite),
                                                     B_in->getLocalViewDevice(Tpetra::Access::ReadOnly),
                                                     zero_x_vector, update_y_vector, DampingFactor_, 1);
    } else if (direction == Tpetra::Backward) {
      KokkosSparse::backward_sweep_gauss_seidel_apply(mtKernelHandle_.getRawPtr(), A_->getLocalNumRows(), A_->getLocalNumCols(),
                                                      kcsr.graph.row_map, kcsr.graph.entries, kcsr.values,
                                                      X_colMap->getLocalViewDevice(Tpetra::Access::ReadWrite),
                                                      B_in->getLocalViewDevice(Tpetra::Access::ReadOnly),
                                                      zero_x_vector, update_y_vector, DampingFactor_, 1);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument,
          prefix << "The 'direction' enum does not have any of its valid "
                    "values: Forward, Backward, or Symmetric.");
    }
    update_y_vector = false;
  }

  if (copyBackOutput) {
    try {
      deep_copy(X, *X_domainMap);  // Copy result back into X.
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, prefix << "deep_copy(X, *X_domainMap) "
                                              "threw an exception: "
                                           << e.what());
    }
  }

  const double dampingFlops      = (DampingFactor_ == STS::one()) ? 0.0 : 1.0;
  const double numVectors        = as<double>(X.getNumVectors());
  const double numGlobalRows     = as<double>(A_->getGlobalNumRows());
  const double numGlobalNonzeros = as<double>(A_->getGlobalNumEntries());
  double ApplyFlops              = NumSweeps_ * numVectors *
                      (2.0 * numGlobalRows + 2.0 * numGlobalNonzeros + dampingFlops);
  if (direction == Tpetra::Symmetric)
    ApplyFlops *= 2.0;
  ApplyFlops_ += ApplyFlops;
}

template <class MatrixType>
std::string Relaxation<MatrixType>::description() const {
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Relaxation\": {";

  os << "Initialized: " << (isInitialized() ? "true" : "false") << ", "
     << "Computed: " << (isComputed() ? "true" : "false") << ", ";

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
  } else if (PrecType_ == Ifpack2::Details::MTGS) {
    os << "MT Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::MTSGS) {
    os << "MT Symmetric Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::GS2) {
    os << "Two-stage Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::SGS2) {
    os << "Two-stage Symmetric Gauss-Seidel";
  } else {
    os << "INVALID";
  }
  if (hasBlockCrsMatrix_)
    os << ", BlockCrs";

  os << ", "
     << "sweeps: " << NumSweeps_ << ", "
     << "damping factor: " << DampingFactor_ << ", ";

  if (PrecType_ == Ifpack2::Details::MTGS || PrecType_ == Ifpack2::Details::MTSGS) {
    os << "\"relaxation: mtgs cluster size\": " << clusterSize_ << ", ";
    os << "\"relaxation: long row threshold\": " << longRowThreshold_ << ", ";
    os << "\"relaxation: symmetric matrix structure\": " << (is_matrix_structurally_symmetric_ ? "true" : "false") << ", ";
    os << "\"relaxation: relaxation: mtgs coloring algorithm\": ";
    switch (mtColoringAlgorithm_) {
      case KokkosGraph::COLORING_DEFAULT:
        os << "DEFAULT";
        break;
      case KokkosGraph::COLORING_SERIAL:
        os << "SERIAL";
        break;
      case KokkosGraph::COLORING_VB:
        os << "VB";
        break;
      case KokkosGraph::COLORING_VBBIT:
        os << "VBBIT";
        break;
      case KokkosGraph::COLORING_VBCS:
        os << "VBCS";
        break;
      case KokkosGraph::COLORING_VBD:
        os << "VBD";
        break;
      case KokkosGraph::COLORING_VBDBIT:
        os << "VBDBIT";
        break;
      case KokkosGraph::COLORING_EB:
        os << "EB";
        break;
      default:
        os << "*Invalid*";
    }
    os << ", ";
  }

  if (PrecType_ == Ifpack2::Details::GS2 ||
      PrecType_ == Ifpack2::Details::SGS2) {
    os << "outer sweeps: " << NumOuterSweeps_ << ", "
       << "inner sweeps: " << NumInnerSweeps_ << ", "
       << "inner damping factor: " << InnerDampingFactor_ << ", ";
  }

  if (DoL1Method_) {
    os << "use l1: " << DoL1Method_ << ", "
       << "l1 eta: " << L1Eta_ << ", ";
  }

  if (A_.is_null()) {
    os << "Matrix: null";
  } else {
    os << "Global matrix dimensions: ["
       << A_->getGlobalNumRows() << ", " << A_->getGlobalNumCols() << "]"
       << ", Global nnz: " << A_->getGlobalNumEntries();
  }

  os << "}";
  return os.str();
}

template <class MatrixType>
void Relaxation<MatrixType>::
    describe(Teuchos::FancyOStream& out,
             const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_EXTREME;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_NONE;

  const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  const int myRank = this->getComm()->getRank();

  //    none: print nothing
  //     low: print O(1) info from Proc 0
  //  medium:
  //    high:
  // extreme:

  if (vl != VERB_NONE && myRank == 0) {
    // Describable interface asks each implementation to start with a tab.
    OSTab tab1(out);
    // Output is valid YAML; hence the quotes, to protect the colons.
    out << "\"Ifpack2::Relaxation\":" << endl;
    OSTab tab2(out);
    out << "MatrixType: \"" << TypeNameTraits<MatrixType>::name() << "\""
        << endl;
    if (this->getObjectLabel() != "") {
      out << "Label: " << this->getObjectLabel() << endl;
    }
    out << "Initialized: " << (isInitialized() ? "true" : "false") << endl
        << "Computed: " << (isComputed() ? "true" : "false") << endl
        << "Parameters: " << endl;
    {
      OSTab tab3(out);
      out << "\"relaxation: type\": ";
      if (PrecType_ == Ifpack2::Details::JACOBI) {
        out << "Jacobi";
      } else if (PrecType_ == Ifpack2::Details::GS) {
        out << "Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::SGS) {
        out << "Symmetric Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::MTGS) {
        out << "MT Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::MTSGS) {
        out << "MT Symmetric Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::GS2) {
        out << "Two-stage Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::Details::SGS2) {
        out << "Two-stage Symmetric Gauss-Seidel";
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
      if (PrecType_ == Ifpack2::Details::MTGS || PrecType_ == Ifpack2::Details::MTSGS) {
        out << "\"relaxation: mtgs cluster size\": " << clusterSize_ << endl;
        out << "\"relaxation: long row threshold\": " << longRowThreshold_ << endl;
        out << "\"relaxation: symmetric matrix structure\": " << (is_matrix_structurally_symmetric_ ? "true" : "false") << endl;
        out << "\"relaxation: relaxation: mtgs coloring algorithm\": ";
        switch (mtColoringAlgorithm_) {
          case KokkosGraph::COLORING_DEFAULT:
            out << "DEFAULT";
            break;
          case KokkosGraph::COLORING_SERIAL:
            out << "SERIAL";
            break;
          case KokkosGraph::COLORING_VB:
            out << "VB";
            break;
          case KokkosGraph::COLORING_VBBIT:
            out << "VBBIT";
            break;
          case KokkosGraph::COLORING_VBCS:
            out << "VBCS";
            break;
          case KokkosGraph::COLORING_VBD:
            out << "VBD";
            break;
          case KokkosGraph::COLORING_VBDBIT:
            out << "VBDBIT";
            break;
          case KokkosGraph::COLORING_EB:
            out << "EB";
            break;
          default:
            out << "*Invalid*";
        }
        out << endl;
      }
      if (PrecType_ == Ifpack2::Details::GS2 || PrecType_ == Ifpack2::Details::SGS2) {
        out << "\"relaxation: inner damping factor\": " << InnerDampingFactor_ << endl;
        out << "\"relaxation: outer sweeps\" : " << NumOuterSweeps_ << endl;
        out << "\"relaxation: inner sweeps\" : " << NumInnerSweeps_ << endl;
      }
    }
    out << "Computed quantities:" << endl;
    {
      OSTab tab3(out);
      out << "Global number of rows: " << A_->getGlobalNumRows() << endl
          << "Global number of columns: " << A_->getGlobalNumCols() << endl;
    }
    if (checkDiagEntries_ && isComputed()) {
      out << "Properties of input diagonal entries:" << endl;
      {
        OSTab tab3(out);
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
    if (isComputed()) {
      out << "Saved diagonal offsets: "
          << (savedDiagOffsets_ ? "true" : "false") << endl;
    }
    out << "Call counts and total times (in seconds): " << endl;
    {
      OSTab tab3(out);
      out << "initialize: " << endl;
      {
        OSTab tab4(out);
        out << "Call count: " << NumInitialize_ << endl;
        out << "Total time: " << InitializeTime_ << endl;
      }
      out << "compute: " << endl;
      {
        OSTab tab4(out);
        out << "Call count: " << NumCompute_ << endl;
        out << "Total time: " << ComputeTime_ << endl;
      }
      out << "apply: " << endl;
      {
        OSTab tab4(out);
        out << "Call count: " << NumApply_ << endl;
        out << "Total time: " << ApplyTime_ << endl;
      }
    }
  }
}

}  // namespace Ifpack2

#define IFPACK2_RELAXATION_INSTANT(S, LO, GO, N) \
  template class Ifpack2::Relaxation<Tpetra::RowMatrix<S, LO, GO, N>>;

#endif  // IFPACK2_RELAXATION_DEF_HPP
