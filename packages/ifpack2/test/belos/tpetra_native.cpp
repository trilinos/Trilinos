// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_computeRowAndColumnOneNorms.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_leftAndOrRightScaleCrsMatrix.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_ParameterXMLFileReader.hpp"
#include "KokkosBlas1_abs.hpp"

#include <algorithm> // std::transform
#include <cctype> // std::toupper
#include <sstream>
#include <functional>

namespace { // (anonymous)

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
deepCopyFillCompleteCrsMatrix (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
{
  using Teuchos::RCP;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "deepCopyFillCompleteCrsMatrix: Input matrix A must be fillComplete.");
  RCP<crs_matrix_type> A_copy (new crs_matrix_type (A.getCrsGraph ()));
  auto A_copy_lcl = A_copy->getLocalMatrixDevice ();
  auto A_lcl = A.getLocalMatrixDevice ();
  Kokkos::deep_copy (A_copy_lcl.values, A_lcl.values);
  A_copy->fillComplete (A.getDomainMap (), A.getRangeMap ());
  return A_copy;
}

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors =
           ! std::is_same<
               typename Kokkos::ArithTraits<
                 typename ViewType1::non_const_value_type
               >::mag_type,
               typename ViewType2::non_const_value_type
             >::value,
         const int rank = ViewType1::rank>
class ElementWiseMultiply {};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseMultiply<ViewType1,
                          ViewType2,
                          IndexType,
                          takeSquareRootsOfScalingFactors,
                          takeAbsoluteValueOfScalingFactors,
                          1> {
public:
  static_assert (ViewType1::rank == 1, "ViewType1 must be a rank-1 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseMultiply (const ViewType1& X,
                       const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    if (takeAbsoluteValueOfScalingFactors) {
      const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
      const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAM::sqrt (scalFactAbs) : scalFactAbs;
      X_(i) = X_(i) * scalFinalVal;
    }
    else {
      const val_type scalFact = scalingFactors_(i);
      const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAT::sqrt (scalFact) : scalFact;
      X_(i) = X_(i) * scalFinalVal;
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseMultiply<ViewType1,
                          ViewType2,
                          IndexType,
                          takeSquareRootsOfScalingFactors,
                          takeAbsoluteValueOfScalingFactors,
                          2> {
public:
  static_assert (ViewType1::rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseMultiply (const ViewType1& X,
                       const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    for (IndexType j = 0; j < static_cast<IndexType> (X_.extent (1)); ++j) {
      if (takeAbsoluteValueOfScalingFactors) {
        const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
        const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAM::sqrt (scalFactAbs) : scalFactAbs;
        X_(i,j) = X_(i,j) * scalFinalVal;
      }
      else {
        const val_type scalFact = scalingFactors_(i);
        const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAT::sqrt (scalFact) : scalFact;
        X_(i,j) = X_(i,j) * scalFinalVal;
      }
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class MultiVectorViewType,
         class ScalingFactorsViewType,
         class IndexType>
void
elementWiseMultiply (const MultiVectorViewType& X,
                     const ScalingFactorsViewType& scalingFactors,
                     const IndexType numRows,
                     const bool takeSquareRootsOfScalingFactors,
                     const bool takeAbsoluteValueOfScalingFactors =
                       ! std::is_same<
                           typename Kokkos::ArithTraits<
                             typename MultiVectorViewType::non_const_value_type
                           >::mag_type,
                           typename ScalingFactorsViewType::non_const_value_type
                         >::value)
{
  using execution_space = typename MultiVectorViewType::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

template<class MultiVectorType, class ScalingFactorsViewType>
void
elementWiseMultiplyMultiVector (MultiVectorType& X,
                                const ScalingFactorsViewType& scalingFactors,
                                const bool takeSquareRootsOfScalingFactors,
                                const bool takeAbsoluteValueOfScalingFactors =
                                  ! std::is_same<
                                      typename Kokkos::ArithTraits<
                                        typename MultiVectorType::scalar_type
                                      >::mag_type,
                                      typename ScalingFactorsViewType::non_const_value_type
                                    >::value)
{
  using device_type = typename MultiVectorType::device_type;
  using dev_memory_space = typename device_type::memory_space;
  using index_type = typename MultiVectorType::local_ordinal_type;

  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.template getLocalView<dev_memory_space> (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseMultiply (X_lcl_1d, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseMultiply (X_lcl, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
}

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors =
           ! std::is_same<
               typename Kokkos::ArithTraits<
                 typename ViewType1::non_const_value_type
               >::mag_type,
               typename ViewType2::non_const_value_type
             >::value,
         const int rank = ViewType1::rank>
class ElementWiseDivide {};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseDivide<ViewType1,
                        ViewType2,
                        IndexType,
                        takeSquareRootsOfScalingFactors,
                        takeAbsoluteValueOfScalingFactors,
                        1> {
public:
  static_assert (ViewType1::rank == 1, "ViewType1 must be a rank-1 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseDivide (const ViewType1& X,
                     const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    if (takeAbsoluteValueOfScalingFactors) {
      const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
      const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAM::sqrt (scalFactAbs) : scalFactAbs;
      X_(i) = X_(i) / scalFinalVal;
    }
    else {
      const val_type scalFact = scalingFactors_(i);
      const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAT::sqrt (scalFact) : scalFact;
      X_(i) = X_(i) / scalFinalVal;
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseDivide<ViewType1,
                        ViewType2,
                        IndexType,
                        takeSquareRootsOfScalingFactors,
                        takeAbsoluteValueOfScalingFactors,
                        2> {
public:
  static_assert (ViewType1::rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseDivide (const ViewType1& X,
                     const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    for (IndexType j = 0; j < static_cast<IndexType> (X_.extent (1)); ++j) {
      if (takeAbsoluteValueOfScalingFactors) {
        const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
        const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAM::sqrt (scalFactAbs) : scalFactAbs;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
      else {
        const val_type scalFact = scalingFactors_(i);
        const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAT::sqrt (scalFact) : scalFact;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class MultiVectorViewType,
         class ScalingFactorsViewType,
         class IndexType>
void
elementWiseDivide (const MultiVectorViewType& X,
                   const ScalingFactorsViewType& scalingFactors,
                   const IndexType numRows,
                   const bool takeSquareRootsOfScalingFactors,
                   const bool takeAbsoluteValueOfScalingFactors =
                     ! std::is_same<
                         typename Kokkos::ArithTraits<
                           typename MultiVectorViewType::non_const_value_type
                         >::mag_type,
                         typename ScalingFactorsViewType::non_const_value_type
                       >::value)
{
  using execution_space = typename MultiVectorViewType::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

template<class MultiVectorType, class ScalingFactorsViewType>
void
elementWiseDivideMultiVector (MultiVectorType& X,
                              const ScalingFactorsViewType& scalingFactors,
                              const bool takeSquareRootsOfScalingFactors,
                              const bool takeAbsoluteValueOfScalingFactors =
                                ! std::is_same<
                                    typename Kokkos::ArithTraits<
                                      typename MultiVectorType::scalar_type
                                    >::mag_type,
                                    typename ScalingFactorsViewType::non_const_value_type
                                  >::value)
{
  using device_type = typename MultiVectorType::device_type;
  using dev_memory_space = typename device_type::memory_space;
  using index_type = typename MultiVectorType::local_ordinal_type;

  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.template getLocalView<dev_memory_space> (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseDivide (X_lcl_1d, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseDivide (X_lcl, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
}

// See example here:
//
// http://en.cppreference.com/w/cpp/string/byte/toupper
std::string stringToUpper (std::string s)
{
  std::transform (s.begin (), s.end (), s.begin (),
                  [] (unsigned char c) { return std::toupper (c); });
  return s;
}

std::vector<std::string>
splitIntoStrings (const std::string& s,
                  const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<std::string> strings;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      strings.push_back (stringToUpper (token));
    }
    last_pos = cur_pos + size_type (1);
  }
  return strings;
}

template<class T>
std::vector<T>
splitIntoValues (const std::string& s,
                 const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<T> values;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      T val {};
      std::istringstream is (token);
      is >> val;
      if (is) {
        values.push_back (val);
      }
    }
    last_pos = cur_pos + size_type (1);
  }
  return values;
}

// Values of command-line arguments.
struct CmdLineArgs {
  std::string matrixFilename;
  std::string rhsFilename;
  std::string solverTypes = "GMRES";
  std::string orthogonalizationMethod = "ICGS";
  std::string convergenceToleranceValues = "1.0e-2";
  std::string maxIterValues = "40";
  std::string restartLengthValues = "20";
  std::string preconditionerTypes = "RELAXATION";
  bool solverVerbose = false;
  bool equilibrate = false;
  bool assumeSymmetric = false;
  bool assumeZeroInitialGuess = true;
  bool useDiagonalToEquilibrate = false;
};

// Read in values of command-line arguments.
bool
getCmdLineArgs (CmdLineArgs& args, int argc, char* argv[])
{
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("matrixFilename", &args.matrixFilename, "Name of Matrix "
                  "Market file with the sparse matrix A");
  cmdp.setOption ("rhsFilename", &args.rhsFilename, "Name of Matrix Market "
                  "file with the right-hand side vector(s) B");
  cmdp.setOption ("solverTypes", &args.solverTypes,
                  "One or more Belos solver types, "
                  "separated by commas");
  cmdp.setOption ("convergenceTolerances", &args.convergenceToleranceValues,
                  "One or more doubles, separated by commas; each value "
                  "is a convergence tolerance to try");
  cmdp.setOption ("orthogonalizationMethod", &args.orthogonalizationMethod,
                  "Orthogonalization method (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("maxIters", &args.maxIterValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum number of solver iterations to try");
  cmdp.setOption ("restartLengths", &args.restartLengthValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum restart length to try (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("preconditionerTypes", &args.preconditionerTypes,
                  "One or more Ifpack2 preconditioner types, "
                  "separated by commas");
  cmdp.setOption ("solverVerbose", "solverQuiet", &args.solverVerbose,
                  "Whether the Belos solver should print verbose output");
  cmdp.setOption ("equilibrate", "no-equilibrate", &args.equilibrate,
                  "Whether to equilibrate the linear system before solving it");
  cmdp.setOption ("assumeSymmetric", "no-assumeSymmetric",
                  &args.assumeSymmetric, "Whether equilibration should assume "
                  "that the matrix is symmetric");
  cmdp.setOption ("assumeZeroInitialGuess", "assumeNonzeroInitialGuess",
                  &args.assumeZeroInitialGuess, "Whether equilibration should "
                  "assume that the initial guess (vector) is zero");
  cmdp.setOption ("useDiagonalToEquilibrate", "useOneNorms",
                  &args.useDiagonalToEquilibrate,
                  "Whether equilibration should use the matrix's diagonal; "
                  "default is to use row and column one norms");

  auto result = cmdp.parse (argc, argv);
  return result == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL;
}

template<class ScalarType>
struct BelosSolverResult {
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType achievedTolerance;
  int numIters;
  bool converged;
  bool lossOfAccuracy;
};

template<class CrsMatrixType>
class BelosIfpack2Solver {
private:
  using scalar_type = typename CrsMatrixType::scalar_type;
  using local_ordinal_type = typename CrsMatrixType::local_ordinal_type;
  using global_ordinal_type = typename CrsMatrixType::global_ordinal_type;
  using node_type = typename CrsMatrixType::node_type;
  using row_matrix_type =
    Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using OP = Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using problem_type = Belos::LinearProblem<scalar_type, MV, OP>;
  using solver_type = Belos::SolverManager<scalar_type, MV, OP>;
  using preconditioner_type =
    Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;

  void createPreconditioner ()
  {
    rightPrec_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    if (precType_ != "NONE") {
      rightPrec_ =
        Ifpack2::Factory::create<row_matrix_type> (precType_, A_);
      rightPrec_->setParameters (* (precParams_));
    }
  }

  void createSolver ()
  {
    Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
    solver_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    solver_ = belosFactory.create (solverType_, solverParams_);
  }

  void setPreconditionerMatrix (const Teuchos::RCP<const row_matrix_type>& A)
  {
    if (precType_ != "NONE" && rightPrec_.get () != nullptr) {
      // Not all Ifpack2 preconditioners have a setMatrix method.  Use
      // it if it exists; otherwise, start over with a new instance.
      using can_change_matrix = Ifpack2::Details::CanChangeMatrix<row_matrix_type>;
      can_change_matrix* rightPrec = dynamic_cast<can_change_matrix*> (rightPrec_.get ());
      if (rightPrec != nullptr) {
        rightPrec_->setMatrix (A);
      }
      else {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
    }
  }

  void initializePreconditioner ()
  {
    if (precType_ != "NONE") {
      if (rightPrec_.get () == nullptr) {
        createPreconditioner ();
      }
      rightPrec_->initialize ();
    }
  }

  void computePreconditioner ()
  {
    if (precType_ != "NONE") {
      if (rightPrec_.get () == nullptr) {
        createPreconditioner ();
        rightPrec_->initialize ();
      }
      rightPrec_->compute ();
    }
  }

  void equilibrateMatrix ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call compute.");
    if (equilibrate_) {
      using Tpetra::computeRowAndColumnOneNorms;
      using Tpetra::leftAndOrRightScaleCrsMatrix;

      equibResult_ = computeRowAndColumnOneNorms (*A_, assumeSymmetric_);
      if (useDiagonalToEquilibrate_) {
        using device_type = typename node_type::device_type;
        using mag_type = typename Kokkos::ArithTraits<scalar_type>::mag_type;
        using view_type = Kokkos::View<mag_type*, device_type>;

        view_type rowDiagAbsVals ("rowDiagAbsVals",
                                  equibResult_.rowDiagonalEntries.extent (0));
        KokkosBlas::abs (rowDiagAbsVals, equibResult_.rowDiagonalEntries);
        view_type colDiagAbsVals ("colDiagAbsVals",
                                  equibResult_.colDiagonalEntries.extent (0));
        KokkosBlas::abs (colDiagAbsVals, equibResult_.colDiagonalEntries);

        leftAndOrRightScaleCrsMatrix (*A_, rowDiagAbsVals, colDiagAbsVals,
                                      true, true, equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        leftAndOrRightScaleCrsMatrix (*A_, equibResult_.rowNorms,
                                      colScalingFactors, true, true,
                                      equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
    } // if equilibrate_
  }

  void
  preScaleRightHandSides (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type>& B) const
  {
    if (equilibrate_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (B, equibResult_.rowDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        const bool takeSquareRootsOfScalingFactors = equibResult_.assumeSymmetric;
        elementWiseDivideMultiVector (B, equibResult_.rowNorms,
                                      takeSquareRootsOfScalingFactors);
      }
    }
  }

  void
  preScaleInitialGuesses (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type>& X) const
  {
    if (equilibrate_ && ! assumeZeroInitialGuess_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseMultiplyMultiVector (X, equibResult_.colDiagonalEntries,
                                        takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult_.assumeSymmetric;
        elementWiseMultiplyMultiVector (X, colScalingFactors,
                                        takeSquareRootsOfScalingFactors);
      }
    }
  }

  void
  postScaleSolutionVectors (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type>& X) const
  {
    if (equilibrate_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (X, equibResult_.colDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult_.assumeSymmetric;
        elementWiseDivideMultiVector (X, colScalingFactors,
                                      takeSquareRootsOfScalingFactors);
      }
    }
  }

public:
  BelosIfpack2Solver () = default;

  BelosIfpack2Solver (const Teuchos::RCP<CrsMatrixType>& A,
                      const std::string& solverType = "GMRES",
                      const std::string& precType = "NONE") :
    A_ (A),
    solverType_ (solverType),
    precType_ (precType),
    equilibrate_ (false),
    useDiagonalToEquilibrate_ (false)
  {}

  void setMatrix (const Teuchos::RCP<const CrsMatrixType>& A)
  {
    if (A_.get () != A.get ()) {
      setPreconditionerMatrix (A);
      // Belos solvers don't deal well with a complete change of the matrix.
      solver_ = Teuchos::null;
    }
    A_ = A;
  }

  void
  setPreconditionerTypeAndParameters (const std::string& precType,
                                      const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    precType_ = precType;
    precParams_ = params;
    if (rightPrec_.get () != nullptr) {
      if (precType_ != precType) {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        rightPrec_->setParameters (*params);
      }
    }
  }

  void
  setSolverTypeAndParameters (const std::string& solverType,
                              const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    solverType_ = solverType;
    solverParams_ = params;
    if (solver_.get () != nullptr) {
      if (solverType_ != solverType) {
        solver_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        solver_->setParameters (params);
      }
    }
  }

  void
  setEquilibrationParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    if (params.get () != nullptr) {
      equilibrate_ = params->get ("Equilibrate", equilibrate_);
      assumeSymmetric_ = params->get ("Assume symmetric", assumeSymmetric_);
      assumeZeroInitialGuess_ = params->get ("Assume zero initial guess",
                                             assumeZeroInitialGuess_);
      useDiagonalToEquilibrate_ = params->get ("Use diagonal to equilibrate",
                                               useDiagonalToEquilibrate_);
    }
  }

  void initialize ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call initialize.");
    // Calling this implies that the matrix's graph has changed.
    // Belos' solvers don't handle that very well, so best practice is
    // to recreate them in this case.
    solver_ = Teuchos::null;
    initializePreconditioner ();
  }

  void compute ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call compute.");
    equilibrateMatrix ();
    // equilibration changes the matrix, so don't compute the
    // preconditioner until after doing that.
    computePreconditioner ();
  }

  BelosSolverResult<scalar_type>
  solve (Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
         Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& B)
  {
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;

    if (solver_.get () == nullptr) {
      createSolver ();
    }
    if (rightPrec_.get () == nullptr) {
      createPreconditioner ();
      initializePreconditioner ();
      computePreconditioner ();
    }

    preScaleRightHandSides (B);
    preScaleInitialGuesses (X);

    RCP<problem_type> problem (new problem_type (A_, rcpFromRef (X), rcpFromRef (B)));
    if (rightPrec_.get () != nullptr) {
      problem->setRightPrec (rightPrec_);
    }
    problem->setProblem ();
    solver_->setProblem (problem);
    const Belos::ReturnType solveResult = solver_->solve ();

    postScaleSolutionVectors (X);

    typename Teuchos::ScalarTraits<scalar_type>::magnitudeType tol {-1.0};
    try {
      tol = solver_->achievedTol ();
    }
    catch (...) {}
    const int numIters = solver_->getNumIters ();
    const bool converged = (solveResult == Belos::Converged);
    const bool lossOfAccuracy = solver_->isLOADetected ();
    return {tol, numIters, converged, lossOfAccuracy};
  }

  int
  getNumIters () {
    return solver_->getNumIters ();
  }

private:
  Teuchos::RCP<CrsMatrixType> A_;
  Teuchos::RCP<solver_type> solver_;
  Teuchos::RCP<preconditioner_type> rightPrec_;
  using equilibration_result_type = decltype (Tpetra::computeRowAndColumnOneNorms (*A_, false));
  equilibration_result_type equibResult_;

  std::string solverType_;
  Teuchos::RCP<Teuchos::ParameterList> solverParams_;
  std::string precType_;
  Teuchos::RCP<Teuchos::ParameterList> precParams_;

  bool equilibrate_;
  bool assumeSymmetric_;
  bool assumeZeroInitialGuess_;
  bool useDiagonalToEquilibrate_;
};

class TpetraInstance {
public:
  TpetraInstance (int* argc, char*** argv) {
    Tpetra::initialize (argc, argv);
  }

  ~TpetraInstance () {
    Tpetra::finalize ();
  }
};

template<class CrsMatrixType, class MultiVectorType>
bool
solveAndReport (BelosIfpack2Solver<CrsMatrixType>& solver,
                const CrsMatrixType& A_original, // before scaling
                MultiVectorType& X,
                MultiVectorType& B,
                const std::string& solverType,
                const std::string& precType,
                const typename MultiVectorType::mag_type convergenceTolerance,
                const int maxIters,
                const int restartLength,
                const CmdLineArgs& args)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  X.putScalar (0.0);

  RCP<ParameterList> solverParams (new ParameterList ("Belos"));
  if (args.solverVerbose) {
    solverParams->set ("Verbosity",
                       Belos::IterationDetails |
                       Belos::FinalSummary |
                       Belos::StatusTestDetails);
  }
  solverParams->set ("Convergence Tolerance", convergenceTolerance);
  solverParams->set ("Maximum Iterations", maxIters);
  if (solverType == "GMRES") {
    solverParams->set ("Num Blocks", restartLength);
    solverParams->set ("Maximum Restarts", restartLength * maxIters);
    solverParams->set ("Orthogonalization", args.orthogonalizationMethod);
  }

  RCP<ParameterList> precParams (new ParameterList ("Ifpack2"));
  if (precType == "RELAXATION") {
    precParams->set ("relaxation: type", "Symmetric Gauss-Seidel");
  }

  RCP<ParameterList> equibParams (new ParameterList ("Equilibration"));
  equibParams->set ("Equilibrate", args.equilibrate);
  equibParams->set ("Assume symmetric", args.assumeSymmetric);
  equibParams->set ("Assume zero initial guess",
                    args.assumeZeroInitialGuess);
  equibParams->set ("Use diagonal to equilibrate",
                    args.useDiagonalToEquilibrate);

  solver.setSolverTypeAndParameters (solverType, solverParams);
  solver.setPreconditionerTypeAndParameters (precType, precParams);
  solver.setEquilibrationParameters (equibParams);

  solver.initialize ();
  solver.compute ();

  // Keep this around for later computation of the explicit residual
  // norm.  If the solver equilibrates, it will modify the original B.
  MultiVectorType R (B, Teuchos::Copy);

  // Compute ||B||_2.
  using mag_type = typename MultiVectorType::mag_type;
  Teuchos::Array<mag_type> norms (R.getNumVectors ());
  R.norm2 (norms ());
  mag_type B_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < B.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    B_norm2_max = norms[j] < B_norm2_max ? B_norm2_max : norms[j];
  }

  // Solve the linear system AX=B.
  auto result = solver.solve (X, B);

  // Compute the actual residual norm ||B - A*X||_2.
  using scalar_type = typename MultiVectorType::scalar_type;
  const scalar_type ONE = Teuchos::ScalarTraits<scalar_type>::one ();
  A_original.apply (X, R, Teuchos::NO_TRANS, -ONE, ONE); // R := -A*X + B
  R.norm2 (norms ());

  mag_type R_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < R.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    R_norm2_max = norms[j] < R_norm2_max ? R_norm2_max : norms[j];
  }

  X.norm2 (norms ());
  mag_type X_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < R.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    X_norm2_max = norms[j] < X_norm2_max ? X_norm2_max : norms[j];
  }

  const int myRank = X.getMap ()->getComm ()->getRank ();
  if (myRank == 0) {
    using std::cout;
    using std::endl;
    cout << "Solver:" << endl
         << "  Solver type: " << solverType << endl
         << "  Preconditioner type: " << precType << endl
         << "  Convergence tolerance: " << convergenceTolerance << endl
         << "  Maximum number of iterations: " << maxIters << endl;
    if (solverType == "GMRES") {
      cout << "  Restart length: " << restartLength << endl
           << "  Orthogonalization method: " << args.orthogonalizationMethod << endl;
    }
    cout << "Results:" << endl
         << "  Converged: " << (result.converged ? "true" : "false") << endl
         << "  Number of iterations: " << result.numIters << endl
         << "  Achieved tolerance: " << result.achievedTolerance << endl
         << "  Loss of accuracy: " << result.lossOfAccuracy << endl
         << "  ||B-A*X||_2: " << R_norm2_max << endl
         << "  ||B||_2: " << B_norm2_max << endl
         << "  ||X||_2: " << X_norm2_max << endl;
    if (B_norm2_max != Kokkos::ArithTraits<mag_type>::zero ()) {
      cout << "  ||B-A*X||_2 / ||B||_2: " << (R_norm2_max / B_norm2_max)
           << endl;
    }
    cout << endl;
  }
  return result.converged;
}

} // namespace (anonymous)

namespace BelosTpetra {
namespace Impl {
  // extern void register_Cg (const bool verbose);
  // extern void register_CgPipeline (const bool verbose);
  // extern void register_CgSingleReduce (const bool verbose);
  extern void register_Gmres (const bool verbose);
  extern void register_GmresPipeline (const bool verbose);
  extern void register_GmresS (const bool verbose);
  extern void register_GmresSingleReduce (const bool verbose);
  extern void register_GmresSstep (const bool verbose);
} // namespace Impl
} // namespace BelosTpetra

int
main (int argc, char* argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::cerr;
  using std::endl;
  using crs_matrix_type = Tpetra::CrsMatrix<>;
  using MV = Tpetra::MultiVector<>;
  // using mag_type = MV::mag_type;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;

  constexpr bool verbose = false;
  // BelosTpetra::Impl::register_Cg (verbose);
  // BelosTpetra::Impl::register_CgPipeline (verbose);
  // BelosTpetra::Impl::register_CgSingleReduce (verbose);
  BelosTpetra::Impl::register_Gmres (verbose);
  BelosTpetra::Impl::register_GmresPipeline (verbose);
  BelosTpetra::Impl::register_GmresS (verbose);
  BelosTpetra::Impl::register_GmresSingleReduce (verbose);
  BelosTpetra::Impl::register_GmresSstep (verbose);

  TpetraInstance tpetraInstance (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();

  // Get command-line arguments.
  CmdLineArgs args;
  const bool gotCmdLineArgs = getCmdLineArgs (args, argc, argv);
  if (! gotCmdLineArgs) {
    if (comm->getRank () == 0) {
      cerr << "Failed to get command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }

  if (args.matrixFilename == "") {
    if (comm->getRank () == 0) {
      cerr << "Must specify sparse matrix filename!" << endl;
    }
    return EXIT_FAILURE;
  }

  std::vector<std::string> solverTypes;
  if (args.solverTypes == "") {
    solverTypes = {"GMRES", "TFQMR", "BICGSTAB"};
  }
  else {
    solverTypes = splitIntoStrings (args.solverTypes);
  }

  std::vector<std::string> preconditionerTypes;
  if (args.preconditionerTypes == "") {
    preconditionerTypes = {"RELAXATION"};
  }
  else {
    preconditionerTypes = splitIntoStrings (args.preconditionerTypes);
  }

  std::vector<int> maxIterValues;
  if (args.maxIterValues == "") {
    maxIterValues = {100};
  }
  else {
    maxIterValues = splitIntoValues<int> (args.maxIterValues);
  }

  std::vector<int> restartLengthValues;
  if (args.restartLengthValues == "") {
    restartLengthValues = {20};
  }
  else {
    restartLengthValues = splitIntoValues<int> (args.restartLengthValues);
  }

  std::vector<double> convergenceToleranceValues;
  if (args.convergenceToleranceValues == "") {
    convergenceToleranceValues = {20};
  }
  else {
    convergenceToleranceValues =
      splitIntoValues<double> (args.convergenceToleranceValues);
  }

  // Read sparse matrix A from Matrix Market file.
  RCP<crs_matrix_type> A =
    reader_type::readSparseFile (args.matrixFilename, comm);
  if (A.get () == nullptr) {
    if (comm->getRank () == 0) {
      cerr << "Failed to load sparse matrix A from file "
        "\"" << args.matrixFilename << "\"!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Read right-hand side vector(s) B from Matrix Market file, or
  // generate B if file not specified.
  RCP<MV> B;
  RCP<MV> X;
  if (args.rhsFilename == "") {
    B = Teuchos::rcp (new MV (A->getRangeMap (), 1));
    X = Teuchos::rcp (new MV (A->getDomainMap (), 1));
    X->putScalar (1.0);
    A->apply (*X, *B);
    X->putScalar (0.0);

    double norm {0.0};
    using host_device_type = Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>;
    Kokkos::View<double*, host_device_type> normView (&norm, B->getNumVectors ());
    B->norm2 (normView);
    if (norm != 0.0) {
      B->scale (1.0 / norm);
    }
  }
  else {
    auto map = A->getRangeMap ();
    B = reader_type::readDenseFile (args.rhsFilename, comm, map);
    if (B.get () == nullptr) {
      if (comm->getRank () == 0) {
        cerr << "Failed to load right-hand side vector(s) from file \""
             << args.rhsFilename << "\"!" << endl;
      }
      return EXIT_FAILURE;
    }
    X = Teuchos::rcp (new MV (A->getDomainMap (), B->getNumVectors ()));
  }

  auto A_original = deepCopyFillCompleteCrsMatrix (*A);

  // Create the solver.
  BelosIfpack2Solver<crs_matrix_type> solver (A);

  // Solve the linear system using various solvers and preconditioners.
  bool success = true;
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  for (std::string solverType : solverTypes) {
    for (std::string precType : preconditionerTypes) {
      for (int maxIters : maxIterValues) {
        for (int restartLength : restartLengthValues) {
          for (double convTol : convergenceToleranceValues) {
            bool converged = 
              solveAndReport (solver, *A_original, *X, *B,
                              solverType,
                              precType,
                              convTol,
                              maxIters,
                              restartLength,
                              args);
            if (!converged)
              success = false;
          }
        }
      }
    }
  }

  if (success) {
    *out << "End Result: TEST PASSED\n";
  }
  else {
    *out << "End Result: TEST FAILED\n";
  }

  return EXIT_SUCCESS;
}
