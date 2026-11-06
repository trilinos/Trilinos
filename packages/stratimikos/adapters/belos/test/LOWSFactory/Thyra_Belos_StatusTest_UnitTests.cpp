// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Operator.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_toString.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"

#include <cmath>


namespace Thyra {


//
// Helper code
//


const std::string matrixFileName = "nos1.mtx";


RCP<const LinearOpBase<double> > getFwdLinearOp()
{
  RCP<const LinearOpBase<double> > fwdLinearOp;
  using crs_type = Tpetra::CrsMatrix<>;
  using op_type = Tpetra::Operator<>;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_type>;
  auto comm = Tpetra::getDefaultComm();
  auto tpetraCrsMatrix = reader_type::readSparseFile(matrixFileName, comm, true);
  fwdLinearOp = createLinearOp(Teuchos::rcp_dynamic_cast<op_type>(tpetraCrsMatrix, true));
  return fwdLinearOp;
}


/** \brief Mock NormInf ReductionFunctional subclass used for unit testing. */
template<class Scalar>
class MockNormInfReductionFunctional : public ReductionFunctional<Scalar> {
protected:

  /** \name Overridded protected functions overridden from ReductionFunctional. */
  //@{

  /** \brief . */
  virtual typename ScalarTraits<Scalar>::magnitudeType
  reduceImpl(const VectorBase<Scalar> &v) const
    { return norm_inf(v); }

  /** \brief . */
  virtual bool isCompatibleImpl( const VectorBase<Scalar> &v ) const
    { return true; }

  //@}

};


/** \brief Non-member constructor.
 *
 * \relates MockNormInfReductionFunctional
 */
template<class Scalar>
RCP<MockNormInfReductionFunctional<Scalar> >
createMockNormReductionFunctional()
{
  return Teuchos::rcp(new MockNormInfReductionFunctional<Scalar>());
}


/** \brief Mock max(NormInf, eps) ReductionFunctional subclass used for unit
 * testing. */
template<class Scalar>
class MockMaxNormInfEpsReductionFunctional : public ReductionFunctional<Scalar> {
protected:

  /** \name Overridded protected functions overridden from ReductionFunctional. */
  //@{

  /** \brief . */
  virtual typename ScalarTraits<Scalar>::magnitudeType
  reduceImpl(const VectorBase<Scalar> &v) const
    {
      typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
      return std::max(norm_inf(v), ScalarTraits<ScalarMag>::eps());
    }

  /** \brief . */
  virtual bool isCompatibleImpl( const VectorBase<Scalar> &v ) const
    { return true; }

  //@}

};


/** \brief Non-member constructor.
 *
 * \relates MockNormInfReductionFunctional
 */
template<class Scalar>
RCP<MockMaxNormInfEpsReductionFunctional<Scalar> >
createMockMaxNormInfEpsReductionFunctional()
{
  return Teuchos::rcp(new MockMaxNormInfEpsReductionFunctional<Scalar>());
}


template<class Scalar>
void runGeneralSolveCriteriaBelosStatusTestCase(
  const SolveCriteria<Scalar> &solveCriteria,
  const Ptr<RCP<const VectorBase<Scalar> > > &x_out,
  const Ptr<RCP<const VectorBase<Scalar> > > &r_out,
  bool &success,
  FancyOStream &out
  )
{

  using Teuchos::describe; using Teuchos::optInArg; using Teuchos::rcpFromRef;
  using Teuchos::toString;

  typedef ScalarTraits<Scalar> ST;

  // A) Set up the linear system
  
  const RCP<const LinearOpBase<Scalar> > fwdOp = getFwdLinearOp();
  out << "\nfwdOp = " << describe(*fwdOp, Teuchos::VERB_MEDIUM) << "\n";
  const RCP<VectorBase<Scalar> > b = createMember(fwdOp->range());
  V_S(b.ptr(), ST::one());
  const RCP<VectorBase<Scalar> > x = createMember(fwdOp->domain());

  // B) Print out the specialized SolveCriteria object

  out << "\nsolveCriteria:\n" << solveCriteria;
  
  // ToDo: Fill in the rest of the fields!

  // C) Solve the system with the given SolveCriteria object

  const int convergenceTestFrequency = 10;

  const RCP<ParameterList> pl = Teuchos::getParametersFromXmlString(
    "<ParameterList name=\"Belos\">"
    "  <Parameter name=\"Solver Type\" type=\"string\" value=\"Pseudo Block GMRES\"/>"
    "  <Parameter name=\"Convergence Test Frequency\" type=\"int\" value=\""+toString(convergenceTestFrequency)+"\"/>"
    "  <ParameterList name=\"Solver Types\">"
    "    <ParameterList name=\"Pseudo Block GMRES\">"
    "      <Parameter name=\"Block Size\" type=\"int\" value=\"1\"/>"
    "      <Parameter name=\"Convergence Tolerance\" type=\"double\" value=\"1e-13\"/>"
    "      <Parameter name=\"Output Frequency\" type=\"int\" value=\""+toString(convergenceTestFrequency)+"\"/>"
    "      <Parameter name=\"Show Maximum Residual Norm Only\" type=\"bool\" value=\"1\"/>"
    "      <Parameter name=\"Maximum Iterations\" type=\"int\" value=\"400\"/>"
    "      <Parameter name=\"Verbosity\" type=\"int\" value=\"1\"/>"
    "    </ParameterList>"
    "  </ParameterList>"
    "</ParameterList>"
    );
  out << "\n\npl:\n" << *pl;

  Thyra::BelosLinearOpWithSolveFactory<Scalar> lowsFactory;
  lowsFactory.setParameterList(pl);
  lowsFactory.setOStream(rcpFromRef(out));
  //lowsFactory.setVerbLevel(Teuchos::VERB_MEDIUM);
  lowsFactory.setVerbLevel(Teuchos::VERB_HIGH);

  // NOTE: To get Belos ouptut to be quite, turn down the Belos "Verbosity"
  // option above.  To get just the StatusTest VerboseObject output, turn up
  // lowsFactory output level.


  const RCP<LinearOpWithSolveBase<Scalar> > lows = linearOpWithSolve<Scalar>(
    lowsFactory, fwdOp);

  V_S(x.ptr(), ST::zero());
  SolveStatus<Scalar> solveStatus = solve<Scalar>(*lows, NOTRANS, *b, x.ptr(),
    optInArg(solveCriteria));
  out << "\nsolveStatus:\n" << solveStatus;

  TEST_COMPARE( solveStatus.achievedTol, <=, solveCriteria.requestedTol );

  // D) Compute the actual residual and return x and r
  
  const RCP<VectorBase<Scalar> > r = b->clone_v();
  fwdOp->apply(NOTRANS, *x, r.ptr(), ST::one(), -ST::one());

  *x_out = x;
  *r_out = r;

}
  


//
// GeneralSolveCriteriaBelosStatusTest Unit Tests
//


TEUCHOS_UNIT_TEST( GeneralSolveCriteriaBelosStatusTest, norm_inf_r_over_norm_inf_r0 )
{
  
  using Teuchos::outArg;

  typedef double Scalar;
  typedef ScalarTraits<Scalar> ST;
  typedef ST::magnitudeType ScalarMag;

  SolveCriteria<Scalar> solveCriteria;
  solveCriteria.solveMeasureType.numerator = SOLVE_MEASURE_NORM_RESIDUAL;
  solveCriteria.numeratorReductionFunc = createMockNormReductionFunctional<Scalar>();
  solveCriteria.solveMeasureType.denominator = SOLVE_MEASURE_NORM_SOLUTION;
  solveCriteria.denominatorReductionFunc = createMockMaxNormInfEpsReductionFunctional<Scalar>();
  solveCriteria.requestedTol = 0.9;

  RCP<const VectorBase<Scalar> > x, r;
  runGeneralSolveCriteriaBelosStatusTestCase(solveCriteria, outArg(x), outArg(r),
    success, out);

  out << "\nChecking convergence ...\n\n";
  
  const ScalarMag r_nrm_inf = norm_inf(*r);
  const ScalarMag x_nrm_inf = norm_inf(*x);
  
  out << "||r||inf = " << r_nrm_inf << "\n";
  out << "||x||inf = " << x_nrm_inf << "\n";
  
  TEST_COMPARE( r_nrm_inf / x_nrm_inf, <=, solveCriteria.requestedTol );

}


TEUCHOS_UNIT_TEST( GeneralSolveCriteriaBelosStatusTest, norm_inf_r_over_1 )
{
  
  using Teuchos::outArg;

  typedef double Scalar;
  typedef ScalarTraits<Scalar> ST;
  typedef ST::magnitudeType ScalarMag;

  SolveCriteria<Scalar> solveCriteria;
  solveCriteria.solveMeasureType.numerator = SOLVE_MEASURE_NORM_RESIDUAL;
  solveCriteria.numeratorReductionFunc = createMockNormReductionFunctional<Scalar>();
  solveCriteria.solveMeasureType.denominator = SOLVE_MEASURE_ONE;
  solveCriteria.requestedTol = 0.9;

  RCP<const VectorBase<Scalar> > x, r;
  runGeneralSolveCriteriaBelosStatusTestCase(solveCriteria, outArg(x), outArg(r),
    success, out);

  out << "\nChecking convergence ...\n\n";
  
  const ScalarMag r_nrm_inf = norm_inf(*r);
  const ScalarMag x_nrm_inf = norm_inf(*x);
  
  out << "||r||inf = " << r_nrm_inf << "\n";
  out << "||x||inf = " << x_nrm_inf << "\n";
  
  TEST_COMPARE( r_nrm_inf, <=, solveCriteria.requestedTol );

}



TEUCHOS_UNIT_TEST( BelosLinearOpWithSolve, pseudoBlockCGConditionEstimate )
{
  /* 1) Test whether a condition number estimate from Belos' PseudoBlockCG is
   *    computed and propagated correctly via the solve status.
   * 2) Reuse the same solver, disable condition number estimation,
   *    and check whether no condition number estimate is returned.
   *    This makes sure that the estimate from the previous run is not exposed.
   * - Do not use a preconditioner.
   * - This test does not verify if the condition number estimate was
   *   computed correctly. It does, however, check whether the estimate
   *   is >= 1 and if it is finite.
   */
  using crs_type = Tpetra::CrsMatrix<double>;
  using map_type = crs_type::map_type;
  using GO = crs_type::global_ordinal_type;
  using op_type = Tpetra::Operator<double>;

  // Create a small symmetric, positive definite, tridiagonal matrix (discrete Laplacian)
  // (  2  -1                  )
  // ( -1   2  -1              )
  // (     -1   2  -1          )
  // (       ..........        )
  // (              -1   2  -1 )
  // (                  -1   2 )
  const GO numRows = 10;
  const auto map = Teuchos::rcp(new map_type(numRows, 0, Tpetra::getDefaultComm()));
  const auto matrix = Teuchos::rcp(new crs_type(map, 3));
  for (const GO row : map->getLocalElementList()) {
    Teuchos::Array<GO> columns;
    Teuchos::Array<double> values;

    // -1 on subdiaonal
    if (row > 0) {
      columns.push_back(row - 1);
      values.push_back(-1.0);
    }

    // 2 on main diagonal
    columns.push_back(row);
    values.push_back(2.0);

    // -1 on superdiagonal
    if (row + 1 < numRows) {
      columns.push_back(row + 1);
      values.push_back(-1.0);
    }

    matrix->insertGlobalValues(row, columns(), values());
  }
  matrix->fillComplete();

  // Create operator from matrix
  const RCP<const LinearOpBase<double> > fwdOp =
    createLinearOp(Teuchos::rcp_dynamic_cast<op_type>(matrix, true));

  // Create right-hand side b of A*x = b with b = ( 1 , ... , 1 )
  const auto b = createMember(fwdOp->range());
  V_S(b.ptr(), 1.0);

  // Initialize solution from A*x = b with x = ( 0 , ... , 0 )
  const auto x = createMember(fwdOp->domain());
  V_S(x.ptr(), 0.0);

  // Set up parameter list for Pseudo Block CG
  const auto paramlist = Teuchos::parameterList();
  paramlist->set("Solver Type", "Pseudo Block CG");
  auto& cgParams = paramlist->sublist("Solver Types").sublist("Pseudo Block CG");
  cgParams.set("Estimate Condition Number", true);
  cgParams.set("Convergence Tolerance", 1e-6);
  cgParams.set("Maximum Iterations", 20);

  // Create solver
  BelosLinearOpWithSolveFactory<double> lowsFactory;
  lowsFactory.setParameterList(paramlist);
  lowsFactory.setOStream(Teuchos::rcpFromRef(out));
  const auto lows = linearOpWithSolve<double>(lowsFactory, fwdOp);

  // Solve A*x = b
  const auto status = solve<double>(*lows, NOTRANS, *b, x.ptr());

  out << "\nCondition number estimation enabled\n";
  out <<   "-----------------------------------\n";
  out << "\nSolve status:\n" << status;
  TEST_EQUALITY(status.solveStatus, SOLVE_STATUS_CONVERGED);

  // The condition number estimate should have been written to parameter list.
  TEST_ASSERT(Teuchos::nonnull(status.extraParameters));
  if (Teuchos::is_null(status.extraParameters)) return;

  // Verify that the two following keys are in the parameter list.
  const char* belosKey = "Belos/Condition Number Estimate";
  const char* genericKey = "Condition Number Estimate";
  TEST_ASSERT(status.extraParameters->isType<double>(belosKey));
  TEST_ASSERT(status.extraParameters->isType<double>(genericKey));
  if (!status.extraParameters->isType<double>(belosKey) ||
      !status.extraParameters->isType<double>(genericKey)) return;

  // Extract the estimate, which is stored in two locations.
  // Make sure the estimates match and check that
  //    1.0 <= condition number estimate < infinity
  const double estimate = status.extraParameters->get<double>(belosKey);
  TEST_EQUALITY(estimate, status.extraParameters->get<double>(genericKey));
  TEST_ASSERT(std::isfinite(estimate));
  TEST_COMPARE(estimate, >=, 1.0);

  // Reinitialize the same LOWS object and reuse Belos' solver manager.
  // Disable condition number estimation.
  // The solver should reset the condition number estimate and not
  // expose the estimate from the previous solve.
  cgParams.set("Estimate Condition Number", false);
  lowsFactory.setParameterList(paramlist);
  initializeOp<double>(lowsFactory, fwdOp, lows.ptr());
  V_S(x.ptr(), 0.0);

  // Solve A*x = b
  const auto status_rerun = solve<double>(*lows, NOTRANS, *b, x.ptr());

  out << "\nCondition number estimation disabled\n";
  out <<   "------------------------------------\n";
  out << "\nSolve status:\n" << status_rerun;
  TEST_EQUALITY(status_rerun.solveStatus, SOLVE_STATUS_CONVERGED);

  // Make sure that condition number estimates have not been returned.
  if (Teuchos::nonnull(status_rerun.extraParameters)) {
    TEST_ASSERT(!status_rerun.extraParameters->isParameter(belosKey));
    TEST_ASSERT(!status_rerun.extraParameters->isParameter(genericKey));
  }
}


} // namespace Thyra
