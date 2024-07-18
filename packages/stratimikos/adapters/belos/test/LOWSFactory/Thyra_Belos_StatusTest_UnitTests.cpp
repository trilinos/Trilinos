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
#include "Thyra_EpetraLinearOp.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_toString.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


//
// Helper code
//


const std::string matrixFileName = "nos1.mtx";


RCP<const LinearOpBase<double> > getFwdLinearOp()
{
  static RCP<const LinearOpBase<double> > fwdLinearOp;
  if (is_null(fwdLinearOp)) {
    Teuchos::RCP<Epetra_CrsMatrix> epetraCrsMatrix;
    EpetraExt::readEpetraLinearSystem( matrixFileName, Epetra_SerialComm(), &epetraCrsMatrix );
    fwdLinearOp = epetraLinearOp(epetraCrsMatrix);
  }
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



} // namespace Thyra
