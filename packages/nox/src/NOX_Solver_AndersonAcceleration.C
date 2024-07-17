// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_AndersonAcceleration.H"    // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Observer.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Assert.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_LineSearch_Generic.H"
#include "NOX_LineSearch_Factory.H"
#include "NOX_SolverStats.hpp"
#include <cmath>

NOX::Solver::AndersonAcceleration::
AndersonAcceleration(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& t,
        const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),                               // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  testPtr(t),
  paramsPtr(p),
  workVec(xGrp->getX().clone(NOX::ShapeCopy)),
  precF(xGrp->getX().clone(NOX::ShapeCopy)),
  oldPrecF(xGrp->getX().clone(NOX::ShapeCopy))
{
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));
  init();
}

// Protected
void NOX::Solver::AndersonAcceleration::init()
{
  // Set up the parameter list
  {
    Teuchos::ParameterList validParams;
    validParams.set("Storage Depth", 2, "max number of previous iterates for which data stored");
    validParams.set("Disable Checks for Unit Testing", false, "If set to true, the check on the storage depth size is disabled so that we can generate some corner cases for unit testing.  WARNING: users should never set this to true!");
    validParams.set("Mixing Parameter", 1.0, "damping factor applied to residuals");
    validParams.set("Reorthogonalization Frequency", 0, "Least-squares problem solved by updating previous QR factorization. Number of iterations between reorthogonalizing columns of Q. Never reorthogonalize if less than 1.");
    validParams.sublist("Preconditioning").set("Precondition", false, "flag for preconditioning");
    validParams.sublist("Preconditioning").set("Recompute Jacobian", false, "set true if preconditioner requires computation of Jacobian");
    validParams.set("Adjust Matrix for Condition Number", false, "If true, the QR matrix will be resized if the condiiton number is greater than the dropTolerance");
    validParams.set("Condition Number Drop Tolerance", 1.0e+12, "If adjusting for condition number, this is the condition number value above which the QR matrix will be resized.");
    validParams.set("Acceleration Start Iteration",1,"The nonlinear iteration where Anderson Acceleration will start. Normally AA starts from the first iteration, but it can be advantageous to delay the start and allow normal picard iteration to get a better initial guess for the AA history.");
    paramsPtr->sublist("Anderson Parameters").validateParametersAndSetDefaults(validParams);
  }

  storeParam = paramsPtr->sublist("Anderson Parameters").get<int>("Storage Depth");
  disableChecksForUnitTesting = paramsPtr->sublist("Anderson Parameters").get<bool>("Disable Checks for Unit Testing");

  if (!disableChecksForUnitTesting) {
    TEUCHOS_TEST_FOR_EXCEPTION((storeParam > solnPtr->getX().length()),std::logic_error,"Error - The \"Storage Depth\" with a value of " << storeParam << " must be less than the number of unknowns in the nonlinear problem which is currently " << solnPtr->getX().length() << ".  This reults in an ill-conditioned F matrix.");
  }

  mixParam = paramsPtr->sublist("Anderson Parameters").get<double>("Mixing Parameter");
  if (storeParam < 0) {
    utilsPtr->out() << "NOX::Solver::AndersonAcceleration::init - "
      << "Storage parameter must be non-negative" << std::endl;
    throw std::runtime_error("NOX Error");
  }
  if ((mixParam < -1.0) || (mixParam == 0) || (mixParam > 1.0)) {
    utilsPtr->out() << "NOX::Solver::AndersonAcceleration::init - "
      << "Mixing parameter must be in [-1,0)U(0,1]" << std::endl;
    throw std::runtime_error("NOX Error");
  }
  orthoFrequency = paramsPtr->sublist("Anderson Parameters").get<int>("Reorthogonalization Frequency");
  precond = paramsPtr->sublist("Anderson Parameters").sublist("Preconditioning").get<bool>("Precondition");
  recomputeJacobian = paramsPtr->sublist("Anderson Parameters").sublist("Preconditioning").get<bool>("Recompute Jacobian");
  adjustForConditionNumber = paramsPtr->sublist("Anderson Parameters").get<bool>("Adjust Matrix for Condition Number");
  dropTolerance = paramsPtr->sublist("Anderson Parameters").get<double>("Condition Number Drop Tolerance");
  accelerationStartIteration = paramsPtr->sublist("Anderson Parameters").get<int>("Acceleration Start Iteration");

  TEUCHOS_TEST_FOR_EXCEPTION((accelerationStartIteration < 1),std::logic_error,"Error - The \"Acceleration Start Iteration\" should be greater than 0!");

  // Initialize
  stepSize = 0.0;
  nIter = 0;
  nStore = 0;
  workVec->init(0.0);
  xMat.resize(0);
  qMat.resize(0);
  for (int ii=0; ii < storeParam; ii++) {
    xMat.push_back(solnPtr->getX().clone(NOX::ShapeCopy));
    qMat.push_back(solnPtr->getX().clone(NOX::ShapeCopy));
  }

  status = NOX::StatusTest::Unconverged;
  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  lineSearchPtr = NOX::LineSearch::
    buildLineSearch(globalDataPtr, paramsPtr->sublist("Line Search"));

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::AndersonAcceleration::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;
  globalDataPtr->getNonConstSolverStatistics()->reset();
  stepSize = 0.0;
  nIter = 0;
  nStore = 0;
  workVec->init(0.0);
  status = NOX::StatusTest::Unconverged;
}

void NOX::Solver::AndersonAcceleration::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  globalDataPtr->getNonConstSolverStatistics()->reset();
  nIter = 0;
  nStore = 0;
  workVec->init(0.0);
  status = NOX::StatusTest::Unconverged;
}

void NOX::Solver::AndersonAcceleration::
reset()
{
  globalDataPtr->getNonConstSolverStatistics()->reset();
  nIter = 0;
  nStore = 0;
  workVec->init(0.0);
  status = NOX::StatusTest::Unconverged;
}

NOX::Solver::AndersonAcceleration::~AndersonAcceleration()
{

}


NOX::StatusTest::StatusType NOX::Solver::AndersonAcceleration::getStatus() const
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::AndersonAcceleration::step()
{
  observer->runPreIterate(*this);
  Teuchos::ParameterList lsParams = paramsPtr->sublist("Direction").sublist("Newton").sublist("Linear Solver");

  // On the first step, do some initializations
  if (nIter == 0) {
    globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearSolves();

    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->out() << "NOX::Solver::AndersonAcceleration::init - "
              << "Unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
    (utilsPtr->isPrintType(NOX::Utils::Warning))) {
      utilsPtr->out() << "Warning: NOX::Solver::AndersonAcceleration::init() - "
              << "The solution passed into the solver (either "
              << "through constructor or reset method) "
              << "is already converged!  The solver wil not "
              << "attempt to solve this system since status is "
              << "flagged as converged." << std::endl;
    }
    printUpdate();

    // First check status
    if (status != NOX::StatusTest::Unconverged) {
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }

    // Apply preconditioner if enabled
    if (precond) {
      if (recomputeJacobian)
        solnPtr->computeJacobian();
      solnPtr->applyRightPreconditioning(false, lsParams, solnPtr->getF(), *oldPrecF);
    }
    else
      *oldPrecF = solnPtr->getF();

    // Copy initial guess to old soln
    *oldSolnPtr = *solnPtr;

    // Compute step then first iterate with a line search.
    workVec->update(mixParam,*oldPrecF);
    bool ok = lineSearchPtr->compute(*solnPtr, stepSize, *workVec, *this);
    if (!ok)
    {
      if (stepSize == 0.0)
      {
        utilsPtr->out() << "NOX::Solver::AndersonAcceleratino::iterate - line search failed" << std::endl;
        status = NOX::StatusTest::Failed;
        observer->runPostIterate(*this);
        printUpdate();
        return status;
      }
      else if (utilsPtr->isPrintType(NOX::Utils::Warning))
        utilsPtr->out() << "NOX::Solver::AndersonAcceleration::iterate - using recovery step for line search" << std::endl;
    }

    // Compute F for the first iterate in case it isn't in the line search
    rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok)
    {
      utilsPtr->out() << "NOX::Solver::AndersonAcceleration::iterate - unable to compute F" << std::endl;
      status = NOX::StatusTest::Failed;
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }

    // Evaluate the current status.
    status = testPtr->checkStatus(*this, checkType);

    //Update iteration count
    nIter++;
    globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Apply preconditioner if enabled
  if (precond) {
    if (recomputeJacobian)
      solnPtr->computeJacobian();
    solnPtr->applyRightPreconditioning(false, lsParams, solnPtr->getF(), *precF);
  }
  else
    *precF = solnPtr->getF();

  // Manage the matrices of past iterates and QR factors
  if (storeParam > 0) {
    if (nIter == accelerationStartIteration) {
      // Initialize
      nStore = 0;
      rMat.shape(0,0);
      oldPrecF->update(1.0, *precF, -1.0);
      qrAdd(*oldPrecF);
      xMat[0]->update(1.0, solnPtr->getX(), -1.0, oldSolnPtr->getX(), 0.0);
    }
    else if (nIter > accelerationStartIteration) {
      if (nStore == storeParam) {
        Teuchos::RCP<NOX::Abstract::Vector> tempPtr = xMat[0];
        for (int ii = 0; ii<nStore-1; ii++)
          xMat[ii] = xMat[ii+1];
        xMat[nStore-1] = tempPtr;
        qrDelete();
      }
      oldPrecF->update(1.0, *precF, -1.0);
      qrAdd(*oldPrecF);
      xMat[nStore-1]->update(1.0, solnPtr->getX(), -1.0, oldSolnPtr->getX(), 0.0);
    }
  }

  // Reorthogonalize 
  if ( (nStore > 1) && (orthoFrequency > 0) )
    if (nIter % orthoFrequency == 0)
      reorthogonalize();

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;
  *oldPrecF = *precF;

  // Adjust for condition number
  if (nStore > 0) {
    Teuchos::LAPACK<int,double> lapack;
    char normType = '1';
    double invCondNum = 0.0;
    int info = 0;
    if ( WORK.size() < static_cast<std::size_t>(4*nStore) )
      WORK.resize(4*nStore);
    if (IWORK.size() < static_cast<std::size_t>(nStore))
      IWORK.resize(nStore);
    lapack.GECON(normType,nStore,rMat.values(),nStore,rMat.normOne(),&invCondNum,&WORK[0],&IWORK[0],&info);
    if (utilsPtr->isPrintType(Utils::Details))
      utilsPtr->out() << "    R condition number estimate ("<< nStore << ") = " << 1.0/invCondNum << std::endl;

    if (adjustForConditionNumber) {
      while ( (1.0/invCondNum > dropTolerance) && (nStore > 1)  ) {
        Teuchos::RCP<NOX::Abstract::Vector> tempPtr = xMat[0];
        for (int ii = 0; ii<nStore-1; ii++)
          xMat[ii] = xMat[ii+1];
        xMat[nStore-1] = tempPtr;
        qrDelete();
        lapack.GECON(normType,nStore,rMat.values(),nStore,rMat.normOne(),&invCondNum,&WORK[0],&IWORK[0],&info);
        if (utilsPtr->isPrintType(Utils::Details))
          utilsPtr->out() << "    Adjusted R condition number estimate ("<< nStore << ") = " << 1.0/invCondNum << std::endl;
      }
    }
  }

  // Solve the least-squares problem.
  Teuchos::SerialDenseMatrix<int,double> gamma(nStore,1), RHS(nStore,1), Rgamma(nStore,1);
  for (int ii = 0; ii<nStore; ii++)
    RHS(ii,0) = precF->innerProduct( *(qMat[ii]) );

  //Back-solve for gamma
  for (int ii = nStore-1; ii>=0; ii--) {
    gamma(ii,0) = RHS(ii,0);
    for (int jj = ii+1; jj<nStore; jj++) {
      gamma(ii,0) -= rMat(ii,jj)*gamma(jj,0);
    }
    gamma(ii,0) /= rMat(ii,ii);
  }

  if (nStore > 0)
    Rgamma.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,mixParam,rMat,gamma,0.0);

  // Compute the step and new solution using the line search.
  workVec->update(mixParam,*precF);
  for (int ii=0; ii<nStore; ii++)
    workVec->update(-gamma(ii,0), *(xMat[ii]), -Rgamma(ii,0),*(qMat[ii]),1.0);
  observer->runPreSolutionUpdate(*workVec,*this);
  bool ok = lineSearchPtr->compute(*solnPtr, stepSize, *workVec, *this);
  observer->runPostSolutionUpdate(*this);
  if (!ok)
  {
    if (stepSize == 0.0)
    {
      utilsPtr->out() << "NOX::Solver::AndersonAcceleration::iterate - line search failed" << std::endl;
      status = NOX::StatusTest::Failed;
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::Solver::AndersonAcceleration::iterate - using recovery step for line search" << std::endl;
  }

  // Compute F for new current solution in case the line search didn't.
  NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utilsPtr->out() << "NOX::Solver::AndersonAcceleration::iterate - unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Update iteration count
  nIter++;
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

  // Evaluate the current status.
  status = testPtr->checkStatus(*this, checkType);

  observer->runPostIterate(*this);

  printUpdate();
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::AndersonAcceleration::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
    step();

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  observer->runPostSolve(*this);

  return status;
}

void NOX::Solver::AndersonAcceleration::qrAdd(NOX::Abstract::Vector& newCol)
{
  // Increment storage depth and resize R factor
  nStore++;
  rMat.reshape(nStore,nStore);

  // Update QR factors
  // Orthogonalize against previous columns once
  for (int ii = 0; ii<nStore-1; ii++) {
    rMat(ii,nStore-1) = qMat[ii]->innerProduct(newCol);
    newCol.update(-rMat(ii,nStore-1),*(qMat[ii]),1.0);
  }
  // Reorthogonalize
  for (int ii = 0; ii<nStore-1; ii++) {
    double alpha = qMat[ii]->innerProduct(newCol);
    newCol.update(-alpha,*(qMat[ii]),1.0);
    rMat(ii,nStore-1) += alpha;
  }
  rMat(nStore-1,nStore-1) = newCol.norm();
  if (!disableChecksForUnitTesting)
    TEUCHOS_TEST_FOR_EXCEPTION((rMat(nStore-1,nStore-1) < 1.0e-16),std::runtime_error,"Error - R factor is singular to machine precision!");
  *(qMat[nStore-1]) = newCol.scale(1.0/rMat(nStore-1,nStore-1));
}

void NOX::Solver::AndersonAcceleration::qrDelete()
{
  // Apply sequence of Givens rotations
  for (int ii = 0; ii<nStore-1; ii++) {
    double temp = sqrt(rMat(ii,ii+1)*rMat(ii,ii+1) + rMat(ii+1,ii+1)*rMat(ii+1,ii+1));
    double c = rMat(ii,ii+1)/temp;
    double s = rMat(ii+1,ii+1)/temp;
    rMat(ii,ii+1) = temp;
    rMat(ii+1,ii+1) = 0;
    for (int jj = ii+2; jj<nStore; jj++) {
      temp = c*rMat(ii,jj) + s*rMat(ii+1,jj);
      rMat(ii+1,jj) = -s*rMat(ii,jj) + c*rMat(ii+1,jj);
      rMat(ii,jj) = temp;
    }
    *workVec = *(qMat[ii]);
    workVec->update(s, *(qMat[ii+1]), c);
    qMat[ii+1]->update(-s, *(qMat[ii]), c);
    *(qMat[ii]) = *workVec;
  }

  // Decrement storage depth and shrink R factor
  nStore--;
  for (int ii=0; ii<nStore; ii++) {
    for (int jj = 0; jj<nStore; jj++)
      rMat(ii,jj) = rMat(ii,jj+1);
  }
  rMat.reshape(nStore,nStore);
}

void NOX::Solver::AndersonAcceleration::reorthogonalize()
{
  // Probably not necessary for fairly small N, but definitely not needed if N=1
  if (nStore > 1) {
    Teuchos::SerialDenseMatrix<int,double> R(nStore,nStore);
    // Reorthogonalize the columns of Q with a modified Gram Schmidt sweep
    for (int ii = 0; ii<nStore; ii++) {
      for (int jj = 0; jj<ii; jj++) {
        R(jj,ii) = qMat[jj]->innerProduct( *(qMat[ii]) );
        qMat[ii]->update(-R(jj,ii),*(qMat[jj]),1.0);
      }
      R(ii,ii) = qMat[ii]->norm();
      if (!disableChecksForUnitTesting)
        TEUCHOS_TEST_FOR_EXCEPTION((R(ii,ii) < 1.0e-16),std::runtime_error,"Error - R factor is singular to machine precision!");
      qMat[ii]->scale(1.0/R(ii,ii));
    }

    // Update the R factor
    Teuchos::SerialDenseMatrix<int,double> rMatCopy(rMat);
    rMat.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,R,rMatCopy,0.0);
  }
}

const NOX::Abstract::Group&
NOX::Solver::AndersonAcceleration::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group&
NOX::Solver::AndersonAcceleration::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

int NOX::Solver::AndersonAcceleration::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList&
NOX::Solver::AndersonAcceleration::getList() const
{
  return *paramsPtr;
}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::AndersonAcceleration::getSolverStatistics() const
{
  return globalDataPtr->getSolverStatistics();
}

// protected
void NOX::Solver::AndersonAcceleration::printUpdate()
{
  double normSoln = 0;

  // Print the status test parameters at each iteration if requested
  if ((status == NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
  }

  // ...But only the print process actually prints the result.
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}

double NOX::Solver::AndersonAcceleration::getStepSize() const
{
  return stepSize;
}
