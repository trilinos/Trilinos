//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_ErrorCheck.H"	// class definition
#include "LOCA_PhaseTransition_ExtendedGroup.H"	// class definition

extern "C" {
extern void post_process(double**, char*, int*, double*, int, int);
}

LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup(
      const Teuchos::RCP<LOCA::GlobalData> gD,
      const Teuchos::RCP<Teuchos::ParameterList>& ptParams_,
      const Teuchos::RCP<LOCA::PhaseTransition::AbstractGroup>& grp_) :
  LOCA::Extended::MultiAbstractGroup(),
  LOCA::MultiContinuation::AbstractGroup(),
  grp(grp_), // Underlying group for regular system of size n
  globalData(gD),
  normF(0)
{
  const char *func = "LOCA::PhaseTransition::ExtendedGroup()";

  // Get second solution part od Extended Solution
  if (!ptParams_->isParameter("Second Solution Vector")) {
    globalData->locaErrorCheck->throwError(func,
                                 "\"Second Solution Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> secondSolution =
    (*ptParams_).INVALID_TEMPLATE_QUALIFIER
    get<Teuchos::RCP<NOX::Abstract::Vector> >("Second Solution Vector");
  
  // Get Parameter componenet of extended vector, ptp
  if (!ptParams_->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
                                 "\"Bifurcation Parameter\" name is not set!");
  }

  std::string bifParamName = ptParams_->get("Bifurcation Parameter",
                                                 "None");
  const ParameterVector& p = grp->getParams();
  bifParamID = p.getIndex(bifParamName);
  double ptp = grp->getParam(bifParamID);


  // Finally, contruct extended vector, and spave for other vectors
  xVector = Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(
                         globalData, grp->getX(), *secondSolution, ptp)) ;
  fVector =  Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*xVector, NOX::ShapeCopy));
  newtonVector =  Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*xVector, NOX::ShapeCopy));
}

LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup(
    const LOCA::PhaseTransition::ExtendedGroup& source, NOX::CopyType type) :
  LOCA::Extended::MultiAbstractGroup(),
  LOCA::MultiContinuation::AbstractGroup(),
  grp(source.grp),
  xVector(Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*source.xVector, type))), 
  fVector(Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*source.fVector, type))), 
  newtonVector(Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*source.newtonVector, type))), 
  globalData(source.globalData)
{
 
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normF = source.normF;
    break;

  case NOX::ShapeCopy:
    normF = 0.0;
    break;

  default:
    std::cerr << "LOCA::PhaseTransition::ExtendedGroup - invalid CopyType for copy constructor." << std::endl;
    throw "LOCA::PhaseTransition:Error";
  }

}

LOCA::PhaseTransition::ExtendedGroup::~ExtendedGroup() 
{
}

void LOCA::PhaseTransition::ExtendedGroup::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

Teuchos::RCP<NOX::Abstract::Group> LOCA::PhaseTransition::ExtendedGroup::
clone(NOX::CopyType type) const 
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp = 
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedGroup(*this, type));
  return newgrp;
}

NOX::Abstract::Group& LOCA::PhaseTransition::ExtendedGroup::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedGroup&> (source));
}

NOX::Abstract::Group& LOCA::PhaseTransition::ExtendedGroup::operator=(const ExtendedGroup& source)
{
  if (this != &source) {

    grp = source.grp;

    // Copy the xVector
    xVector = source.xVector;

    // Update the isValidVectors
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    
    // Only copy vectors that are valid
    if (isValidF) {
      fVector = source.fVector;
      normF = source.normF;
    }

    if (isValidNewton)
      newtonVector = source.newtonVector;
    
  }

  return *this;
}

void LOCA::PhaseTransition::ExtendedGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&> (y));
}

void LOCA::PhaseTransition::ExtendedGroup::setX(const LOCA::PhaseTransition::ExtendedVector& y) 
{
  resetIsValid();
  *xVector = y;
}

void LOCA::PhaseTransition::ExtendedGroup::computeX(const NOX::Abstract::Group& exgrp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const ExtendedGroup& trgrp = dynamic_cast<const ExtendedGroup&> (exgrp);
  const LOCA::PhaseTransition::ExtendedVector& trd = dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&> (d);
  computeX(trgrp, trd, step); 
}

void LOCA::PhaseTransition::ExtendedGroup::computeX(const ExtendedGroup& exgrp,
        const LOCA::PhaseTransition::ExtendedVector& d, double step) 
{
  resetIsValid();
  xVector->update(1.0, *exgrp.xVector, step, d);
}

NOX::Abstract::Group::ReturnType LOCA::PhaseTransition::ExtendedGroup::computeF() 
{
  if (isValidF) 
    return NOX::Abstract::Group::Ok;

  grp->setParam(bifParamID, xVector->PTP());

  grp->setX(*xVector->X1());
  grp->computeF();
  *fVector->X1() = grp->getF();
  double omega1 = grp->computeFreeEnergy();

  grp->setX(*xVector->X2());
  grp->computeF();
  *fVector->X2() = grp->getF();
  double omega2 = grp->computeFreeEnergy();

  fVector->PTP() = omega1 - omega2;

  normF = fVector->norm();

  isValidF = true;
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType LOCA::PhaseTransition::ExtendedGroup::computeJacobian() 
{
  // To save memory, only compute Jacobians at withing applyJacobianInverse. 
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType LOCA::PhaseTransition::ExtendedGroup::computeNewton(Teuchos::ParameterList& p) 
{
  if (isNewton())
    return NOX::Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid F" << std::endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << std::endl;
    throw "NOX Error";
  }

  NOX::Abstract::Group::ReturnType status = applyJacobianInverse(p, *fVector, *newtonVector);
  isValidNewton = (status == NOX::Abstract::Group::Ok);

  // Scale soln by -1
  newtonVector->scale(-1.0);

  // Return solution
  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::PhaseTransition::ExtendedGroup::applyJacobian(const NOX::Abstract::Vector& input, 
				  NOX::Abstract::Vector& result) const
{
  const LOCA::PhaseTransition::ExtendedVector& lapackinput = 
     dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&> (input);
  LOCA::PhaseTransition::ExtendedVector& lapackresult =
     dynamic_cast<LOCA::PhaseTransition::ExtendedVector&> (result);
  return applyJacobian(lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
LOCA::PhaseTransition::ExtendedGroup::applyJacobian(const LOCA::PhaseTransition::ExtendedVector& input,
                                          LOCA::PhaseTransition::ExtendedVector& result) const
{
  std::cout << "ERROR:  Apply Jacobian not implemented for ExtendedGroup !!!!" << std::endl;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  const LOCA::PhaseTransition::ExtendedVector& lapackinput =
    dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&> (input);
  LOCA::PhaseTransition::ExtendedVector& lapackresult =
    dynamic_cast<LOCA::PhaseTransition::ExtendedVector&> (result); 
  return applyJacobianInverse(p, lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const LOCA::PhaseTransition::ExtendedVector& input, 
					 LOCA::PhaseTransition::ExtendedVector& result) const 
{
  // Algorithm from Equations 18-26 of Salinger&Frink JChemPhys (2003).
  // This implementation does not assume that input=f so it can solve
  // other RHS's for arclength or bifurcations of PhaseTransitions
  
  const double eps=1.0e-7;
  double perturb = eps * (eps + fabs(xVector->PTP()));
  double pertParam = xVector->PTP() + perturb;
  perturb = pertParam - xVector->PTP();  //improves roundoff accuracy

  // temporary space...
  LOCA::PhaseTransition::ExtendedVector bdvec(result);
  Teuchos::RCP<NOX::Abstract::Vector> fVec = result.X1()->clone();;
  Teuchos::RCP<NOX::Abstract::Vector> fPertVec = result.X1()->clone();;
 
  // First matrix block...

  grp->setX(*xVector->X1());

  //perturb parameter
  grp->setParam(bifParamID, pertParam);
  grp->computeF();
  *fPertVec = grp->getF();
  double omega1pert = grp->computeFreeEnergy();

  //unperturb parameter, compute df/dp and load Jacobian
  grp->setParam(bifParamID, xVector->PTP());
  grp->computeJacobian(); // computes resid as well
  grp->computeF();
  *fVec = grp->getF();
  fPertVec->update(-1.0, *fVec, 1.0);
  fPertVec->scale(1.0/perturb);
  double omega1 = grp->computeFreeEnergy();

  // Do two solves with same matrix for "a" and "b" vectors, Equations 18,19
  grp->applyJacobianInverse(p, *input.X1(), *result.X1());
  grp->applyJacobianInverse(p, *fPertVec, *bdvec.X1());

  // Second matrix block... same but with X2

  grp->setX(*xVector->X2());

  //perturb parameter
  grp->setParam(bifParamID, pertParam);
  grp->computeF();
  *fPertVec = grp->getF();
  double omega2pert = grp->computeFreeEnergy();

  //unperturb parameter, compute df/dp and load Jacobian
  grp->setParam(bifParamID, xVector->PTP());
  grp->computeJacobian(); // computes resid as well
  grp->computeF();
  *fVec = grp->getF();
  fPertVec->update(-1.0, *fVec, 1.0);
  fPertVec->scale(1.0/perturb);
  double omega2 = grp->computeFreeEnergy();

  // Do two solves with same matrix for "c" and "d" vectors, Equations 20,21
  grp->applyJacobianInverse(p, *input.X2(), *result.X2());
  grp->applyJacobianInverse(p, *fPertVec, *bdvec.X2());

  // Compute contributions of equal-energy constraint equation
  
  double g = input.PTP(); // may not be omega1 - omega2 !
  double dgdp = (omega1pert - omega1 - omega2pert + omega2) / perturb;
  
  //For each of 4 terms: calc perturbation, calc perturbed energy, diff
  //Equation 22
  perturb = eps * xVector->X1()->norm() / (eps + result.X1()->norm());
  fPertVec->update(1.0, *xVector->X1(), perturb, *result.X1(), 0.0);
  grp->setX(*fPertVec);
  double dOmdx1a = (grp->computeFreeEnergy() - omega1) / perturb;
  
  perturb = eps * xVector->X2()->norm() / (eps + result.X2()->norm());
  fPertVec->update(1.0, *xVector->X2(), perturb, *result.X2(), 0.0);
  grp->setX(*fPertVec);
  double dOmdx2c = (grp->computeFreeEnergy() - omega2) / perturb;
  
  perturb = eps * xVector->X1()->norm() / (eps + bdvec.X1()->norm());
  fPertVec->update(1.0, *xVector->X1(), perturb, *bdvec.X1(), 0.0);
  grp->setX(*fPertVec);
  double dOmdx1b = (grp->computeFreeEnergy() - omega1) / perturb;
  
  perturb = eps * xVector->X2()->norm() / (eps + bdvec.X2()->norm());
  fPertVec->update(1.0, *xVector->X2(), perturb, *bdvec.X2(), 0.0);
  grp->setX(*fPertVec);
  double dOmdx2d = (grp->computeFreeEnergy() - omega2) / perturb;

  // Equations 23 and 24
  double delta_param = ( g - dOmdx1a + dOmdx2c) / (dgdp + dOmdx1b - dOmdx2d);
  bdvec.PTP() = 0.0;
  result.update(delta_param, bdvec, 1.0);
  result.PTP() = -delta_param; // Negative sign added, and somehow works

  return NOX::Abstract::Group::Ok;
}

bool LOCA::PhaseTransition::ExtendedGroup::isF() const 
{   
  return isValidF;
}

bool LOCA::PhaseTransition::ExtendedGroup::isJacobian() const 
{  
  return isValidJacobian;
}

bool LOCA::PhaseTransition::ExtendedGroup::isNewton() const 
{   
  return isValidNewton;
}

const NOX::Abstract::Vector& LOCA::PhaseTransition::ExtendedGroup::getX() const 
{
  return *xVector;
}

const NOX::Abstract::Vector& LOCA::PhaseTransition::ExtendedGroup::getF() const 
{  
  return *fVector;
}

double LOCA::PhaseTransition::ExtendedGroup::getNormF() const
{
  return normF;
}

const NOX::Abstract::Vector& LOCA::PhaseTransition::ExtendedGroup::getNewton() const 
{
  return *newtonVector;
}

const NOX::Abstract::Vector& LOCA::PhaseTransition::ExtendedGroup::getGradient() const 
{
  std::cout << "ERROR: GRADIENT VECTOR NOT CALCULATEED IN TRAMONTO_GROUP!! " << std::endl;
  return *newtonVector;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getXPtr() const 
{
  return xVector;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getFPtr() const 
{  
  return fVector;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getNewtonPtr() const 
{
  return newtonVector;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getGradientPtr() const 
{
  std::cout << "ERROR: GRADIENT VECTOR NOT CALCULATEED IN TRAMONTO_GROUP!! " << std::endl;
  return newtonVector;
}


void LOCA::PhaseTransition::ExtendedGroup::print() const
{
/*
  std::cout << "x = " << *xVector << "\n";

  if (isValidF) {
    std::cout << "F(x) = " << *fVector << "\n";
    std::cout << "|| F(x) || = " << normF << "\n";
  }
  else
    std::cout << "F(x) has not been computed" << "\n";
*/
  
  std::cout << std::endl;
}

void  LOCA::PhaseTransition::ExtendedGroup::setParams(const LOCA::ParameterVector& p)
{ grp->setParams(p);}

void  LOCA::PhaseTransition::ExtendedGroup::setParam(std::string paramID, double val)
{ 
  resetIsValid();
  grp->setParam(paramID, val);
}

void  LOCA::PhaseTransition::ExtendedGroup::setParam(int paramID, double val)
{ grp->setParam(paramID, val); }

const LOCA::ParameterVector&  LOCA::PhaseTransition::ExtendedGroup::getParams() const
{ return grp->getParams(); }
double  LOCA::PhaseTransition::ExtendedGroup::getParam(std::string paramID) const
{  return grp->getParam(paramID); }
double  LOCA::PhaseTransition::ExtendedGroup::getParam(int paramID) const
{  return grp->getParam(paramID); }

void
LOCA::PhaseTransition::ExtendedGroup::setParamsMulti(
                          const std::vector<int>& paramIDs,
                          const NOX::Abstract::MultiVector::DenseMatrix& vals)
{

  grp->setParamsMulti(paramIDs, vals);
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (paramIDs[i] == bifParamID)
      grp->setParam(bifParamID, vals(i,0));
  }
}

NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::computeDfDpMulti(
                                            const std::vector<int>& paramIDs,
                                            NOX::Abstract::MultiVector& dfdp,
                                            bool isValid_F)
{
   std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  return  NOX::Abstract::Group::BadDependency;
}


void  LOCA::PhaseTransition::ExtendedGroup::printSolution(const NOX::Abstract::Vector& sol_,
      const double param) const
{ 

  const NOX::Abstract::Vector& sol =
    *((dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(sol_)).X1());
  grp->printSolution(sol, param);

  const NOX::Abstract::Vector& sol2 =
    *((dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(sol_)).X2());
  double ptp = (dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(sol_)).PTP();
  grp->printSolution(sol2, ptp);
}

void  LOCA::PhaseTransition::ExtendedGroup::printSolution(const double param) const
{
  printSolution(*xVector, param);
}
