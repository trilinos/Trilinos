// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface.H"        // class data members
#include "NOX_Parameter_List.H"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "NOX_Epetra_SharedOperator.H"
#include "Epetra_Vector.h"
#include "LOCA_Utils.H"

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J) :
  NOX::Epetra::Group(printParams, par, i, x, J),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  eigenvalCounter(0)
{
}

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J, 
			   Epetra_Operator& M) :
  NOX::Epetra::Group(printParams, par, i, x, J, M),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  eigenvalCounter(0)
{
}

LOCA::Epetra::Group::Group(const LOCA::Epetra::Group& source, 
			   NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  params(source.params),
  userInterface(source.userInterface),
  tmpVectorPtr2(0),
  eigenvalCounter(source.eigenvalCounter)
{
}

LOCA::Epetra::Group::~Group() 
{
  delete tmpVectorPtr2;
}

NOX::Abstract::Group* 
LOCA::Epetra::Group::clone(NOX::CopyType type) const 
{
  return new Group(*this, type);
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Epetra::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Epetra::Group& source)
{
  params = source.params;
  NOX::Epetra::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);
  return *this;
}

void 
LOCA::Epetra::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::Epetra::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::Epetra::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeF() 
{
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);
  
  return NOX::Epetra::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::computeJacobian() 
{
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);

  return NOX::Epetra::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::Epetra::Group::getParams() const 
{
  return params;
}

NOX::Epetra::Interface& 
LOCA::Epetra::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::Epetra::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector& x_,
                                   const double conParam) const
{
  userInterface.printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::augmentJacobianForHomotopy(double conParamValue)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == 0)
    tmpVectorPtr = new Epetra_Vector(xVector.getEpetraVector());
  if (tmpVectorPtr2 == 0)
    tmpVectorPtr2 = new Epetra_Vector(xVector.getEpetraVector());

  tmpVectorPtr2->PutScalar(1.0-conParamValue);

  // See if it is an Epetra_CrsMatrix
  Epetra_CrsMatrix* testCrs = 0;
  testCrs = dynamic_cast<Epetra_CrsMatrix*>
            (&(sharedJacobian.getOperator(this)));
  if (testCrs != 0) {

    testCrs->Scale(conParamValue);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  Epetra_VbrMatrix* testVbr = 0;
  testVbr = dynamic_cast<Epetra_VbrMatrix*>(&(sharedJacobian.getOperator(this)));
  if (testVbr != 0) {
    
    testVbr->Scale(conParamValue);
    testVbr->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testVbr->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;
  }

  // Otherwise this alg won't work
  cout << "ERROR: LOCA::Epetra::Group::augmentJacobianForHomotopy() - "
       << "the Jacobian must be either an Epetra_CrsMatrix or an "
       << "Epetra_VbrMatrix!" << endl;
  throw "LOCA Error";

  return LOCA::Abstract::Group::Ok;
}

#ifdef HAVE_LOCA_ANASAZI
#include "AnasaziLOCAInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#endif

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeEigenvalues(NOX::Parameter::List& params)
{
#ifdef HAVE_LOCA_ANASAZI

  NOX::Parameter::List& aList = params.sublist("LOCA").sublist("Stepper").sublist("Anasazi");
  // The following are Parameter for the Anasazi Eigensolver
  int blksz =   aList.getParameter("Block Size", 1);       //  The block size
  int length =  aList.getParameter("Arnoldi Size", 30);      //  The maximum length of the Arnoldi factorization
  int nev =     aList.getParameter("NEV", 4);         //  The number of requested eigenvalues
  double tol =  aList.getParameter("Tol", 1.0e-7); //  Tolerance for the converged eigenvalues
  int step =    aList.getParameter("Convergence Check", 1);        //  This checks convergence every so many steps
  int restart = aList.getParameter("Restarts",1);    //  This is the number of restarts allowed
  int freq =    aList.getParameter("Frequency",1);    // How often to recalculate eigenvalues
  int debug =   aList.getParameter("Debug Level",1);  // Anasazi Debug level
  string which="LM";   //  Which eigenvalues are of interest.

  // Check if eigenvalues are requested this continuation step
  if (eigenvalCounter++%freq != 0) {
    if (Utils::doPrint(Utils::StepperIteration)) 
      cout <<"\tAnasazi Eigensolver not requested this continuation step." << endl;
    return LOCA::Abstract::Group::Ok;
  }
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\n" << Utils::fill(64,'=')
         << "\nAnasazi Eigensolver starting with block size " << blksz
	 << "\n" << endl;
  }

  // Create updated Jacobian matrix
  computeJacobian();

  // Create the operator and initial vector
  AnasaziLOCAMat<double> Amat( params, *this );
  AnasaziLOCAVec<double> ivec( xVector, blksz );
  ivec.MvRandom();

  // Initialize the solver
  Anasazi::BlockArnoldi<double> LOCABlockArnoldi( Amat, ivec, tol, nev, length,
                                      blksz, which, step, restart );

  // Print out debugging information on single proc
  LOCABlockArnoldi.setDebugLevel(debug);

  // Solve the problem to the specified tolerance
  LOCABlockArnoldi.solve();

  // Look at the solutions once if debug=0
  if (Utils::doPrint(Utils::StepperIteration)) 
    if (debug == 0) LOCABlockArnoldi.currentStatus();

  // Obtain the eigenvalues / eigenvectors
  double * evalr = LOCABlockArnoldi.getEvals();
  double * evali = LOCABlockArnoldi.getiEvals();

  if (Utils::doPrint(Utils::StepperIteration)) {
    cout<<"Untransformed eigenvalues (since the operator was the Jacobian inverse)"<<endl;
    cout.precision(7);
    for (int i=0; i<nev; i++) {
      double mag=evalr[i]*evalr[i]+evali[i]*evali[i];
      cout<<"Eigenvalue "<<i<<" : "<<evalr[i]/mag<<"  "<<-evali[i]/mag<<" i"<<endl;
    }
  }

  /* Comment out Eigenvector extraction for now
  AnasaziLOCAVec<double> evecR( xVector, nev );
  LOCABlockArnoldi.getEvecs( evecR );
  AnasaziLOCAVec<double> evecI( xVector, nev );
  LOCABlockArnoldi.getiEvecs( evecI );
  */

  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\nAnasazi Eigensolver finished.\n" 
         << Utils::fill(64,'=') << "\n" << endl;
  }

  return LOCA::Abstract::Group::Ok;
#else
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\nWarning: LOCA::Epetra::Group::computeEigenvalues:\n\t"
         <<  "Anasazi Eigensolver requested but not compiled in!" << endl;
  }
  return LOCA::Abstract::Group::Ok;
#endif
}

