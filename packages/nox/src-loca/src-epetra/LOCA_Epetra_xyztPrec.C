//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_xyztPrec.H"

//! Constructor
LOCA::Epetra::xyztPrec::
xyztPrec(EpetraExt::BlockCrsMatrix &jacobian_, 
	 EpetraExt::BlockVector &solution_,
	 NOX::Parameter::List &precPrintParams_, 
	 NOX::Parameter::List &precLSParams_, 
	 const Teuchos::RefCountPtr<Epetra_Comm> globalComm_) :  
  jacobian(jacobian_),
  solution(solution_),
  printParams(precPrintParams_), 
  lsParams(precLSParams_), 
  globalComm(globalComm_),
  linSys(),
  label("LOCA::Epetra::xyztPrec")
{

  string prec = lsParams.getParameter("XYZTPreconditioner","None");
  if (prec == "Global") {
    cout << "LOCA::Epetra::xyztPrec = Global" << endl;
    // Create the Linear System
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = 
      Teuchos::rcp(&((NOX::Epetra::Interface::Required&)*this),false);
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = 
      Teuchos::rcp(&((NOX::Epetra::Interface::Jacobian&)*this),false);
    Teuchos::RefCountPtr<Epetra_Operator> A =
      Teuchos::rcp(&jacobian,false);
    linSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
							iReq, iJac, 
							A, solution));

    Teuchos::RefCountPtr<const Epetra_Operator> Asolve =
      Teuchos::rcp(&jacobian,false);
    linSys->setJacobianOperatorForSolve(Asolve);
  }
  else if (prec == "None") {
    cout << "LOCA::Epetra::xyztPrec = None" << endl;
  }
  else {
    string errorMessage = "Option for \"Preconditioner\" is invalid!";
    throwError("LOCA::Epetra::xyztPrec::xyztPrec", errorMessage);
  }

  cout << "LOCA::Epetra::xyztPrec constructor completed." << endl;
}

LOCA::Epetra::xyztPrec::
~xyztPrec()
{
}


int LOCA::Epetra::xyztPrec::
SetUseTranspose(bool UseTranspose)
{
  // Disable this option
  return false;
}

int LOCA::Epetra::xyztPrec::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Not implemented: Throw an error!
  cout << "ERROR: LOCA::Epetra::xyztPrec::Apply() - "
       << "method is NOT implemented!!  " << endl;
  throw "LOCA Error";
  return false;
}

int LOCA::Epetra::xyztPrec::
ApplyInverse(const Epetra_MultiVector& input,
	     Epetra_MultiVector& result) const
{
  cout << "LOCA::Epetra::xyztPrec::ApplyInverse called - ";

  string prec = lsParams.getParameter("XYZTPreconditioner","None");

  if (prec == "Global") {
    cout << "Global Preconditioner" << endl;

    Epetra_MultiVector& input_EMV = const_cast<Epetra_MultiVector&>(input);
    Epetra_Vector& input_EV = dynamic_cast<Epetra_Vector&>(input_EMV);
    Epetra_Vector& result_EV = dynamic_cast<Epetra_Vector&>(result);

    const Teuchos::RefCountPtr<Epetra_Vector> input_EV_RCP = Teuchos::rcp(&input_EV, false);
    Teuchos::RefCountPtr<Epetra_Vector> result_EV_RCP = Teuchos::rcp(&result_EV, false);

    const NOX::Epetra::Vector* input_NEV = 
      new const NOX::Epetra::Vector(input_EV_RCP, NOX::Epetra::Vector::CreateView, NOX::DeepCopy);
    NOX::Epetra::Vector* result_NEV = 
      new NOX::Epetra::Vector(result_EV_RCP, NOX::Epetra::Vector::CreateView, NOX::DeepCopy);

#if 0 
    bool stat = linSys->applyJacobianInverse(lsParams, *input_NEV, *result_NEV);
    cout << "LOCA::Epetra::xyztPrec::ApplyInverse - linSys->applyJacobianInverse = " << stat << endl;
    return (stat) ? 0 : 1;
#else
    //cout << "result before linSys->applyRightPreconditioning:\n" << endl;
    //result_NEV->print(cout);
    bool stat = linSys->applyRightPreconditioning(false,lsParams,*input_NEV,*result_NEV);
    cout << "result after  linSys->applyRightPreconditioning:\n" << endl;
    //result_NEV->print(cout);
    //cout << "LOCA::Epetra::xyztPrec::ApplyInverse - linSys->applyRightPreconditioning = " << stat << endl;
    return (stat) ? 0 : 1;
#endif
  }    
  else {
    cout << "No Preconditioner" << endl;
    // Do nothing
    return 0;
  }
}

double LOCA::Epetra::xyztPrec::
NormInf() const
{
  // Not implemented: Throw an error!
  cout << "ERROR: LOCA::Epetra::xyztPrec::NormInf() - "
       << "method is NOT implemented!!  " << endl;
  throw "LOCA Error";
  return 0.0;
}

bool LOCA::Epetra::xyztPrec::
UseTranspose() const
{
  // Not implemented: Throw an error!
  cout << "ERROR: LOCA::Epetra::xyztPrec::UseTranspose() - "
       << "method is NOT implemented!!  " << endl;
  throw "LOCA Error";
  return false;
}

bool LOCA::Epetra::xyztPrec::
HasNormInf() const
{
  // NormInf is not implemented
  return false;
}

const char* LOCA::Epetra::xyztPrec::
Label() const
{
   return label.c_str();
}

const Epetra_Comm& LOCA::Epetra::xyztPrec::
Comm() const
{
  return jacobian.Comm();
}

const Epetra_Map& LOCA::Epetra::xyztPrec::
OperatorDomainMap () const
{
  return jacobian.OperatorDomainMap();
}

const Epetra_Map& LOCA::Epetra::xyztPrec::
OperatorRangeMap () const
{
  return jacobian.OperatorRangeMap();
}

bool LOCA::Epetra::xyztPrec::
computeF(const Epetra_Vector&, Epetra_Vector&, 
	 NOX::Epetra::Interface::Required::FillType)
{
  cout << "LOCA::Epetra::xyztPrec::computeF called." << endl;
  return false;
}

bool LOCA::Epetra::xyztPrec::
computeJacobian(const Epetra_Vector&, Epetra_Operator&)
{
  cout << "LOCA::Epetra::xyztPrec::computeJacobian called." << endl;
  return false;
}

bool LOCA::Epetra::xyztPrec::
computePreconditioner(const Epetra_Vector& x,
		      Epetra_Operator& Prec,
		      NOX::Parameter::List* p)
{
  cout << "LOCA::Epetra::xyztPrec::computePreconditioner called - ";
  string prec = lsParams.getParameter("XYZTPreconditioner","None");
  if (prec == "Global") {
    cout << "Global Preconditioner" << endl;
    linSys->destroyPreconditioner();
    linSys->createPreconditioner(x, lsParams, false);
    return true;
  }    
  else {
    cout << "No Preconditioner" << endl;
    return true;
  }
}

void LOCA::Epetra::xyztPrec::
throwError(const string& functionName, const string& errorMsg) const
{
  cout << "LOCA::Epetra::xyztPrec" << functionName 
	 << " - " << errorMsg << endl;
  throw "LOCA Error";
}

