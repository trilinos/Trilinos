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
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"

//#define FILE_DEBUG 1

//! Constructor
LOCA::Epetra::xyztPrec::
xyztPrec(EpetraExt::BlockCrsMatrix &jacobian_, 
	 Epetra_Vector &splitVec_, 
	 Epetra_RowMatrix &splitJac_,
	 Epetra_RowMatrix &splitMass_,
	 EpetraExt::BlockVector &solution_,
         EpetraExt::BlockVector& solutionOverlap_,
	 Epetra_Import &overlapImporter_,
	 NOX::Parameter::List &precPrintParams_, 
	 NOX::Parameter::List &precLSParams_, 
	 const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm_) :  
  jacobian(jacobian_),
  splitVec(0),
  splitRes(0),
  splitVecOld(0),
  splitJac(dynamic_cast<Epetra_CrsMatrix &>(splitJac_)),
  splitMass(dynamic_cast<Epetra_CrsMatrix &>(splitMass_)),
  solution(solution_),
  solutionOverlap(solutionOverlap_),
  overlapImporter(overlapImporter_),
  printParams(precPrintParams_), 
  lsParams(precLSParams_), 
  globalComm(globalComm_),
  linSys(),
  jacobianBlock(0),
  massBlock(0)
{

  string prec = lsParams.getParameter("XYZTPreconditioner","None");

  if (prec == "Global") {
    label = "LOCA::Epetra::xyztPrec::Global";
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
  else if (prec == "Sequential") {
    label = "LOCA::Epetra::xyztPrec::Sequential";
    cout << "LOCA::Epetra::xyztPrec = Sequential" << endl;
    
    // Create the single block jacobian
    jacobianBlock = new Epetra_CrsMatrix(splitJac);
    jacobianBlock->PutScalar(0.0);
    jacobian.ExtractBlock(*jacobianBlock, 0, 0);
    //cout << "jacobianBlock = \n" << *jacobianBlock << endl;

    // Create the single block jacobian
    massBlock = new Epetra_CrsMatrix(splitMass);
    massBlock->PutScalar(0.0);
    //cout << "jacobianBlock = \n" << *jacobianBlock << endl;

    // Create the Linear System
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = 
      Teuchos::rcp(&((NOX::Epetra::Interface::Required&)*this),false);
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = 
      Teuchos::rcp(&((NOX::Epetra::Interface::Jacobian&)*this),false);
    Teuchos::RefCountPtr<Epetra_Operator> A =
      Teuchos::rcp(jacobianBlock,false);
    linSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
							iReq, iJac, 
							A, solution));
    Asolve = Teuchos::rcp(jacobianBlock,false);
    splitVec = new Epetra_Vector(splitVec_);
    splitRes = new Epetra_Vector(splitVec_);
    splitVecOld = new Epetra_Vector(splitVec_);
  }
  else if (prec == "None") {
    label = "LOCA::Epetra::xyztPrec::None";
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
  if (jacobianBlock) delete jacobianBlock;
  if (massBlock) delete massBlock;
  if (splitVec) delete splitVec;
  if (splitRes) delete massBlock;
  if (splitVecOld) delete splitVecOld;

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
  string prec = lsParams.getParameter("XYZTPreconditioner","None");

  if (prec == "None") {
    return 0;
  }

  if (prec == "Global") {
    // convert multivectors to vectors
    Epetra_MultiVector& input_EMV = const_cast<Epetra_MultiVector&>(input);
    Epetra_Vector& input_EV = dynamic_cast<Epetra_Vector&>(input_EMV);
    Epetra_Vector& result_EV = dynamic_cast<Epetra_Vector&>(result);
    
    // use referenced counted pointers to vectors
    const Teuchos::RefCountPtr<Epetra_Vector> input_EV_RCP = Teuchos::rcp(&input_EV, false);
    Teuchos::RefCountPtr<Epetra_Vector> result_EV_RCP = Teuchos::rcp(&result_EV, false);
    
    // create views of vectors 
    const NOX::Epetra::Vector input_NEV(input_EV_RCP, NOX::Epetra::Vector::CreateView);
    NOX::Epetra::Vector result_NEV(result_EV_RCP, NOX::Epetra::Vector::CreateView);
    
    // apply preconditioner as specified in lsParams
    bool stat = linSys->applyRightPreconditioning(false, lsParams, input_NEV, result_NEV);

    return (stat) ? 0 : 1;
  }    
  else if (prec == "Sequential") {

    // convert multivectors to vectors
    Epetra_MultiVector& input_EMV = const_cast<Epetra_MultiVector&>(input);
    Epetra_Vector& input_EV = dynamic_cast<Epetra_Vector&>(input_EMV);
    Epetra_Vector& result_EV = dynamic_cast<Epetra_Vector&>(result);
    
    // residual needs to be in BlockVector and contain values from input
    EpetraExt::BlockVector residual(solution);    // inherit offsets from solution
    residual.Epetra_Vector::operator=(input_EV);  // populate with input values

#ifdef FILE_DEBUG
      cout << "residual = \n" << residual << endl;
      char fname[255];
#endif

    // Assumes one time domain currently

    for (int i=0; i < globalComm->NumTimeSteps(); i++) {
      // Extract jacobian block to use in solve
      //   diagonal is column 0 on first step, colum one on all other)
      if (globalComm->FirstTimeStepOnDomain() == 0 && i==0)
           jacobian.ExtractBlock(*jacobianBlock, 0, 0);
      else
           jacobian.ExtractBlock(*jacobianBlock, i, 1);
      linSys->setJacobianOperatorForSolve(Asolve);

#ifdef FILE_DEBUG
      sprintf(fname,"jacBlock_%d.m",i);
      EpetraExt::RowMatrixToMatlabFile(fname,*jacobianBlock);
      cout << "jacobianBlock = \n" << *jacobianBlock << endl;
#endif

      // get x and residual (r) corresponding to current block
      residual.ExtractBlockValues(*splitRes, jacobian.RowIndex(i));
      solution.ExtractBlockValues(*splitVec, jacobian.RowIndex(i));

#ifdef FILE_DEBUG
      // debugging output
      sprintf(fname,"splitRes_%d.m",i);
      EpetraExt::VectorToMatlabFile(fname,*splitRes);
      cout << "splitRes = \n" << *splitRes << endl;
      sprintf(fname,"splitVec_%d.m",i);
      EpetraExt::VectorToMatlabFile(fname,*splitVec);
      cout << "splitVec = \n" << *splitVec << endl;
#endif

      // Create a new preconditioner for the single block
      linSys->destroyPreconditioner();
      linSys->createPreconditioner(*splitVec, lsParams, false);

      int blockRowOld = jacobian.RowIndex(i) - 1;  //Hardwired for -1 in stencil
      if (blockRowOld >= 0)  {
	// update RHS with mass matrix * solution from previous time step
	solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	jacobian.ExtractBlock(*massBlock, i, 0);
#ifdef FILE_DEBUG
	sprintf(fname,"massBlock_%d.m",i);
	EpetraExt::RowMatrixToMatlabFile(fname,*massBlock);
#endif
	massBlock->Multiply(false, *splitVecOld, *splitVec);

	splitRes->Update(-1.0, *splitVec, 1.0);
      }	
      
      // use referenced counted pointers to vectors
      const Teuchos::RefCountPtr<Epetra_Vector> input_EV_RCP = Teuchos::rcp(splitRes, false);
      Teuchos::RefCountPtr<Epetra_Vector> result_EV_RCP = Teuchos::rcp(splitVec, false);
    
      // copied from GLOBAL prec, but only need to do once in constructor for SEQUENTIAL
      const NOX::Epetra::Vector input_NEV(input_EV_RCP, NOX::Epetra::Vector::CreateView);
      NOX::Epetra::Vector result_NEV(result_EV_RCP, NOX::Epetra::Vector::CreateView);
//result_NEV.init(0.0);  //Zero init guess might be better?

      // solve the problem
      bool stat = linSys->applyJacobianInverse(lsParams, input_NEV, result_NEV);
      solution.LoadBlockValues(result_NEV.getEpetraVector(), jacobian.RowIndex(i));

#ifdef FILE_DEBUG
      sprintf(fname,"splitResAfterSolve_%d.m",i);
      EpetraExt::VectorToMatlabFile(fname,*splitRes);
      cout << "splitResAfterSolve = \n" << *splitRes << endl;
#endif

    }
    
    result = solution;
    return 0;
  }    
  else {
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
  else if (prec == "Sequential") {
    cout << "Sequential Preconditioner" << endl;
    // The next two lines should be performed for each block in ApplyInverse
    //linSys->destroyPreconditioner();
    //linSys->createPreconditioner(x, lsParams, false);
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

