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
	 Epetra_CrsMatrix &splitJac_,
	 Epetra_CrsMatrix &splitMass_,
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
  splitJac(splitJac_),
  splitMass(splitMass_),
  solution(solution_),
  solutionOverlap(solutionOverlap_),
  overlapImporter(overlapImporter_),
  printParams(precPrintParams_), 
  lsParams(precLSParams_), 
  globalComm(globalComm_),
  linSys(std::vector<NOX::Epetra::LinearSystemAztecOO*>(globalComm_->NumTimeStepsOnDomain())),
  jacobianBlock(std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> >(globalComm_->NumTimeStepsOnDomain())),
  massBlock(std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> >(globalComm_->NumTimeStepsOnDomain()))
{

  string prec = lsParams.getParameter("XYZTPreconditioner","None");

  Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = 
    Teuchos::rcp(&((NOX::Epetra::Interface::Required&)*this),false);
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = 
    Teuchos::rcp(&((NOX::Epetra::Interface::Jacobian&)*this),false);

  if (prec == "Global") {
    label = "LOCA::Epetra::xyztPrec::Global";
    cout << "LOCA::Epetra::xyztPrec = Global" << endl;

    // Create the Linear System
    linSys[0] = new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
				iReq, iJac, Teuchos::rcp(&jacobian,false), solution);

    linSys[0]->setJacobianOperatorForSolve(Teuchos::rcp(&jacobian,false));
  }
  else if (prec == "SequentialLowMemory" || prec == "ParallelLowMemory" || 
           prec == "Sequential" || prec == "Parallel") {
    if (prec == "SequentialLowMemory") {
      label = "LOCA::Epetra::xyztPrec::SequentialLowMemory";
      cout << "LOCA::Epetra::xyztPrec = SequentialLowMemory" << endl;
    }
    else if (prec == "ParallelLowMemory") {
      label = "LOCA::Epetra::xyztPrec::ParallelLowMemory";
      cout << "LOCA::Epetra::xyztPrec = ParallelLowMemory" << endl;
    }
    if (prec == "Sequential") {
      label = "LOCA::Epetra::xyztPrec::Sequential";
      cout << "LOCA::Epetra::xyztPrec = Sequential" << endl;
    }
    else if (prec == "Parallel") {
      label = "LOCA::Epetra::xyztPrec::Parallel";
      cout << "LOCA::Epetra::xyztPrec = Parallel" << endl;
    }
    
    // Create temporary space for Epetra vectors and NOX vies of same space
    splitVecOld = new Epetra_Vector(splitJac.RowMap());
    splitVec = new Epetra_Vector(splitJac.RowMap());
    splitRes = new Epetra_Vector(splitJac.RowMap());
    splitRes_NEV = new NOX::Epetra::Vector(Teuchos::rcp(splitRes, false),
		                           NOX::Epetra::Vector::CreateView);
    splitVec_NEV = new NOX::Epetra::Vector(Teuchos::rcp(splitVec, false),
		                           NOX::Epetra::Vector::CreateView);

    // Create the Linear System
    int imax = 1;
    if (prec == "Sequential" || prec == "Parallel") 
       imax = globalComm->NumTimeStepsOnDomain();

    for (int i=0; i < imax; i++) {
      jacobianBlock[i] = Teuchos::rcp(new Epetra_CrsMatrix(splitJac));
      jacobianBlock[i]->PutScalar(0.0);
      massBlock[i] = Teuchos::rcp(new Epetra_CrsMatrix(splitMass));
      massBlock[i]->PutScalar(0.0);
      //jacobian.ExtractBlock(*jacobianBlock[i], 0, 0);
      linSys[i] = new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
			iReq, iJac, jacobianBlock[i], solution);
    }
	
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
  if (splitVec) delete splitVec;
  if (splitRes) delete splitRes;
  if (splitVecOld) delete splitVecOld;
  if (splitRes_NEV) delete splitRes_NEV;
  if (splitVec_NEV) delete splitVec_NEV;

  string prec = lsParams.getParameter("XYZTPreconditioner","None");

  int imax=1;
  if (prec == "Sequential" || prec == "Parallel") 
     imax =  globalComm->NumTimeStepsOnDomain();

  for (int i=0; i < imax; i++) {
    linSys[i]->destroyPreconditioner();
    delete linSys[i];
  }
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
    // create NOX vectors as views of input and results vectors 
    const NOX::Epetra::Vector input_NEV(Teuchos::rcp( const_cast<Epetra_Vector*>(input(0)), false),
		                        NOX::Epetra::Vector::CreateView);
    NOX::Epetra::Vector result_NEV(Teuchos::rcp(result(0), false), NOX::Epetra::Vector::CreateView);
    
    // apply preconditioner as specified in lsParams
    bool stat = linSys[0]->applyRightPreconditioning(false, lsParams, input_NEV, result_NEV);

    //    return (stat) ? 0 : 1;
  }    
  else if (prec == "SequentialLowMemory" || prec == "ParallelLowMemory") {
    
    // residual needs to be in BlockVector and contain values from input
    EpetraExt::BlockVector residual(solution);    // inherit offsets from solution
    residual.Epetra_Vector::operator=(*input(0));  // populate with input values
    
    // Set up loop logic for SequentialLowMemory (block Gauss Seydel over time domains)
    // or parallel (Jacobi by time domains)
    int sequentialDomains, innerCheck;
    if (prec=="SequentialLowMemory") {
      sequentialDomains = globalComm->NumSubDomains();
      innerCheck = globalComm->SubDomainRank();
    }
    else if (prec=="ParallelLowMemory") {
      sequentialDomains = 1;
      innerCheck = 0;  
      //In General:  innerCheck = globalComm->SubDomainRank() % sequentialDomains;
    }
    
    for (int isd=0; isd < sequentialDomains; isd++) {
      
      // Communicate data from other time domains into overlapped vector
      // This serves as a barrier as well.
      solutionOverlap.Import(solution, overlapImporter, Insert);
      
      //Work only on active time domain
      if ( isd == innerCheck ) {
	
	for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {
	  bool isFirstGlobalTimeStep = (globalComm->FirstTimeStepOnDomain() + i == 0);
	  // Extract jacobian block to use in solve
	  //   diagonal is column 0 on first step, colum one on all other)
	  if (isFirstGlobalTimeStep)
	    jacobian.ExtractBlock(*jacobianBlock[0], 0, 0);
	  else
	    jacobian.ExtractBlock(*jacobianBlock[0], i, 1);
	  linSys[0]->setJacobianOperatorForSolve(jacobianBlock[0]);
	  
	  // get x and residual (r) corresponding to current block
	  residual.ExtractBlockValues(*splitRes, jacobian.RowIndex(i));
	  solution.ExtractBlockValues(*splitVec, jacobian.RowIndex(i));

	  // Create a new preconditioner for the single block
	  linSys[0]->destroyPreconditioner();
	  linSys[0]->createPreconditioner(*splitVec, lsParams, false);
	  
	  if (!isFirstGlobalTimeStep)  {
	    // update RHS with mass matrix * solution from previous time step 
	    // (which will be from previous time domain if i==0)
	    if (i==0) 
	      solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    else      
	      solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    jacobian.ExtractBlock(*massBlock[0], i, 0);
	    massBlock[0]->Multiply(false, *splitVecOld, *splitVec);
	    splitRes->Update(-1.0, *splitVec, 1.0);
	  }	
	  
	  // solve the problem, and put solution into slution vector
	  bool stat = linSys[0]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
	  solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
	  
	}
      }
    }
    
    result = solution;
  }    
  else if (prec == "Sequential" || prec == "Parallel") {
    
    // residual needs to be in BlockVector and contain values from input
    EpetraExt::BlockVector residual(solution);    // inherit offsets from solution
    residual.Epetra_Vector::operator=(*input(0));  // populate with input values
    
    // Set up loop logic for SequentialLowMemory (block Gauss Seydel over time domains)
    // or parallel (Jacobi by time domains)
    int sequentialDomains, innerCheck;
    if (prec=="Sequential") {
      sequentialDomains = globalComm->NumSubDomains();
      innerCheck = globalComm->SubDomainRank();
    }
    else if (prec=="Parallel") {
      sequentialDomains = 1;
      innerCheck = 0;  
      //In General:  innerCheck = globalComm->SubDomainRank() % sequentialDomains;
    }
    
    for (int isd=0; isd < sequentialDomains; isd++) {
      
      // Communicate data from other time domains into overlapped vector
      // This serves as a barrier as well.
      solutionOverlap.Import(solution, overlapImporter, Insert);
      
      //Work only on active time domain
      if ( isd == innerCheck ) {
	
	for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {

	  bool isFirstGlobalTimeStep = (globalComm->FirstTimeStepOnDomain() + i == 0);

	  // This line is needed or else Aztec error -11
	  linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
	  
	  // get x and residual (r) corresponding to current block
	  residual.ExtractBlockValues(*splitRes, jacobian.RowIndex(i));
	  solution.ExtractBlockValues(*splitVec, jacobian.RowIndex(i));

	  if (!isFirstGlobalTimeStep)  {
	    // update RHS with mass matrix * solution from previous time step 
	    // (which will be from previous time domain if i==0)
	    if (i==0) 
	      solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    else      
	      solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    massBlock[i]->Multiply(false, *splitVecOld, *splitVec);
	    splitRes->Update(-1.0, *splitVec, 1.0);
	  }	
	  
	  // solve the problem, and put solution into slution vector
	  bool stat = linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
	  solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
	}
      }
    }
    
    result = solution;
  }    

  return 0;
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
    linSys[0]->destroyPreconditioner();
    linSys[0]->createPreconditioner(x, lsParams, false);
    return true;
  }    
  else if (prec == "SequentialLowMemory") {
    cout << "SequentialLowMemory Preconditioner" << endl;
    return true;
  }    
  else if (prec == "ParallelLowMemory") {
    cout << "ParallelLowMemory Preconditioner" << endl;
    return true;
  }    
  else if (prec == "Sequential"  || prec == "Parallel") {
    for (int i=0; i< globalComm->NumTimeStepsOnDomain(); i++ ) {
      if (globalComm->FirstTimeStepOnDomain() + i == 0)
        jacobian.ExtractBlock(*jacobianBlock[i], 0, 0);
      else {
        jacobian.ExtractBlock(*jacobianBlock[i], i, 1);
        jacobian.ExtractBlock(*massBlock[i], i, 0);
      }
      //is following line needed???
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
      linSys[i]->destroyPreconditioner();
      linSys[i]->createPreconditioner(x, lsParams, false);
    }
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

