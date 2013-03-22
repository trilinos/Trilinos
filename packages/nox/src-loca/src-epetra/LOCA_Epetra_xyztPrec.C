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

#include "LOCA_Epetra_xyztPrec.H"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_Map.h"

//#define FILE_DEBUG 1

//! Constructor
LOCA::Epetra::xyztPrec::
xyztPrec(EpetraExt::BlockCrsMatrix &jacobian_, 
	 Epetra_CrsMatrix &splitJac_,
	 EpetraExt::BlockVector &solution_,
         EpetraExt::BlockVector &solutionOverlap_,
	 Epetra_Import &overlapImporter_,
	 Teuchos::ParameterList &precPrintParams_, 
	 Teuchos::ParameterList &precLSParams_, 
	 const Teuchos::RCP<EpetraExt::MultiComm> globalComm_) :  
  jacobian(jacobian_),
  splitJac(splitJac_),
  solution(solution_),
  solutionOverlap(solutionOverlap_),
  overlapImporter(overlapImporter_),
  printParams(precPrintParams_), 
  lsParams(precLSParams_), 
  globalComm(globalComm_),
  linSys(std::vector<NOX::Epetra::LinearSystemAztecOO*>(globalComm_->NumTimeStepsOnDomain())),
  jacobianBlock(std::vector<Teuchos::RCP<Epetra_CrsMatrix> >(1 + globalComm_->NumTimeStepsOnDomain())),
  massBlock(std::vector<Teuchos::RCP<Epetra_CrsMatrix> >(1 + globalComm_->NumTimeStepsOnDomain())),
  diagBlockSubdiag(std::vector<Teuchos::RCP<Epetra_Vector> >(globalComm_->NumTimeStepsOnDomain())),
  isPeriodic(precLSParams_.get("Periodic",false))
{

  std::string prec = lsParams.get("XYZTPreconditioner","None");

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = 
    Teuchos::rcp(&((NOX::Epetra::Interface::Required&)*this),false);
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = 
    Teuchos::rcp(&((NOX::Epetra::Interface::Jacobian&)*this),false);

  if (prec == "Global") {
    label = "LOCA::Epetra::xyztPrec::Global";
    // TODO: pass in globalData and use output stream
    //cout << "LOCA::Epetra::xyztPrec = Global" << std::endl;

    // Create the Linear System
    linSys[0] = new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
				iReq, iJac, Teuchos::rcp(&jacobian,false), solution);

    linSys[0]->setJacobianOperatorForSolve(Teuchos::rcp(&jacobian,false));
  }
  else if ( prec == "Sequential" || prec == "Parallel" ||
            prec == "BlockDiagonal" || prec == "Parareal" || prec == "BDSDT") {
    if (prec == "Sequential") {
      label = "LOCA::Epetra::xyztPrec::Sequential";
      // TODO: pass in globalData and use output stream
      //cout << "LOCA::Epetra::xyztPrec = Sequential" << std::endl;
    }
    else if (prec == "Parallel") {
      label = "LOCA::Epetra::xyztPrec::Parallel";
      // TODO: pass in globalData and use output stream
      //cout << "LOCA::Epetra::xyztPrec = Parallel" << std::endl;
    }
    else if (prec == "BlockDiagonal") {
      label = "LOCA::Epetra::xyztPrec::BlockDiagonal";
      // TODO: pass in globalData and use output stream
      //cout << "LOCA::Epetra::xyztPrec = BlockDiagonal" << std::endl;
    }
    else if (prec == "BDSDT") {
      label = "LOCA::Epetra::xyztPrec::BDSDT";
      // TODO: pass in globalData and use output stream
      //cout << "LOCA::Epetra::xyztPrec = BDSDT" << std::endl;
    }
    else if (prec == "Parareal") {
      label = "LOCA::Epetra::xyztPrec::Parareal";
      // TODO: pass in globalData and use output stream
      //cout << "LOCA::Epetra::xyztPrec = Parareal" << std::endl;
    }
    
    // Create temporary space for Epetra vectors and NOX vies of same space
    splitVecOld = Teuchos::rcp(new Epetra_Vector(splitJac.RowMap()));
    splitVec = Teuchos::rcp(new Epetra_Vector(splitJac.RowMap()));
    splitRes = Teuchos::rcp(new Epetra_Vector(splitJac.RowMap()));
    splitRes_NEV = Teuchos::rcp(new NOX::Epetra::Vector(splitRes,
		                           NOX::Epetra::Vector::CreateView));
    splitVec_NEV = Teuchos::rcp(new NOX::Epetra::Vector(splitVec,
		                           NOX::Epetra::Vector::CreateView));
    residual = Teuchos::rcp(new EpetraExt::BlockVector(solution));

    // Create the Linear System
    int imax = globalComm->NumTimeStepsOnDomain();
    if (prec == "Parareal") imax++;

    for (int i=0; i < imax; i++) {
      jacobianBlock[i] = Teuchos::rcp(new Epetra_CrsMatrix(splitJac));
      jacobianBlock[i]->PutScalar(0.0);
      if (prec != "BlockDiagonal") {
        massBlock[i] = Teuchos::rcp(new Epetra_CrsMatrix(splitJac));
        massBlock[i]->PutScalar(0.0);
      }
      linSys[i] = new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
			iReq, iJac, jacobianBlock[i], solution);
      if (prec == "BDSDT") {
        diagBlockSubdiag[i] = Teuchos::rcp(new Epetra_Vector(*splitVec));
      }
    }
	
  }
  else if (prec == "None") {
    label = "LOCA::Epetra::xyztPrec::None";
    // TODO: pass in globalData and use output stream
    //cout << "LOCA::Epetra::xyztPrec = None" << std::endl;
  }
  else {
    std::string errorMessage = "Option for \"Preconditioner\" is invalid!";
    throwError("LOCA::Epetra::xyztPrec::xyztPrec", errorMessage);
  }
  // TODO: pass in globalData and use output stream
  //cout << "LOCA::Epetra::xyztPrec constructor completed." << std::endl;
}

LOCA::Epetra::xyztPrec::
~xyztPrec()
{
  std::string prec = lsParams.get("XYZTPreconditioner","None");

  int imax=0;
  if (prec == "Global") 
     imax =  1;
  else if (prec == "Sequential" || prec == "Parallel" || prec == "BlockDiagonal" || prec == "BDSDT") 
     imax =  globalComm->NumTimeStepsOnDomain();
  else if (prec == "Parareal") 
     imax =  1 + globalComm->NumTimeStepsOnDomain();

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
  std::cout << "ERROR: LOCA::Epetra::xyztPrec::Apply() - "
       << "method is NOT implemented!!  " << std::endl;
  throw "LOCA Error";
  return false;
}

int LOCA::Epetra::xyztPrec::
ApplyInverse(const Epetra_MultiVector& input,
	     Epetra_MultiVector& result) const
{

  std::string prec = lsParams.get("XYZTPreconditioner","None");

  if (prec == "None") {
    return 0;
  }
  
  if (prec == "Global") {
    // create NOX vectors as views of input and results vectors 
    const NOX::Epetra::Vector input_NEV(Teuchos::rcp( const_cast<Epetra_Vector*>(input(0)), false),
		                        NOX::Epetra::Vector::CreateView);
    NOX::Epetra::Vector result_NEV(Teuchos::rcp(result(0), false), NOX::Epetra::Vector::CreateView);
    
    // apply preconditioner as specified in lsParams
    linSys[0]->applyRightPreconditioning(false, lsParams, input_NEV, result_NEV);

    //    return (stat) ? 0 : 1;
  }    
  else if (prec == "Sequential") {
    
    // residual needs to be in BlockVector and contain values from input
    solution.PutScalar(0.0);
    residual->Epetra_Vector::operator=(*input(0));  // populate with input values

    for (int isd=0; isd < globalComm->NumSubDomains(); isd++) {
      
      // Communicate data from other time domains into overlapped vector
      // This serves as a barrier as well.
      solutionOverlap.Import(solution, overlapImporter, Insert);
      
      //Work only on active time domain
      if ( isd == globalComm->SubDomainRank() ) {
	
	for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {

	  bool isFirstGlobalTimeStep = (globalComm->FirstTimeStepOnDomain() + i == 0);

	  // This line is needed or else Aztec error -11
	  linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
	  
	  // get x and residual (r) corresponding to current block
	  residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

	  if (!isFirstGlobalTimeStep)  {
	    // update RHS with mass matrix * solution from previous time step 
	    // (which will be from previous time domain if i==0)
	    //   NOTE:  "- 1" can be changed to jacobian.Stencil(i)[0]  
	    if (i==0) 
	      solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    else      
	      solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    massBlock[i]->Multiply(false, *splitVecOld, *splitVec);
	    splitRes->Update(-1.0, *splitVec, 1.0);
	  }	
	  
	  // solve the problem, and put solution into slution vector
          splitVec->PutScalar(0.0);
	  linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
	  solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
	}
      }
    }

    result = solution;
  }    
  else if (prec == "Parallel") {
    
    // residual needs to be in BlockVector and contain values from input
    solution.PutScalar(0.0);
    residual->Epetra_Vector::operator=(*input(0));  // populate with input values
    
    for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {
    
      // This line is needed or else Aztec error -11
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
  
      // get x and residual (r) corresponding to current block
      residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

      if (i > 0)  {
        // update RHS with mass matrix * solution from previous time step 
        // First step on each time domain doesn't feel previous time domain
        solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
        massBlock[i]->Multiply(false, *splitVecOld, *splitVec);
        splitRes->Update(-1.0, *splitVec, 1.0);
      }	
  
      // solve the problem, and put solution into slution vector
      splitVec->PutScalar(0.0);
      linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
      solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
    }
   
    result = solution;
  }    
  else if (prec == "BlockDiagonal") {
    
    // residual needs to be in BlockVector and contain values from input
    solution.PutScalar(0.0);
    residual->Epetra_Vector::operator=(*input(0));  // populate with input values
    
    for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {
    
      // This line is needed or else Aztec error -11
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
  
      // get x and residual (r) corresponding to current block
      residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

      // solve the problem, and put solution into slution vector
      splitVec->PutScalar(0.0);
      linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
      solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
    }
   
    result = solution;
  }    
  else if (prec == "BDSDT") {
    
    // 1. Sequential solve with I on diag blocks and diagonal matrix on subdiagonal
    
    // residual needs to be in BlockVector and contain values from input
    solution.PutScalar(0.0);
    residual->Epetra_Vector::operator=(*input(0));  // populate with input values

    for (int isd=0; isd < globalComm->NumSubDomains(); isd++) {
      
      // Communicate data from other time domains into overlapped vector
      // This serves as a barrier as well.
      solutionOverlap.Import(solution, overlapImporter, Insert);
      
      //Work only on active time domain
      if ( isd == globalComm->SubDomainRank() ) {
	
	for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {

	  bool isFirstGlobalTimeStep = (globalComm->FirstTimeStepOnDomain() + i == 0);

	  // get x and residual (r) corresponding to current block
	  residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

	  if (!isFirstGlobalTimeStep)  {
	    // update RHS with mass matrix * solution from previous time step 
	    // (which will be from previous time domain if i==0)
	    //   NOTE:  "- 1" can be changed to jacobian.Stencil(i)[0]  
	    if (i==0) 
	      solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    else      
	      solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	    splitRes->Multiply(-1.0, *diagBlockSubdiag[i], *splitVecOld, 1.0);
	  }	
	  
	  // Update this block of the solution vector
	  solution.LoadBlockValues(*splitRes, jacobian.RowIndex(i));
	}
      }
    }

    // Move result of first preconditioning step into initial vector of next

    // 2. BlockDiagonal Preconsitioner witrh this modified vector
    
    for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {
    
      // This line is needed or else Aztec error -11
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
  
      // get x and residual (r) corresponding to current block
      // updated residual in "solution" after first prec step
      solution.ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

      // solve the problem, and put solution into slution vector
      splitVec->PutScalar(0.0);
      linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
      solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
    }
   
    result = solution;
  }    
  else if (prec == "Parareal") {
    
    // residual needs to be in BlockVector and contain values from input
    solution.PutScalar(0.0);
    residual->Epetra_Vector::operator=(*input(0));  // populate with input values
    
    int N = globalComm->NumTimeStepsOnDomain();

    // Do sequential over all time domains, just 1 step per domain
    for (int isd=0; isd < globalComm->NumSubDomains(); isd++) {
      
      // Communicate data from other time domains into overlapped vector
      // This serves as a barrier as well.
      if (isd>0) solutionOverlap.Import(solution, overlapImporter, Insert);
      
      //Work only on active time domain, and skip the last domain if not periodic
      if ( isd == globalComm->SubDomainRank() && 
           ((isd != globalComm->NumSubDomains()-1) && !isPeriodic)) {
	
	  // Get extra matrix for 1 big time step per time domain
	  linSys[N]->setJacobianOperatorForSolve(jacobianBlock[N]);
	  
	  // get and residual for last time step on this domain
	  residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(N-1));

	  if (isd > 0)  {
	    // Get solution vector from end of previous time domain,
	    // Don't do this for first domain, even if periodic, because 
	    //   we are doing a xssequential sweep
	    solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(0) - 1));
	    massBlock[N]->Multiply(false, *splitVecOld, *splitVec);
	    splitRes->Update(-1.0, *splitVec, 1.0);
	  }	
	  
	  // solve the problem, and put solution into last solution vector on this domain
          splitVec->PutScalar(0.0);
	  linSys[N]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
	  solution.LoadBlockValues(*splitVec, jacobian.RowIndex(N-1));
	}
    }
    // Final communication of big Sequential steps
    solutionOverlap.Import(solution, overlapImporter, Insert);

    // Now launch into parallel algorithm, but with non-zero initial guess
    
    for (int i=0; i < globalComm->NumTimeStepsOnDomain(); i++) {
    
      // This line is needed or else Aztec error -11
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
  
      // get x and residual (r) corresponding to current block
      residual->ExtractBlockValues(*splitRes, jacobian.RowIndex(i));

      bool isFirstGlobalTimeStep = (globalComm->FirstTimeStepOnDomain() + i == 0);
      if (!isFirstGlobalTimeStep || isPeriodic)  {
        // update RHS with mass matrix * solution from previous time step 
        // (which will be from previous time domain if i==0)
        if (i==0) {
          //solutionOverlap.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
	  int offsetToPrevSolution = jacobian.Stencil(i)[0];
cout << "In Parareal, OFFSET  " << offsetToPrevSolution << " for time domain " <<  globalComm->SubDomainRank() << std::endl;
          solutionOverlap.ExtractBlockValues(*splitVecOld, jacobian.RowIndex(i) + offsetToPrevSolution);
	}
        else      
          solution.ExtractBlockValues(*splitVecOld, (jacobian.RowIndex(i) - 1));
        massBlock[i]->Multiply(false, *splitVecOld, *splitVec);
        splitRes->Update(-1.0, *splitVec, 1.0);
      }	
  
      // solve the problem, and put solution into slution vector
      splitVec->PutScalar(0.0);
      linSys[i]->applyJacobianInverse(lsParams, *splitRes_NEV, *splitVec_NEV);
      solution.LoadBlockValues(*splitVec, jacobian.RowIndex(i));
    }

    result = solution;
  }    

  return 0;
}

double LOCA::Epetra::xyztPrec::
NormInf() const
{
  // Not implemented: Throw an error!
  std::cout << "ERROR: LOCA::Epetra::xyztPrec::NormInf() - "
       << "method is NOT implemented!!  " << std::endl;
  throw "LOCA Error";
  return 0.0;
}

bool LOCA::Epetra::xyztPrec::
UseTranspose() const
{
  // Not implemented: Throw an error!
  std::cout << "ERROR: LOCA::Epetra::xyztPrec::UseTranspose() - "
       << "method is NOT implemented!!  " << std::endl;
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
	 const NOX::Epetra::Interface::Required::FillType)
{
  //cout << "LOCA::Epetra::xyztPrec::computeF called." << std::endl;
  return false;
}

bool LOCA::Epetra::xyztPrec::
computeJacobian(const Epetra_Vector&, Epetra_Operator&)
{
  //cout << "LOCA::Epetra::xyztPrec::computeJacobian called." << std::endl;
  return false;
}

bool LOCA::Epetra::xyztPrec::
computePreconditioner(const Epetra_Vector& x,
		      Epetra_Operator& Prec,
		      Teuchos::ParameterList* p)
{
  //cout << "LOCA::Epetra::xyztPrec::computePreconditioner called - ";
  std::string prec = lsParams.get("XYZTPreconditioner","None");
  if (prec == "Global") {
    //cout << "Global Preconditioner" << std::endl;
    linSys[0]->destroyPreconditioner();
    linSys[0]->createPreconditioner(x, lsParams, false);
    return true;
  }    
  else if (prec == "Sequential"  || prec == "Parallel" ||
           prec == "BlockDiagonal" || prec == "Parareal" || prec == "BDSDT") {
    for (int i=0; i< globalComm->NumTimeStepsOnDomain(); i++ ) {
      if (globalComm->FirstTimeStepOnDomain() + i == 0 && !isPeriodic)
        jacobian.ExtractBlock(*jacobianBlock[i], 0, 0);
      else {
        jacobian.ExtractBlock(*jacobianBlock[i], i, 1);
        if (prec != "BlockDiagonal") jacobian.ExtractBlock(*massBlock[i], i, 0);
      }
      //is following line needed???
      linSys[i]->setJacobianOperatorForSolve(jacobianBlock[i]);
      linSys[i]->destroyPreconditioner();
      linSys[i]->createPreconditioner(x, lsParams, false);
    }
    if (prec == "Parareal") {
      int N = globalComm->NumTimeStepsOnDomain();
      jacobian.ExtractBlock(*jacobianBlock[N], N-1, 1);
      jacobian.ExtractBlock(*massBlock[N], N-1, 0);

      //Transform to matrices for N steps
      (*massBlock[N]).Scale(1.0 / (double)N);
      EpetraExt::MatrixMatrix::Add(*massBlock[N], false, (double)(N-1), *jacobianBlock[N], 1.0);

      linSys[N]->setJacobianOperatorForSolve(jacobianBlock[N]);
      linSys[N]->destroyPreconditioner();
      linSys[N]->createPreconditioner(x, lsParams, false);
    }
    if (prec == "BDSDT") {
      for (int i=0; i< globalComm->NumTimeStepsOnDomain(); i++ ) {
	int j=i-1;
        if (globalComm->FirstTimeStepOnDomain() + i != 0) {
	   if (j<0) j=0; // Jac is from previous block, but don't go off proc to get it
          (*jacobianBlock[j]).ExtractDiagonalCopy(*splitVec);
          (*massBlock[i]).ExtractDiagonalCopy(*splitRes);
	  (*diagBlockSubdiag[i]).ReciprocalMultiply(1.0, *splitVec, *splitRes, 0.0);
	}
      }
    }
    return true;
  }
  else {
    //cout << "No Preconditioner" << std::endl;
    return true;
  }
}

void LOCA::Epetra::xyztPrec::
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  std::cout << "LOCA::Epetra::xyztPrec" << functionName 
	 << " - " << errorMsg << std::endl;
  throw "LOCA Error";
}

