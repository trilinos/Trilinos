//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "LOCA_Epetra_Interface_xyzt.H"
  
#ifdef HAVE_NOX_EPETRAEXT

#include "EpetraExt_MatrixMatrix.h"

LOCA::Epetra::Interface::xyzt::
xyzt(
       const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> &iReq_,
       const Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> &iJac_,
       const Teuchos::RefCountPtr<LOCA::Epetra::Interface::MassMatrix> &iMass_,
       const Epetra_MultiVector &splitMultiVec_, 
       const Teuchos::RefCountPtr<Epetra_RowMatrix> &splitJac_,
       const Teuchos::RefCountPtr<Epetra_RowMatrix> &splitMass_, 
       const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> &globalComm_,
       Teuchos::ParameterList *precPrintParams_,
       Teuchos::ParameterList *precLSParams_) :
  iReq(iReq_), 
  iJac(iJac_), 
  iMass(iMass_),
  splitJac(splitJac_), 
  splitMass(splitMass_), 
  globalComm(globalComm_),
  splitVec(*(splitMultiVec_(0))),
  splitRes(*(splitMultiVec_(0))), 
  splitVecOld(*(splitMultiVec_(0))), 
  jacobian(0), 
  solution(0),
  solutionOverlap(0), 
  overlapImporter(0), 
  timeStepsOnTimeDomain(splitMultiVec_.NumVectors()), 
  numTimeDomains(globalComm_->NumSubDomains()),
  timeDomain(globalComm_->SubDomainRank()), 
  conStep(0),
  rowStencil(0),
  rowIndex(0), 
  precPrintParams(precPrintParams_), 
  precLSParams(precLSParams_)
{
   if (precLSParams)
     isPeriodic = precLSParams_->get("Periodic",false);
   else
     isPeriodic = false;

   if (globalComm->MyPID()==0) {
     // TODO: pass in globalData and use output stream
     cout  << "--------------XYZT Partition Info---------------"
           << "\n\tNumProcs               = " << globalComm->NumProc()
           << "\n\tSpatial Decomposition  = " << splitMultiVec_.Comm().NumProc()
           << "\n\tNumber of Time Domains = " << numTimeDomains
           << "\n\tTime Steps on Domain 0 = " << timeStepsOnTimeDomain
           << "\n\tNumber of Time Steps   = " << globalComm->NumTimeSteps();
    if (isPeriodic) cout   << "\n\t-->Solving for a Periodic Orbit!" ;
    cout   << "\n-----------------------------------------------" << endl;
    }

   // Construct global block matrix graph from split jacobian and stencil
   // Each block has identical sparsity, and assumes mass matrix's sparsity 
   // is a subset of the Jacobian's

   rowStencil = new std::vector< std::vector<int> >(timeStepsOnTimeDomain);
   rowIndex = new std::vector<int>;
   for (int i=0; i < timeStepsOnTimeDomain; i++) {
     if (timeDomain!=0 || i!=0)  
       (*rowStencil)[i].push_back(-1);
     else if (isPeriodic)   
       (*rowStencil)[i].push_back(globalComm->NumTimeSteps()-1);
     (*rowStencil)[i].push_back(0);
     (*rowIndex).push_back(i + globalComm->FirstTimeStepOnDomain());
   }

   jacobian = new EpetraExt::BlockCrsMatrix(*splitJac, *rowStencil, 
					    *rowIndex, *globalComm);

   // Construct global solution vector, the overlap vector, 
   //and importer between them
   solution = new EpetraExt::BlockVector(splitJac->RowMatrixRowMap(), 
					 jacobian->RowMap());
   solutionOverlap = new EpetraExt::BlockVector(splitJac->RowMatrixRowMap(), 
						jacobian->ColMap());
  
   overlapImporter = new Epetra_Import(solutionOverlap->Map(), solution->Map());


   // Load initial guess into block solution vector
   for (int i=0; i < timeStepsOnTimeDomain; i++) 
           solution->LoadBlockValues(*(splitMultiVec_(i)), (*rowIndex)[i]);

   // Create preconditioner
   if (precLSParams != 0) {
     //Preconditioner needs CrsMatrix, must convert VBR or others
     
     splitJacCrs = dynamic_cast<Epetra_CrsMatrix *>(splitJac.get());
     if (splitJacCrs == NULL) {
        cout << "CAST OF splitJacCrs failed!, constructing CRS matrix " << endl;

        std::vector< std::vector<int> > row(1); row[0].push_back(0);
	std::vector<int> col; col.push_back(0);
        splitJacCrs = (Epetra_CrsMatrix *)
	   new EpetraExt::BlockCrsMatrix(*splitJac, row, col, 
					 splitJac->Comm());
        splitMassCrs = (Epetra_CrsMatrix *)
	   new EpetraExt::BlockCrsMatrix(*splitMass, row, col, 
					 splitMass->Comm());
     }
     else splitMassCrs = dynamic_cast<Epetra_CrsMatrix *>(splitMass.get());

     preconditioner = 
       new LOCA::Epetra::xyztPrec(*jacobian, *splitJacCrs, *splitMassCrs, 
				  *solution, *solutionOverlap, *overlapImporter,
				  *precPrintParams, *precLSParams, globalComm);
     if (preconditioner != 0) {
       // TODO: pass in globalData and use output stream
       //cout << "LOCA::Epetra::Interface::xyzt - " 
       //    << "preconditioner created successfully" << endl;
     }
   }
   // TODO: pass in globalData and use output stream
   //cout << "Ending xyzt constructor" << endl;
}

LOCA::Epetra::Interface::xyzt::
~xyzt()
{
  delete solution;
  delete solutionOverlap;
  delete overlapImporter;
  delete jacobian;
  delete rowStencil;
  delete rowIndex;
}

bool LOCA::Epetra::Interface::xyzt::
computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  bool stat = true;

  // Copy owned parts of vector from vector with global map 
  // to one with split map
  solution->Epetra_Vector::operator=(x);
  solutionOverlap->Import(*solution, *overlapImporter, Insert);

  EpetraExt::BlockVector residual(*solution);

  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    if (i==0 && timeDomain==0 && !isPeriodic) {
      iMass->setOldSolutionFirstStep();
    }
    else {
      int blockRowOld = (*rowIndex)[i] + (*rowStencil)[i][0];
      solutionOverlap->ExtractBlockValues(splitVecOld, blockRowOld);
      iMass->setOldSolution(splitVecOld, 
			    i + globalComm->FirstTimeStepOnDomain());
    }

    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);

    splitRes.PutScalar(0.0);
    stat = stat && iReq->computeF(splitVec,  splitRes, fillFlag);

    residual.LoadBlockValues(splitRes, (*rowIndex)[i]);
  }

  // F does not know it is a clone of a block vector, 
  // so must copy values from residual
  // -- can be fixed? Maybe make residual a view of F instad of copying.
  F = residual;

  return stat;
}

bool LOCA::Epetra::Interface::xyzt::
computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  bool stat = true;
  jacobian->PutScalar(0.0);

  // Copy owned parts of vector from vector with global map 
  // to one with split map
  solution->Epetra_Vector::operator=(x);
  solutionOverlap->Import(*solution, *overlapImporter, Insert);

  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    if (i==0 && timeDomain==0 && !isPeriodic) {
      //First time step gets static Xold, not part of block solution vector
      iMass->setOldSolutionFirstStep();
    }
    else {
      int blockRowOld = (*rowIndex)[i] + (*rowStencil)[i][0];
      solutionOverlap->ExtractBlockValues(splitVecOld, blockRowOld);
      iMass->setOldSolution(splitVecOld, 
			    i + globalComm->FirstTimeStepOnDomain());
    }
  
    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);
    stat =  stat && iJac->computeJacobian( splitVec, Jac );

    // Hardwired for -1 0 stencil meaning [M J]
    if (i==0 && timeDomain==0 && !isPeriodic) {
      jacobian->LoadBlock(*splitJac, i, 0);
    }
    else {
      jacobian->LoadBlock(*splitJac, i, 1);
      stat = stat && iMass->computeMassMatrix( splitVecOld );
      jacobian->LoadBlock(*splitMass, i, 0);
    }
  }

  return stat;
}

void LOCA::Epetra::Interface::xyzt::
setParameters(const LOCA::ParameterVector& params)
{
   iReq->setParameters(params);
}

void LOCA::Epetra::Interface::xyzt::
printSolution(const Epetra_Vector& x, double conParam)
{
  solution->Epetra_Vector::operator=(x);

  // Barriers force printing of all time steps in order
  for (int j=0; j<timeDomain; j++) globalComm->Barrier();
  for (int i=0; i < timeStepsOnTimeDomain; i++) {
    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);
    // Pass indexing data for possible application use in 
    // output naming convention
    iMass->dataForPrintSolution(conStep, (*rowIndex)[i], 
				globalComm->NumTimeSteps());
    iReq->printSolution(splitVec, conParam);
  }
  for (int j=timeDomain; j<numTimeDomains-1; j++) globalComm->Barrier();

  conStep++; // Counter for continuation step, used for printing
}

EpetraExt::BlockVector& LOCA::Epetra::Interface::xyzt::
getSolution()
{
  return  *solution; 
}

EpetraExt::BlockCrsMatrix& LOCA::Epetra::Interface::xyzt::
getJacobian()
{
  return  *jacobian; 
}

LOCA::Epetra::xyztPrec& LOCA::Epetra::Interface::xyzt::
getPreconditioner()
{
  return  *preconditioner; 
}

void LOCA::Epetra::Interface::xyzt::
throwError(const string& functionName, const string& errorMsg) const
{
  cout << "LOCA::Epetra::Interface::xyzt::" << functionName 
	 << " - " << errorMsg << endl;
  throw "LOCA Error";
}

#endif // HAVE_NOX_EPETRAEXT
