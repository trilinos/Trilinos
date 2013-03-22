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

#include "LOCA_Epetra_Interface_MultiPoint.H"
  
#ifdef HAVE_NOX_EPETRAEXT

LOCA::Epetra::Interface::MultiPoint::
MultiPoint(
       const Teuchos::RCP<LOCA::Epetra::Interface::Required> &iReq_,
       const Teuchos::RCP< NOX::Epetra::Interface::Jacobian> &iJac_,
       const Epetra_MultiVector &splitMultiVec_, 
       const Teuchos::RCP<Epetra_RowMatrix> &splitJac_,
       const Teuchos::RCP<EpetraExt::MultiComm> &globalComm_) :
  iReq(iReq_),
  iJac(iJac_),
  splitJac(splitJac_), 
  globalComm(globalComm_),
  splitVec(*(splitMultiVec_(0))),
  splitRes(*(splitMultiVec_(0))), 
  jacobian(0), 
  solution(0),
  solutionOverlap(0), 
  overlapImporter(0), 
  timeStepsOnTimeDomain(splitMultiVec_.NumVectors()), 
  numTimeDomains(globalComm_->NumSubDomains()),
  timeDomain(globalComm_->SubDomainRank()), 
  conStep(0),
  rowStencil(0),
  rowIndex(0)
{

   if (globalComm->MyPID()==0) {
     // TODO: pass in globalData and use output stream
     std::cout  << "----------MultiPoint Partition Info------------"
           << "\n\tNumProcs              = " << globalComm->NumProc()
           << "\n\tSpatial Decomposition = " << splitMultiVec_.Comm().NumProc()
           << "\n\tNumber of Domains     = " << numTimeDomains
           << "\n\tSteps on Domain 0     = " << timeStepsOnTimeDomain
           << "\n\tTotal Number of Steps = " << globalComm->NumTimeSteps();
    std::cout   << "\n-----------------------------------------------" << std::endl;
    }

   // Construct global block matrix graph from split jacobian and stencil,
   // which is just diagonal in this case

   rowStencil = new std::vector< std::vector<int> >(timeStepsOnTimeDomain);
   rowIndex = new std::vector<int>;
   for (int i=0; i < timeStepsOnTimeDomain; i++) {
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
}

LOCA::Epetra::Interface::MultiPoint::
~MultiPoint()
{
  delete solution;
  delete solutionOverlap;
  delete overlapImporter;
  delete jacobian;
  delete rowStencil;
  delete rowIndex;
}

bool LOCA::Epetra::Interface::MultiPoint::
computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  bool stat = true;

  // Copy owned parts of vector from vector with global map 
  // to one with split map
  solution->Epetra_Vector::operator=(x);
  solutionOverlap->Import(*solution, *overlapImporter, Insert);

  EpetraExt::BlockVector residual(*solution);

  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);

    splitRes.PutScalar(0.0);

    // Pass step index so application can adjust parameters
    iReq->setMultiPointParameter((*rowIndex)[i]); 
    stat = iReq->computeF(splitVec, splitRes, fillFlag);

    residual.LoadBlockValues(splitRes, (*rowIndex)[i]);
  }

  // F does not know it is a clone of a block vector, 
  // so must copy values from residual
  // -- can be fixed? Maybe make residual a view of F instad of copying.
  F = residual;

  return stat;
}

bool LOCA::Epetra::Interface::MultiPoint::
computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  bool stat = true;
  jacobian->PutScalar(0.0);

  // Copy owned parts of vector from vector with global map 
  // to one with split map
  solution->Epetra_Vector::operator=(x);

  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);

    // Pass step index so application can adjust parameters
    iReq->setMultiPointParameter((*rowIndex)[i]); 
    stat =  iJac->computeJacobian(splitVec, *splitJac );

    jacobian->LoadBlock(*splitJac, i, 0);
  }

  return stat;
}

void LOCA::Epetra::Interface::MultiPoint::
setParameters(const LOCA::ParameterVector& params)
{
   iReq->setParameters(params);
}

void LOCA::Epetra::Interface::MultiPoint::
printSolution(const Epetra_Vector& x, double conParam)
{
  solution->Epetra_Vector::operator=(x);

  // Barriers force printing of all time steps in order
  for (int j=0; j<timeDomain; j++) globalComm->Barrier();
  for (int i=0; i < timeStepsOnTimeDomain; i++) {
    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);
    // Pass indexing data for possible application use in 
    // output naming convention
    iReq->dataForPrintSolution(conStep, (*rowIndex)[i], 
	 		       globalComm->NumTimeSteps());
    iReq->printSolution(splitVec, conParam);
  }
  for (int j=timeDomain; j<numTimeDomains-1; j++) globalComm->Barrier();

  conStep++; // Counter for continuation step, used for printing
}

EpetraExt::BlockVector& LOCA::Epetra::Interface::MultiPoint::
getSolution()
{
  return  *solution; 
}

EpetraExt::BlockCrsMatrix& LOCA::Epetra::Interface::MultiPoint::
getJacobian()
{
  return  *jacobian; 
}

void LOCA::Epetra::Interface::MultiPoint::
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  std::cout << "LOCA::Epetra::Interface::MultiPoint::" << functionName 
	 << " - " << errorMsg << std::endl;
  throw "LOCA Error";
}

#endif // HAVE_NOX_EPETRAEXT
