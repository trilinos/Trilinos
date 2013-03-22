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

#include "LOCA_Epetra_Interface_xyzt.H"
  
#ifdef HAVE_NOX_EPETRAEXT

#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_VbrMatrix.h"

LOCA::Epetra::Interface::xyzt::
xyzt( const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> &interface_,
       const Epetra_MultiVector &splitMultiVec_, 
       const Teuchos::RCP<Epetra_RowMatrix> &splitJac_,
       const Teuchos::RCP<EpetraExt::MultiComm> &globalComm_,
       const Epetra_Vector &initialCondVec_, 
       double dt_,
       Teuchos::ParameterList *precPrintParams_,
       Teuchos::ParameterList *precLSParams_) :
  interface(interface_),
  splitJac(splitJac_), 
  globalComm(globalComm_),
  splitVec(*(splitMultiVec_(0))),
  splitRes(*(splitMultiVec_(0))), 
  splitVecOld(*(splitMultiVec_(0))), 
  initialCondVec(initialCondVec_),
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
  precLSParams(precLSParams_),
  savedSplitMassForFloquet(0),
  isCrsMatrix(true),
  floquetFillFlag(false),
  dt(dt_)
{
   if (precLSParams)
     isPeriodic = precLSParams_->get("Periodic",false);
   else
     isPeriodic = false;

   if (globalComm->MyPID()==0) {
     // TODO: pass in globalData and use output stream
     std::cout  << "--------------XYZT Partition Info---------------"
           << "\n\tNumProcs               = " << globalComm->NumProc()
           << "\n\tSpatial Decomposition  = " << splitMultiVec_.Comm().NumProc()
           << "\n\tNumber of Time Domains = " << numTimeDomains
           << "\n\tTime Steps on Domain 0 = " << timeStepsOnTimeDomain
           << "\n\tNumber of Time Steps   = " << globalComm->NumTimeSteps();
    if (isPeriodic) std::cout   << "\n\t-->Solving for a Periodic Orbit!" ;
    std::cout   << "\n-----------------------------------------------" << std::endl;
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
        isCrsMatrix = false;
        std::cout << "CAST OF splitJacCrs failed!, constructing CRS matrix " << std::endl;

        std::vector< std::vector<int> > row(1); row[0].push_back(0);
	std::vector<int> col; col.push_back(0);
        splitJacCrs = (Epetra_CrsMatrix *)
	   new EpetraExt::BlockCrsMatrix(*splitJac, row, col, 
					 splitJac->Comm());
     }

     preconditioner = 
       Teuchos::rcp(new LOCA::Epetra::xyztPrec(*jacobian, *splitJacCrs, *solution,
                                  *solutionOverlap, *overlapImporter,
                                  *precPrintParams, *precLSParams, globalComm));
   }
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

    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);

    /* get solution vector at previous step */
    if (i==0 && timeDomain==0 && !isPeriodic) {
      splitVecOld = initialCondVec;
    }
    else {
      int blockRowOld = (*rowIndex)[i] + (*rowStencil)[i][0];
      solutionOverlap->ExtractBlockValues(splitVecOld, blockRowOld);
    }

    /* calc xDot using two vectors -- generalize later */
    Epetra_Vector& vecDot = splitRes; //use Res as temp space for xdot
    vecDot.Update(1.0/dt, splitVec, -1.0/dt, splitVecOld, 0.0);
    interface->setXdot(vecDot, dt*(i+globalComm->FirstTimeStepOnDomain()));

    splitRes.PutScalar(0.0);
    stat = stat && interface->computeF(splitVec,  splitRes, fillFlag);

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

    solution->ExtractBlockValues(splitVec, (*rowIndex)[i]);

    /* get solution vector at previous step */
    if (i==0 && timeDomain==0 && !isPeriodic) {
      splitVecOld = initialCondVec;
    }
    else {
      int blockRowOld = (*rowIndex)[i] + (*rowStencil)[i][0];
      solutionOverlap->ExtractBlockValues(splitVecOld, blockRowOld);
    }

    /* calc xDot using two vectors -- generalize later */
    Epetra_Vector& vecDot = splitRes; //use Res as temp space for xdot
    vecDot.Update(1.0/dt, splitVec, -1.0/dt, splitVecOld, 0.0);
    interface->setXdot(vecDot, dt*(i+globalComm->FirstTimeStepOnDomain()));
  
    stat =  stat && interface->computeShiftedMatrix(1.0, 1.0/dt, splitVec, *splitJac );

    // Hardwired for -1 0 stencil meaning [M J]
    if (i==0 && timeDomain==0 && !isPeriodic) {
      jacobian->LoadBlock(*splitJac, i, 0);
    }
    else {
      jacobian->LoadBlock(*splitJac, i, 1);
      stat = stat && interface->computeShiftedMatrix(0.0, -1.0/dt,  splitVecOld, *splitJac );
      // Floquet fills, save first mass matrix instead of loading it
      if (i==0 && timeDomain==0 && isPeriodic && floquetFillFlag)
              (*savedSplitMassForFloquet) = *(splitJac.get()); 
      else jacobian->LoadBlock(*splitJac, i, 0); //Normal case
    }
  }

  return stat;
}

void LOCA::Epetra::Interface::xyzt::
setParameters(const LOCA::ParameterVector& params)
{
   interface->setParameters(params);
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
    interface->dataForPrintSolution(conStep, (*rowIndex)[i], 
				globalComm->NumTimeSteps());
    interface->printSolution(splitVec, conParam);
  }
  for (int j=timeDomain; j<numTimeDomains-1; j++) globalComm->Barrier();

  conStep++; // Counter for continuation step, used for printing
}

void LOCA::Epetra::Interface::xyzt::
setFloquetFillFlag(bool fff)
{
   if (fff) {
     floquetFillFlag = true;
     // allocate saved Mass matrix
     if (isCrsMatrix)
       savedSplitMassForFloquet = new Epetra_CrsMatrix(dynamic_cast<Epetra_CrsMatrix&>(*(splitJac.get())));
     else
       savedSplitMassForFloquet = new Epetra_VbrMatrix(dynamic_cast<Epetra_VbrMatrix&>(*(splitJac.get())));
   }
   else {
     floquetFillFlag = false;
     // free the saved Mass matrix
     delete savedSplitMassForFloquet;
   }
}

void LOCA::Epetra::Interface::xyzt::
beginFloquetOperatorApplication(Epetra_Vector& v)
{
  if (!isPeriodic)
    std::cout << "\n\n\t Must be periodic for Floquet theory\n" << std::endl;
  // Take perturbation in final time step, apply the saved
  // mass matrix, and put as forcing of first time step

  // Create BlockVectors for manipulation of single blocks
  solution->Epetra_Vector::operator=(v);
  solutionOverlap->Import(*solution, *overlapImporter, Insert);

  // Set vector to all zeros -- 
  solution->PutScalar(0.0);

  // only process perturbation to last vector, and place in first vector
  if (timeDomain == 0) {
    int blockRowOld = (*rowIndex)[0] + (*rowStencil)[0][0];
    solutionOverlap->ExtractBlockValues(splitVecOld, blockRowOld);
    savedSplitMassForFloquet->Multiply(false, splitVecOld, splitVec);
    splitVec.Scale(-1.0); // subtract over to RHS
    solution->LoadBlockValues(splitVec, 0);
  }

  v = *solution;
}

void LOCA::Epetra::Interface::xyzt::
finishFloquetOperatorApplication(Epetra_Vector& v)
{
  // zero out all component vectors except the final one...
  solution->Epetra_Vector::operator=(v);

  if ( timeDomain == numTimeDomains-1 )
    solution->ExtractBlockValues(splitVec, globalComm->NumTimeSteps()-1);

  solution->PutScalar(0.0);

  if ( timeDomain == numTimeDomains-1 )
    solution->LoadBlockValues(splitVec, globalComm->NumTimeSteps()-1);

   v = *solution;
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
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  std::cout << "LOCA::Epetra::Interface::xyzt::" << functionName 
	 << " - " << errorMsg << std::endl;
  throw "LOCA Error";
}

#endif // HAVE_NOX_EPETRAEXT
