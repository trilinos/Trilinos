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

#include "LOCA_EpetraNew_Interface_xyzt.H"
#include "EpetraExt_MatrixMatrix.h"
  
#ifdef HAVE_MPI
#ifdef HAVE_NOX_EPETRAEXT

using namespace  LOCA::EpetraNew::Interface;

xyzt::xyzt(NOX::EpetraNew::Interface::Required &iReq_, NOX::EpetraNew::Interface::Jacobian &iJac_,
       LOCA::EpetraNew::Interface::MassMatrix &iMass_, Epetra_Vector &splitVec_,
       Epetra_CrsMatrix &splitJac_, Epetra_CrsMatrix &splitMass_, Epetra_MpiComm &globalComm_,
       int replica_) : iReq(iReq_), iJac(iJac_), iMass(iMass_), splitVec(splitVec_),
       splitJac(splitJac_), splitMass(splitMass_), globalComm(globalComm_), replica(replica_),
       splitRes(splitVec_), solution(0), solutionOverlap(0), overlapImporter(0),
       jacobian(0), stencil(0), numReplicas(-1)
{
  numReplicas = globalComm.NumProc() / splitVec.Comm().NumProc();
    cout << " In xyzt constructor: NumReplicas =  " << numReplicas << endl;

  // Allocate and fill stencil for time integratror: no subdiagonal for first block row

   if (replica==0) {
     stencil = new std::vector<int>(1);
     (*stencil)[0]=0;
   } else {
     stencil = new std::vector<int>(2);
     (*stencil)[0]=-1;
     (*stencil)[1]=0;
   }

   cout << "xyztCalling  EpetraExt::BlockCrsMatrix" << endl;

   // Construct global block matrix graph from split jacobian and stencil
   // Each block has identical sparsity, and assumes mass matrix's sparsity 
   // is a subset of the Jacobian's
   jacobian = new EpetraExt::BlockCrsMatrix(splitJac.Graph(), *stencil, replica, globalComm);

   cout << "xyztCalling  EpetraExt::BlockVector" << endl;
   // Construct global solution vector, and fill with initial guess
   solution = new EpetraExt::BlockVector(splitVec.Map(), jacobian->RowMap());

   // ??
   solutionOverlap = new EpetraExt::BlockVector(splitVec.Map(), jacobian->ColMap(), stencil->size());
   overlapImporter = new Epetra_Import(solutionOverlap->Map(), solution->Map());


   solution->Block() = splitVec;

   cout << "Ending xyzt constructor" << endl;
}

xyzt::~xyzt()
{
  delete solution;
  delete solutionOverlap;
  delete overlapImporter;
  delete jacobian;
  delete stencil;
}

bool xyzt::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{

  // Copy owned parts of vector from vector with global map to one with split map
  solution->Block() = x;
  splitVec = F;

  // Set old soltuuon based  on solution from previous block row
  setOldSolution();


  bool stat = iReq.computeF(solution->Block(),  splitVec, fillFlag);

  // Copy residual vector
  F = splitVec;

  return stat;
}

bool xyzt::computeJacobian(const Epetra_Vector& x)
{

  solution->Block() = x;
  setOldSolution();
  
  bool stat =  iJac.computeJacobian( solution->Block() );

  // Put jacobian on diagonal
  if (replica==0) EpetraExt::MatrixMatrix::Add(splitJac, false, 1.0, jacobian->Block(0), 0.0);
  else {
    // Put jacobian on diagonal
    EpetraExt::MatrixMatrix::Add(splitJac, false, 1.0, jacobian->Block(0), 0.0);

    cout << " NEED TO PASS X_OLD TO MASS MATRIX FILL " << endl;
    stat = iMass.computeMassMatrix( solutionOverlap->Block(1) );
    // Put mass matrix on sub diagonal
    EpetraExt::MatrixMatrix::Add(splitJac, false, 1.0, jacobian->Block(1), 0.0);


    /*
 int NumMyRows = jacobian->Graph().NumMyRows();
 vector<double*> Values( NumMyRows );
 vector<int*> Indices( NumMyRows );
 vector<int> NumValues( NumMyRows );

 int j=0;
  jacobian->ExtractMyRowView(j, NumValues[j], Values[j], Indices[j] );

  cout << " Jacobian " << j << endl;
  for (int i=0; i<NumValues[j]; i++)
	  cout << "XX  index " << Indices[j][i] << "GID " << jacobian->GCID(Indices[j][i]) << " val " << Values[j][i] << endl;
	  */
  }


  /*
cout << "splitJac  \n" << splitJac << endl;
for (int i=0; i<1.0e7; i++) {replica++; replica--;} //pause
cout << "jacobian \n" << *jacobian << endl;
for (int i=0; i<1.0e7; i++) {replica++; replica--;} //pause
cout << "jacobian0: " << replica << "\n" <<jacobian->Block(0) << endl;
for (int i=0; i<1.0e7; i++) {replica++; replica--;} //pause
*/

  return stat;
}

EpetraExt::BlockVector& xyzt::getSolution()
{
  return  *solution; 
}

EpetraExt::BlockCrsMatrix& xyzt::getJacobian()
{
  return  *jacobian; 
}

void xyzt::setOldSolution()
{

  // This needs to be run on all processors
  solutionOverlap->Import(*solution, *overlapImporter, Insert);

  if (replica > 0) {

    cout << "Old SOlution block number -- needs fixing " << endl;
    iMass.setOldSolution(solutionOverlap->Block(1));
  }
 
}

#endif
#endif
