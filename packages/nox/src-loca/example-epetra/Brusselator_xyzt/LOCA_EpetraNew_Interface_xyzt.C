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
       Epetra_Vector &splitVec_, Epetra_CrsMatrix &splitJac_, Epetra_MpiComm &globalComm_,
       int replica_) : iReq(iReq_), iJac(iJac_), splitVec(splitVec_), splitJac(splitJac_),
       globalComm(globalComm_), replica(replica_),
       splitRes(splitVec_), solution(0), jacobian(0), stencil(0), numReplicas(-1)
{
  numReplicas = globalComm.NumProc() / splitVec.Comm().NumProc();
    cout << " In xyzt constructor: NumReplicas =  " << numReplicas << endl;

  // Allocate and fill stencil for time integratror

   std::vector<int> sten(1);
   sten[0]=0;

   cout << "xyztCalling  EpetraExt::BlockCrsMatrix" << endl;

   // Construct global block matrix graph from split matrix and stencil
   jacobian = new EpetraExt::BlockCrsMatrix(splitJac.Graph(), sten, replica, globalComm);

   cout << "xyztCalling  EpetraExt::BlockVector" << endl;
   // Construct global solution vector, and fill with initial guess
   solution = new EpetraExt::BlockVector(splitVec.Map(), jacobian->RowMap());

   solution->Block() = splitVec;

   cout << "Ending xyzt constructor" << endl;
}

xyzt::~xyzt()
{
  delete solution;
  delete jacobian;
  delete stencil;
}

bool xyzt::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{

  // Copy owned parts of vector from vector with global map to one with split map
  solution->Block() = x;
  splitVec = F;

  bool stat = iReq.computeF(solution->Block(),  splitVec, fillFlag);

  // Copy residual vector
  F = splitVec;

  return stat;
}

bool xyzt::computeJacobian(const Epetra_Vector& x)
{
  /*
  const Epetra_Vector& splitX =  (dynamic_cast<const EpetraExt::BlockVector&>(x)).Block();

  // This will load Jacobian for this time step into splitJac
  bool stat =  iJac.computeJacobian(splitX);
  */

  solution->Block() = x;
  bool stat =  iJac.computeJacobian( solution->Block() );

  EpetraExt::MatrixMatrix::Add(splitJac, false, 1.0, jacobian->Block(0), 0.0);

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

#endif
#endif
