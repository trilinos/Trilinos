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
  
using namespace  LOCA::EpetraNew::Interface;

xyzt::xyzt(NOX::EpetraNew::Interface::Required &iReq_, NOX::EpetraNew::Interface::Jacobian &iJac_,
       Epetra_Vector &splitVec_, Epetra_RowMatrix &splitJac_, Epetra_MpiComm &globalComm_,
       int replica_) : iReq(iReq_), iJac(iJac_), splitVec(splitVec_), splitJac(splitJac_),
       globalComm(globalComm_), replica(replica_), solution(0), jacobian(0)
{
  cout << " In xyzt constructor: global rank " << globalComm.MyPID() 
       << " and total procs " << globalComm.NumProc()<< endl;
  cout << "                      local  rank " << splitVec.Comm().MyPID() 
       << " and local procs " << splitVec.Comm().NumProc()<< endl;
    
   // Bogus redundant calcs for now
//   solution = new Epetra_Vector(splitVec);
//   jacobian = new Epetra_RowMatrix(splitJac);
}

xyzt::~xyzt()
{
//  delete solution;
//  delete jacobian;
}

bool xyzt::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  return iReq.computeF(x, F, fillFlag);
}

bool xyzt::computeJacobian(const Epetra_Vector& x)
{
  return iJac.computeJacobian(x);
}

Epetra_Vector& xyzt::getSolution()
{
    return  splitVec;
//  return  solution; 
}

Epetra_RowMatrix& xyzt::getJacobian()
{
  return  splitJac; 
  //return  jacobian; 
}
