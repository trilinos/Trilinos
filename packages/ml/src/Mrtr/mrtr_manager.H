/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#ifndef MRTR_MANAGER_H
#define MRTR_MANAGER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
// mrtr includes
#include "mrtr_segment.H"
#include "mrtr_node.H"
#include "mrtr_interface.H"

using namespace std;

namespace MRTR
{

class Manager 
{
public:
  
  // ctor
  explicit Manager(Epetra_Comm& comm, int outlevel);
  
  // dtor
  virtual ~Manager();
  
  // get Epetra_Comm
  Epetra_Comm& Comm() const {return comm_;}

  // get Printlevel (0-10)
  int OutLevel() { return outlevel_; }

  // add an interface to the mortar manager
  bool AddInterface(MRTR::Interface& interface);

  // printout this mortar manager
  bool Print() const;
  
  // outputlevel
  int OutLevel() const { return outlevel_; }

  // set the rowmap of the input matrix of the application
  bool SetInputMap(Epetra_Map* map);
  
  // set the uncoupled inputmatrix
  bool SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy = false);
  
  // do integration of all interfaces
  bool Mortar_Integrate();
  
  // make the saddle point problem
  Epetra_CrsMatrix* MakeSaddleProblem();
  
  // get number of interfaces
  int Ninterfaces() { return interface_.size(); }

private:  
  // don't want = operator and copy-ctor
  Manager operator = (const Manager& old);
  Manager(MRTR::Manager& old);

  // choose the dofs for lagrange multipliers
  Epetra_Map* LagrangeMultiplierDofs();
  
  // modified version of the epetraext matrixmatrixadd
  int MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB );

private:

  int                        outlevel_;       // output level (0-10)
  Epetra_Comm&               comm_;           // communicator (global, contains ALL procs)
  Epetra_Map*                inputmap_;       // the rowmap of the input matrix
  bool                       inputmatrixisdeep_; // flagindicating whether the input matrix is owned by this class
  Epetra_CrsMatrix*          inputmatrix_;    // the uncoupled matrix on input
  Epetra_Map*                constraintsmap_; // the rowmap of M and D (both of them use comm_)
  Epetra_CrsMatrix*          D_;              // the coupling matrix D
  Epetra_CrsMatrix*          M_;              // the coupling matrix M
  Epetra_Map*                saddlemap_;      // the rowmap of the saddlepointproblem
  Epetra_CrsMatrix*          saddlematrix_;   // the matrix of the saddle problem;
  
  map<int,MRTR::Interface*>  interface_;      // the interfaces
};

} // namespace MRTR

// << operator
std::ostream& operator << (ostream& os, const MRTR::Manager& node);

#endif // MRTR_MANAGER_H
#endif // TRILINOS_PACKAGE
