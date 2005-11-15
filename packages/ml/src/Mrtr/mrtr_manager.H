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

#ifndef MOERTEL_MANAGER_H
#define MOERTEL_MANAGER_H

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
#include "Teuchos_RefCountPtr.hpp"
// mrtr includes
#include "mrtr_segment.H"
#include "mrtr_node.H"
#include "mrtr_interface.H"

using namespace std;
using namespace Teuchos;

namespace MOERTEL
{

class Manager 
{
public:

  enum DimensionType
  {
    manager_none,
    manager_2D,
    manager_3D
  };
  
  // ctor
  explicit Manager(Epetra_Comm& comm, int outlevel);
  
  // dtor
  virtual ~Manager();
  
  // get Epetra_Comm
  Epetra_Comm& Comm() const {return comm_;}

  // get Printlevel (0-10)
  int OutLevel() { return outlevel_; }

  // set the dimension type
  void SetDimension(MOERTEL::Manager::DimensionType type) { dimensiontype_ = type; return; }
  
  // get the dimension type
  MOERTEL::Manager::DimensionType Dimension() { return dimensiontype_; }
  
  // add an interface to the mortar manager
  bool AddInterface(MOERTEL::Interface& interface);

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
  Manager(MOERTEL::Manager& old);

  // do integration of all interfaces 2D case
  bool Mortar_Integrate_2D();
  
  // do integration of all interfaces 3D case
  bool Mortar_Integrate_3D();
  
  // build the map for the saddle point problem
  bool BuildSaddleMap();

  // choose the dofs for lagrange multipliers
  RefCountPtr<Epetra_Map> LagrangeMultiplierDofs();
  
  // modified version of the epetraext matrixmatrixadd
  int MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB );

  // automatically choose mortar side (called when mortarside==-2 on any interface)
  bool ChooseMortarSide();
  bool ChooseMortarSide(vector<RefCountPtr<MOERTEL::Interface> >& inter);
  
private:

  int                           outlevel_;            // output level (0-10)
  Epetra_Comm&                  comm_;                // communicator (global, contains ALL procs)
  DimensionType                 dimensiontype_;       // problem dimension

  RefCountPtr<Epetra_Map>       inputmap_;            // the rowmap of the input matrix
  RefCountPtr<Epetra_CrsMatrix> inputmatrix_;         // the uncoupled matrix on input
  RefCountPtr<Epetra_Map>       constraintsmap_;      // the rowmap of M and D (both of them use comm_)
  RefCountPtr<Epetra_CrsMatrix> D_;                   // the coupling matrix D
  RefCountPtr<Epetra_CrsMatrix> M_;                   // the coupling matrix M
  RefCountPtr<Epetra_Map>       saddlemap_;           // the rowmap of the saddlepointproblem
  RefCountPtr<Epetra_CrsMatrix> saddlematrix_;        // the matrix of the saddle problem;
  
  map<int,RefCountPtr<MOERTEL::Interface> >  interface_; // the interfaces
};

} // namespace MOERTEL

// << operator
std::ostream& operator << (ostream& os, const MOERTEL::Manager& node);

#endif // MOERTEL_MANAGER_H
#endif // TRILINOS_PACKAGE
