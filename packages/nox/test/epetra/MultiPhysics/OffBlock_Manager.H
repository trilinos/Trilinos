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
                                                                                
#ifndef _NOX_OFFBLOCK_MANAGER_H
#define _NOX_OFFBLOCK_MANAGER_H


#include "NOX.H"
// Currently, this class can only be used with builds which include EpetraExt !!
#ifdef HAVE_NOX_EPETRAEXT

#include "NOX_Epetra.H"

#include "Problem_Interface.H"
#include "Problem_Manager.H"
#include "GenericEpetraProblem.H"

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif  

#include <vector>

class OffBlock_Manager : GenericEpetraProblem { 

  friend class Problem_Manager;

public:

  enum VectorType { SOLUTION, GROUP_F };

  // Constructor
  OffBlock_Manager(Problem_Manager& problemMan_, Epetra_CrsGraph& graph_,
                   int probEqId, int probVarsId);

  // Destructor
  virtual ~OffBlock_Manager();

  // Return the underlying Group 
  NOX::EpetraNew::Group& getGroup();

  // Return the underlying CrsMatrix 
  Epetra_CrsMatrix& getMatrix();


  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(NOX::EpetraNew::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix);

private:

  // Creates the extra objects needed to fill off-diagoanl blocks using FDC
  void createFDCobjects();

private:

  int problemEqId;      // Number of registered Problems
  int problemVarId;     // Number of registered Problems

  Problem_Manager &myManager;  //Reference to the overall problem manager

  Problem_Interface *offBlockInterface; // A problem interface to 
  // this off-block problem manager.  

  // These allow off-diagonal block contributions representing
  // inter-problem coupling computed via FDC
  EpetraExt::CrsGraph_MapColoring* mapColoring;
  Epetra_MapColoring* colorMap;
  EpetraExt::CrsGraph_MapColoringIndex* colorMapIndex;
  vector<Epetra_IntVector>* columnSet;
  NOX::EpetraNew::FiniteDifferenceColoring *matrixOperator;
  NOX::EpetraNew::LinearSystemAztecOO* linearSystem;
  NOX::EpetraNew::Group* group;    
};

#endif

#endif
