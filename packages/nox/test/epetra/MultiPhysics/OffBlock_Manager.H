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

class OffBlock_Manager : GenericEpetraProblem 
{

  friend class Problem_Manager;

public:

  enum VectorType { SOLUTION, GROUP_F };

  static int idToFind;

  // Constructor
  OffBlock_Manager(Problem_Manager& problemMan_, Epetra_CrsGraph& graph_,
                   int probEqId, int probVarsId);

  // Destructor
  virtual ~OffBlock_Manager()
  { };

  // Return the underlying Group 
  Teuchos::RefCountPtr<NOX::Epetra::Group> getGroup();

  // Return the underlying CrsMatrix 
  Epetra_CrsMatrix& getMatrix();

  // Return the rowMapVec
  Teuchos::RefCountPtr<Epetra_Vector>  getRowMapVec() const;

  // Return the problemEqId
  int getProblemEqId() const;

  // Return the problemVarId
  int getProblemVarId() const;

  // Convert block indices to composite indices
  void convertBlockRowIndicesToComposite( int numIndices, int * block, int * composite);
  void convertBlockColIndicesToComposite( int numIndices, int * block, int * composite);


  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(NOX::Epetra::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector);

  // A helper needed to set a target for a predicate
  void setIDtoFind( int id )
  { idToFind = id; }

  // And here is the predicate
  bool isMember()
  { return problemVarId == idToFind; }

private:

  // Creates the extra objects needed to fill off-diagoanl blocks using FD
  void createFDobjects( bool useColoring = false );

  // Creates a graph with indices shifted as needed for coloring
  Epetra_CrsGraph & createBlockGraphFromComposite(Epetra_CrsGraph &);

private:

  int problemEqId;      // Problem Id of off-diagonal block equations
  int problemVarId;     // Problem Id of off-diagonal block variables

  Teuchos::RefCountPtr<Problem_Interface> offBlockInterface; // A problem interface to 
  // this off-block problem manager.  

  // These allow off-diagonal block contributions representing
  // inter-problem coupling computed via FDC
  Teuchos::RefCountPtr<Epetra_Map>                              blockRowMap     ;
  Teuchos::RefCountPtr<Epetra_Map>                              blockColMap     ;
  Teuchos::RefCountPtr<Epetra_Vector>                           rowMapVec       ;
  Teuchos::RefCountPtr<Epetra_CrsGraph>                         graph           ;
  Teuchos::RefCountPtr<NOX::Epetra::Vector>                     noxVec          ;
  Teuchos::RefCountPtr<EpetraExt::CrsGraph_MapColoring>         mapColoring     ;
  Teuchos::RefCountPtr<Epetra_MapColoring>                      colorMap        ;
  Teuchos::RefCountPtr<EpetraExt::CrsGraph_MapColoringIndex>    colorMapIndex   ;
  Teuchos::RefCountPtr<vector<Epetra_IntVector> >               columnSet       ;
  Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference>           matrixOperator  ;
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO>        linearSystem    ;
  Teuchos::RefCountPtr<NOX::Epetra::Group>                      group           ;
  
  // These are needed to shift indices for coloring and 
  // block-based condition number estimates
  map< int, int > rowBlockToComposite;
  map< int, int > rowCompositeToBlock;
  map< int, int > colBlockToComposite;
  map< int, int > colCompositeToBlock;

};

#endif

#endif
