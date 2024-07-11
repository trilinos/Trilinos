// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  Teuchos::RCP<NOX::Epetra::Group> getGroup();

  // Return the underlying CrsMatrix
  Epetra_CrsMatrix& getMatrix();

  // Return the rowMapVec
  Teuchos::RCP<Epetra_Vector>  getRowMapVec() const;

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

  Teuchos::RCP<Problem_Interface> offBlockInterface; // A problem interface to
  // this off-block problem manager.

  // These allow off-diagonal block contributions representing
  // inter-problem coupling computed via FDC
  Teuchos::RCP<Epetra_Map>                              blockRowMap     ;
  Teuchos::RCP<Epetra_Map>                              blockColMap     ;
  Teuchos::RCP<Epetra_Vector>                           rowMapVec       ;
  Teuchos::RCP<Epetra_CrsGraph>                         graph           ;
  Teuchos::RCP<NOX::Epetra::Vector>                     noxVec          ;
  Teuchos::RCP<EpetraExt::CrsGraph_MapColoring>         mapColoring     ;
  Teuchos::RCP<Epetra_MapColoring>                      colorMap        ;
  Teuchos::RCP<EpetraExt::CrsGraph_MapColoringIndex>    colorMapIndex   ;
  Teuchos::RCP<vector<Epetra_IntVector> >               columnSet       ;
  Teuchos::RCP<NOX::Epetra::FiniteDifference>           matrixOperator  ;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO>        linearSystem    ;
  Teuchos::RCP<NOX::Epetra::Group>                      group           ;

  // These are needed to shift indices for coloring and
  // block-based condition number estimates
  map< int, int > rowBlockToComposite;
  map< int, int > rowCompositeToBlock;
  map< int, int > colBlockToComposite;
  map< int, int > colCompositeToBlock;

};

#endif

#endif
