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
                                                                                
#include "NOX.H"
#include "NOX_Epetra.H"
// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

#include "OffBlock_Manager.H"
#include "Problem_Manager.H"
#include "GenericEpetraProblem.H"
#include "Xfer_Operator.H"

// Headers needed for Coloring
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// Header for Timing info
#include "Epetra_Time.h"

OffBlock_Manager::OffBlock_Manager(Problem_Manager& problemMan_, 
		Epetra_CrsGraph& graph_, int probEqId, int probVarId) :
  GenericEpetraProblem(graph_.Comm(), 0),
  myManager(problemMan_),
  problemEqId(probEqId),
  problemVarId(probVarId),
  mapColoring(0),
  colorMap(0),
  colorMapIndex(0),
  columnSet(0),
  matrixOperator(0),
  linearSystem(0),
  group(0)
{
 
  // Set our graph (held in base class)
  AA = &graph_;

  // Create a problem interface to the manager
  offBlockInterface = new Problem_Interface(*this);

  // Use our graph to create FDC objects needed for off-diagonal block fills
  createFDCobjects();

  // Reset number of dofs (overwrites base constructor assignment)
  NumMyNodes = 0;
}

OffBlock_Manager::~OffBlock_Manager()
{
  delete AA; AA = 0;
  delete A; A = 0;

  // Iterate over each problem and destroy/free the necessary objects

  /*
  map<int, GenericEpetraProblem*>::iterator iter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator last = Problems.end();

  map<int, NOX::EpetraNew::Group*>::iterator GroupsIter = Groups.begin();   
  map<int, OffBlock_Interface*>::iterator InterfacesIter = Interfaces.begin();
  map<int, NOX::Solver::Manager*>::iterator SolversIter = Solvers.begin();

#ifdef HAVE_NOX_EPETRAEXT
#ifdef USE_FD
  map<int, EpetraExt::CrsGraph_MapColoring*>::iterator 
	  TmpMapColoringsIter = TmpMapColorings.begin();
  map<int, Epetra_MapColoring*>::iterator 
	  ColorMapsIter = ColorMaps.begin();
  map<int, EpetraExt::CrsGraph_MapColoringIndex*>::iterator 
	  ColorMapIndexSetsIter = ColorMapIndexSets.begin();
  map<int, vector<Epetra_IntVector>*>::iterator 
	  ColumnsSetsIter = ColumnsSets.begin();
  map<int, Epetra_Operator*>::iterator 
	  MatrixOperatorsIter = MatrixOperators.begin();
#endif
#endif

  while( iter != last)
  {
    delete (SolversIter++)->second;
    delete (GroupsIter++)->second;
#ifdef HAVE_NOX_EPETRAEXT
#ifdef USE_FD
    delete (MatrixOperatorsIter++)->second;
    delete (TmpMapColoringsIter++)->second;
    delete (ColorMapsIter++)->second;
    delete (ColorMapIndexSetsIter++)->second;
    //delete *ColumnsSetsIter++;
#endif
#endif
    iter++; // Problems are owned by the app driver (Example.C)
  }
  */
}

// These methods are needed to allow inheritance from GenericEpetraProblem base

bool OffBlock_Manager::evaluate(
              NOX::EpetraNew::Interface::Required::FillType flag,
              const Epetra_Vector *solnVector,
              Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix)
{
  // Determine if fill call is valid
  if (rhsVector == 0 || flag != NOX::EpetraNew::Interface::Required::FD_Res) {
    cout << "ERROR: Either invalid RHS vector or call made from other than "
         << "NOX::EpetraNew::FiniteDifference to OffBlock fill !!" << endl;
    throw "OffBlock_Manager ERROR";
  }

  GenericEpetraProblem &problemEq = myManager.getProblem(problemEqId),
                       &problemVar = myManager.getProblem(problemVarId);

  // Copy relevant part of incoming solnVector into solution vector of
  // problemVarId
  Epetra_Vector &probVarSoln = problemVar.getSolution();
  myManager.copyCompositeToVector(*solnVector, problemVarId, probVarSoln);
  problemEq.doTransfer(); // This does all transfers for this problem
  myManager.setGroupX(problemEqId);
  myManager.computeGroupF(problemEqId);

  const Epetra_Vector &probEqGrpF = dynamic_cast<const NOX::Epetra::Vector&>
    (myManager.getGroup(problemEqId).getF()).getEpetraVector();
  myManager.copyVectorToComposite(*rhsVector, problemEqId, probEqGrpF);

  return true;
}

NOX::EpetraNew::Group& OffBlock_Manager::getGroup()
{
  if( !group ) {
    cout << "ERROR: Unable to get off-block Group for "
         << "dependence of problem " << problemEqId << " on problem "
         << problemVarId << " !!" << endl;
    throw "Problem_Manager ERROR";
  }

  return( *group );
}

Epetra_CrsMatrix& OffBlock_Manager::getMatrix()
{
  if( !matrixOperator ) {
    cout << "ERROR: Unable to get FDColoring underlying matrix for "
         << "dependence of problem " << problemEqId << " on problem "
         << problemVarId << " !!" << endl;
    throw "Problem_Manager ERROR";
  }

  return( matrixOperator->getUnderlyingMatrix() );
}

void OffBlock_Manager::createFDCobjects()
{
#ifdef HAVE_NOX_EPETRAEXT
  // Create a timer for performance
  Epetra_Time colorTime(*Comm);
  bool verbose = false;

  Epetra_CrsGraph &graph = *AA;

  colorTime.ResetStartTime();

  // Just a dummy for now, but needs to be hooked up correctly.
  Epetra_Vector &compositeVec = myManager.getCompositeSoln();

  mapColoring = new EpetraExt::CrsGraph_MapColoring(verbose);
  colorMap = &(*mapColoring)(graph);
  colorMapIndex = new EpetraExt::CrsGraph_MapColoringIndex(*colorMap);
  columnSet = &(*colorMapIndex)(graph);

  if (MyPID == 0) {
    printf("\n\tTime to color Jacobian # %d (%d) --> %e sec. \n",
                problemEqId,problemVarId,colorTime.ElapsedTime());
    cout << "\nUsing " << colorMap->NumColors() << " colors for "
         << graph.NumMyRows() << " unknowns\n" << endl;
  }

  // Now setup each FDC Jacobian as its own group/linearsystem
  matrixOperator = new NOX::EpetraNew::FiniteDifferenceColoring(
    *offBlockInterface, compositeVec, graph, *colorMap, *columnSet );

  NOX::EpetraNew::Interface::Required& reqInt = 
    dynamic_cast<NOX::EpetraNew::Interface::Required&>(*offBlockInterface);
  NOX::EpetraNew::Interface::Jacobian& jacInt = 
    dynamic_cast<NOX::EpetraNew::Interface::Jacobian&>(*matrixOperator);

  // Here we create a linear system solely for the sake of filling an
  // off-diagonal block Jacobian contribution using FDC.  The nlParams and
  // statusTest are irrelevant and so are taken as that of the overall
  // Problem_Manager object.
  NOX::EpetraNew::LinearSystemAztecOO* tmpLinSys = 
    new NOX::EpetraNew::LinearSystemAztecOO(
      myManager.nlParams->sublist("Printing"),
      myManager.nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      reqInt,
      jacInt,
      *matrixOperator,
      compositeVec );

  NOX::Epetra::Vector tmpNOXVec(compositeVec);

  group = new NOX::EpetraNew::Group(
    myManager.nlParams->sublist("Printing"),
    reqInt,
    tmpNOXVec,
    *tmpLinSys);

//    cout << "\n\t\tHERE IS THE GRAPH ----:" << endl;
//    graph.Print(cout);
#else
    if(MyPID==0)
      cout << "ERROR: Cannot use EpetraExt with this build !!" << endl;
    exit(0);
#endif
}
