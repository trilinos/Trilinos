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
                                                                                

#include "NOX.H"

// Currently, this class can only be used with builds which include EpetraExt !!
#ifdef HAVE_NOX_EPETRAEXT

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

// Simple macro for turning on debug code
#undef DEBUG_BLOCKGRAPH
#ifdef ENABLE_DEBUG_BLOCKGRAPH
  #define DEBUG_BLOCKGRAPH(a) a
#else
  #define DEBUG_BLOCKGRAPH(a)
#endif 

int OffBlock_Manager::idToFind = -1;

OffBlock_Manager::OffBlock_Manager(Problem_Manager& problemMan_, 
		Epetra_CrsGraph& graph_, int probEqId, int probVarId) :
  GenericEpetraProblem(graph_.Comm(), 0),
  problemEqId(probEqId),
  problemVarId(probVarId),
  blockRowMap(0),
  blockColMap(0)
{
 
  setManager(&problemMan_);

  // Assign a meaningful name
  GenericEpetraProblem &problemEq = myManager->getProblem(problemEqId),
                       &problemVar = myManager->getProblem(problemVarId);
  
  string myName = "OffBlock " + problemEq.getName() + " wrt " + problemVar.getName();
  setName(myName);

  // Set our graph (held in base class) after converting from incoming
  // global indices to shifted block indices
  //AA = &graph_;
  AA = &( createBlockGraphFromComposite(graph_) );

  // Create a problem interface to the manager
  offBlockInterface = Teuchos::rcp(new Problem_Interface(*this));

  // Use our graph to create FDC objects needed for off-diagonal block fills
  createFDobjects( true );

  // Reset number of dofs (overwrites base constructor assignment)
  NumMyNodes = 0;

}

OffBlock_Manager::~OffBlock_Manager()
{
  delete AA            ; AA             = 0;

  delete blockRowMap   ; blockRowMap    = 0;
  delete blockColMap   ; blockColMap    = 0;

}

// These methods are needed to allow inheritance from GenericEpetraProblem base

bool OffBlock_Manager::evaluate(
              NOX::Epetra::Interface::Required::FillType flag,
              const Epetra_Vector *solnVector,
              Epetra_Vector *rhsVector)
{
  // Determine if fill call is valid
  if (rhsVector == 0 || flag != NOX::Epetra::Interface::Required::FD_Res) {
    cout << "ERROR: Either invalid RHS vector or call made from other than "
         << "NOX::Epetra::FiniteDifference to OffBlock fill !!" << endl;
    throw "OffBlock_Manager ERROR";
  }

  GenericEpetraProblem &problemEq = myManager->getProblem(problemEqId),
                       &problemVar = myManager->getProblem(problemVarId);

  // Copy relevant part of incoming solnVector into solution vector of
  // problemVarId
  Epetra_Vector &probVarSoln = *(problemVar.getSolution());
  //myManager->copyCompositeToVector(*solnVector, problemVarId, probVarSoln);
  probVarSoln = *solnVector;
  problemEq.doTransfer(); // This does all transfers for this problem
  //problemEq.outputSolutionStatus(cout);
  myManager->setGroupX(problemEqId);
  myManager->computeGroupF(problemEqId);

  const Epetra_Vector & probEqGrpF = dynamic_cast<const NOX::Epetra::Vector&>
    (myManager->getGroup(problemEqId).getF()).getEpetraVector();
  //myManager->copyVectorToComposite(*rhsVector, problemEqId, probEqGrpF);
  *rhsVector = probEqGrpF;

  return true;
}

Teuchos::RefCountPtr<NOX::Epetra::Group> OffBlock_Manager::getGroup()
{
  if( !group.get() ) {
    cout << "ERROR: Unable to get off-block Group for "
         << "dependence of problem " << problemEqId << " on problem "
         << problemVarId << " !!" << endl;
    throw "Problem_Manager ERROR";
  }

  return( group );
}

Epetra_CrsMatrix& OffBlock_Manager::getMatrix()
{
  if( !matrixOperator.get() ) {
    cout << "ERROR: Unable to get FDColoring underlying matrix for "
         << "dependence of problem " << problemEqId << " on problem "
         << problemVarId << " !!" << endl;
    throw "Problem_Manager ERROR";
  }

  return( matrixOperator->getUnderlyingMatrix() );
}

Teuchos::RefCountPtr<Epetra_Vector> OffBlock_Manager::getRowMapVec() const
{
  if( !rowMapVec.get() ) {
    cout << "ERROR: Unable to get Row Map Vector for OffBlock " << getName() << endl;
    throw "Problem_Manager ERROR";
  }

  return( rowMapVec );
}

int OffBlock_Manager::getProblemEqId() const
{
  return problemEqId;
}

int OffBlock_Manager::getProblemVarId() const
{
  return problemVarId;
}

void OffBlock_Manager::convertBlockRowIndicesToComposite(int numIndices, 
                         int * blockIndices, int * compositeIndices)
{ 
  for( int i = 0 ; i < numIndices; ++i )
  {
    // Add a check for index validity ? RWH
    compositeIndices[i] = rowBlockToComposite[ blockIndices[i] ];
  }
}

void OffBlock_Manager::convertBlockColIndicesToComposite(int numIndices, 
                         int * blockIndices, int * compositeIndices)
{ 
  for( int i = 0 ; i < numIndices; ++i )
  {
    // Add a check for index validity, RWH
    compositeIndices[i] = colBlockToComposite[ blockIndices[i] ];
  }
}

void OffBlock_Manager::createFDobjects( bool useColoring )
{
  //Teuchos::RefCountPtr<Epetra_CrsGraph> graph = Teuchos::rcp( AA, false );
  graph = Teuchos::rcp( AA );

  DEBUG_BLOCKGRAPH( cout << "OffBlock_Manager::createFDobjects() : incoming graph --> \n" 
                         << *graph << endl << "\n\tDone." << endl;)

  // Note: We use a vector corresponding to compositeSoln
  //Epetra_Vector & compositeVec = myManager->getCompositeSoln();
  rowMapVec = Teuchos::rcp( new Epetra_Vector(graph->RowMap()) );

  noxVec = Teuchos::rcp( new NOX::Epetra::Vector(*rowMapVec) );

  if( !useColoring )
  {
    // Now setup each FD Jacobian as its own group/linearsystem
    matrixOperator = Teuchos::rcp( new NOX::Epetra::FiniteDifference(
      	myManager->nlParams->sublist("Printing"),
        offBlockInterface, 
  	*noxVec, 
  	//rowMapVec, 
//  	compositeVec, 
  	graph) );
  
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> jacInt = matrixOperator;
  
    // Here we create a linear system solely for the sake of filling an
    // off-diagonal block Jacobian contribution using FDC.  The nlParams and
    // statusTest are irrelevant and so are taken as that of the overall
    // Problem_Manager object.
    linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        myManager->nlParams->sublist("Printing"),
        myManager->nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
        offBlockInterface,
        jacInt,
        matrixOperator,
  	*rowMapVec) );
//  	compositeVec, 
  
    //noxVec = new NOX::Epetra::Vector(*rowMapVec);
//    NOX::Epetra::Vector tmpNOXVec(compositeVec);
  
    group = Teuchos::rcp( new NOX::Epetra::Group(
      myManager->nlParams->sublist("Printing"),
      offBlockInterface,
      *noxVec,
      linearSystem) );

    DEBUG_BLOCKGRAPH(   
      cout << "OffBlock_Manager::createFDobjects .... " << myName << endl
           << "---------------------------------------------------------------"
           << "graph :" << *graph << endl
           << "---------------------------------------------------------------"
           << "rowMapVec :" << *rowMapVec << endl
           << "---------------------------------------------------------------" << endl;)
  }
  else // use FDC
  {

#ifdef HAVE_NOX_EPETRAEXT
    // Create a timer for performance
    Epetra_Time colorTime(*Comm);
    // GREEDY dies when doing a transpose for some reason, RWH
    //EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
    //  EpetraExt::CrsGraph_MapColoring::GREEDY;
    EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
      EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;
    int reordering = 0;
    bool useParallel = false;
    bool distance1 = false;
    int verbose = 0;
  
    DEBUG_BLOCKGRAPH( cout << "OffBlock_Manager::createFDobjects() : incoming graph --> \n" << *graph << endl;)
  
    colorTime.ResetStartTime();
  
    // Just a dummy for now, but needs to be hooked up correctly.
    //Epetra_Vector &compositeVec = *(myManager->getCompositeSoln());
  
    mapColoring   = Teuchos::rcp( new EpetraExt::CrsGraph_MapColoring(algType, reordering, distance1, verbose) );
    colorMap      = Teuchos::rcp( & ((*mapColoring)(*graph)) );
    colorMapIndex = Teuchos::rcp( new EpetraExt::CrsGraph_MapColoringIndex(*colorMap) );
    columnSet     = Teuchos::rcp( (& (*colorMapIndex)( *graph )));
  
    if (MyPID == 0) {
      printf("\n\tTime to color Jacobian # %d (%d) --> %e sec. \n",
                  problemEqId,problemVarId,colorTime.ElapsedTime());
      cout << "\nUsing " << colorMap->NumColors() << " colors for "
           << graph->NumMyRows() << " unknowns\n" << endl;
    }
  
    // Now setup each FDC Jacobian as its own group/linearsystem
    matrixOperator = Teuchos::rcp( new NOX::Epetra::FiniteDifferenceColoring(
        myManager->nlParams->sublist("Printing"),
      	offBlockInterface, 
  	*noxVec, 
  	graph, 
  	colorMap, 
  	columnSet,
  	useParallel,
  	distance1 ) );
  
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> jacInt = matrixOperator;
  
    // Here we create a linear system solely for the sake of filling an
    // off-diagonal block Jacobian contribution using FDC.  The nlParams and
    // statusTest are irrelevant and so are taken as that of the overall
    // Problem_Manager object.
    linearSystem = linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        myManager->nlParams->sublist("Printing"),
        myManager->nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
        offBlockInterface,
        jacInt,
        matrixOperator,
        *rowMapVec ) );
  
    NOX::Epetra::Vector tmpNOXVec(*rowMapVec);
  
    group = Teuchos::rcp( new NOX::Epetra::Group(
      myManager->nlParams->sublist("Printing"),
      offBlockInterface,
      tmpNOXVec,
      linearSystem) );
#else
      if(MyPID==0)
        cout << "ERROR: Cannot use EpetraExt with this build !!" << endl;
      exit(0);
#endif
  }
}

Epetra_CrsGraph & OffBlock_Manager::createBlockGraphFromComposite(Epetra_CrsGraph & globalGraph)
{
  // here we simply go through the graph and assign contiguous indices for
  // all unique and meaningful entries
  
  int rowCount = 0;
  int colCount = 0;
  int   numCols;
  int * indices;
  for( int lRow = 0; lRow < globalGraph.NumMyRows(); ++lRow )
  {
    globalGraph.ExtractMyRowView( lRow, numCols, indices ); // indices are local
    if( numCols > 0 )
    {
      rowBlockToComposite[ rowCount ] = globalGraph.GRID( lRow );
      rowCompositeToBlock[ globalGraph.GRID( lRow ) ] = rowCount++ ;
      DEBUG_BLOCKGRAPH( cout << "[" << lRow << "] ";)
      for( int col = 0; col < numCols; ++col )
      {
        DEBUG_BLOCKGRAPH( cout << globalGraph.GCID(indices[col]) << "  ";)
        if( indices[col] >= colCount )
        {
          colBlockToComposite[ colCount ] = globalGraph.GCID( indices[col]  );
          colCompositeToBlock[ globalGraph.GCID( indices[col]  ) ] = colCount++ ;
          
        }
      }
    }
    DEBUG_BLOCKGRAPH( cout << endl;)
  }

  // Now create the block-sized row and column maps
  int numGlobalRows = -1;
  int numGlobalCols = -1;

  blockRowMap = new Epetra_Map ( numGlobalRows, rowBlockToComposite.size(), 0, globalGraph.Comm() );
  blockColMap = new Epetra_Map ( numGlobalCols, colBlockToComposite.size(), 0, globalGraph.Comm() );
 
  DEBUG_BLOCKGRAPH( cout << "\n----> Block-sized Row Map : " << *blockRowMap << endl;)
  DEBUG_BLOCKGRAPH( cout << "\n----> Block-sized Col Map : " << *blockColMap << endl;)

  Epetra_CrsGraph * blockGraph = new Epetra_CrsGraph( Copy, *blockRowMap, *blockColMap, 0, false);

  // Finally, construct a block-sized graph corresponding to a shifted globalGraph
  for( int lRow = 0; lRow < globalGraph.NumMyRows(); ++lRow )
  {
    globalGraph.ExtractMyRowView( lRow, numCols, indices ); // indices are local
    if( numCols > 0 )
    {
      int blockLocalRow = rowCompositeToBlock[ globalGraph.GRID( lRow ) ];
      vector<int> newIndices;
      newIndices.resize(numCols);
      for( int col = 0; col < numCols; ++col )
      {
        newIndices[col] = colCompositeToBlock[ globalGraph.GCID( indices[col]  ) ];
      }

      blockGraph->InsertMyIndices( blockLocalRow, numCols, &(newIndices[0]) );

    }
    DEBUG_BLOCKGRAPH( cout << endl;)
  }

  
  DEBUG_BLOCKGRAPH( cout << "\n----> Block-Graph : " << *blockGraph << endl;)
  
  blockGraph->FillComplete();

  return *blockGraph; 
}

#endif 
