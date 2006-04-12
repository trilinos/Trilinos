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
                                                                                
#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"

#include "Problem_Manager.H"
#include "ConvDiff_PDE.H"

// Constructor - creates the Epetra objects (maps and vectors) 
ConvDiff_PDE::ConvDiff_PDE(
            Epetra_Comm& Comm_           ,
            double peclet_               ,
            double radiation_            ,
            double kappa_                ,
            double xmin_                 ,
            double xmax_                 ,
            double Tleft_                ,
            double Tright_               ,
            int NumGlobalUnknowns_       ,
            string name_                   ) :
  GenericEpetraProblem(Comm_, NumGlobalUnknowns_, name_),
  peclet        (peclet_    )   ,
  radiation     (radiation_ )   ,
  kappa         (kappa_     )   ,
  xmin          (xmin_      )   ,
  xmax          (xmax_      )   ,
  Tleft         (Tleft_     )   ,
  Tright        (Tright_    )
{
  // Create mesh and solution vectors

  // We first initialize the mesh and then the solution since the latter
  // can depend on the mesh.
  xptr = new Epetra_Vector(*StandardMap);
  double Length = xmax - xmin;
  dx = Length / ( (double) NumGlobalNodes - 1 );

  for( int i = 0; i < NumMyNodes; ++i ) 
    (*xptr)[i]=xmin + dx*((double) StandardMap->MinMyGID()+i);

  // Create extra vector needed for this transient problem
  oldSolution = new Epetra_Vector(*StandardMap);

  // Next we create and initialize (using default provided) the solution vector
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  initializeSolution();

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 0);
  generateGraph();

#ifdef DEBUG
  AA->Print(cout);
#endif

  // Create a matrix using the graph just created - this creates a
  // static graph so we can refill the new matirx after FillComplete()
  // is called.
  A = Teuchos::rcp(new Epetra_CrsMatrix (Copy, *AA));new Epetra_CrsMatrix (Copy, *AA);
  A->FillComplete();

  // Create the Importer needed for FD coloring
  ColumnToOverlapImporter = new Epetra_Import(A->ColMap(),*OverlapMap);
}

//-----------------------------------------------------------------------------

// Destructor
ConvDiff_PDE::~ConvDiff_PDE()
{
  delete AA; AA = 0;
  delete xptr; xptr = 0;
  delete oldSolution; oldSolution = 0;
  delete ColumnToOverlapImporter; ColumnToOverlapImporter = 0;
}

//-----------------------------------------------------------------------------

// Initialize method
void 
ConvDiff_PDE::initialize()
{
  // Use this method to compute and output the analytic solution
  Epetra_Vector * exactSolution = new Epetra_Vector(*initialSolution);
  Epetra_Vector & x = *xptr;

  if( 1.e-20 < fabs(peclet) )
  {
    for( int i = 0; i < NumMyNodes; ++i ) 
      (*exactSolution)[i] = (Tright - Tleft*exp(peclet) + (Tleft - Tright)*exp(peclet*(x[i]-xmin))) /
                            ( 1.0 - exp(peclet) );
  }
  else
  {
    for( int i = 0; i < NumMyNodes; ++i ) 
      (*exactSolution)[i] = (Tright - Tleft)*(x[i] - xmax) + Tright;
  }

  ostringstream sval; 
  sval << myId << flush; 
  std::string fileName = "analytic_" + sval.str(); 

  ofstream outFile(fileName.c_str());
  if( !outFile )
  {
    std::string msg = "ERROR: Could not open file \"" + fileName + "\"";
    throw msg;
  }

  for( int i = 0; i < NumMyNodes; ++i ) 
    outFile << i << "  " << x[i] << "  " << (*exactSolution)[i] << endl;

  delete exactSolution;
}

//-----------------------------------------------------------------------------

// Set initialSolution to desired initial condition
void 
ConvDiff_PDE::initializeSolution(double val)
{
  // Aliases for convenience
  Epetra_Vector& soln = *initialSolution;
  Epetra_Vector& x = *xptr;

  soln.PutScalar(val);
  
  *oldSolution = soln;
} 

//-----------------------------------------------------------------------------

// Matrix and Residual Fills
bool 
ConvDiff_PDE::evaluate(
                    NOX::Epetra::Interface::Required::FillType flag,
		    const Epetra_Vector* soln, 
		    Epetra_Vector* tmp_rhs)
{
  // Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  if (tmp_rhs != 0) 
  {
    fillF = true;
    rhs = tmp_rhs;
  }
  else 
  {
    fillMatrix = true;
  }
  
  int numDep = depProblems.size();

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  vector<Epetra_Vector*> dep(numDep);
  for( int i = 0; i < numDep; ++i)
    dep[i] = new Epetra_Vector(*OverlapMap);

  Epetra_Vector xvec(*OverlapMap);

  // Export Solution to Overlap vector
  // If the vector to be used in the fill is already in the Overlap form,
  // we simply need to map on-processor from column-space indices to
  // OverlapMap indices. Note that the old solution is simply fixed data that
  // needs to be sent to an OverlapMap (ghosted) vector.  The conditional
  // treatment for the current soution vector arises from use of
  // FD coloring in parallel.

  uold.Import(*oldSolution, *Importer, Insert);

  for( int i = 0; i < numDep; ++i )
    (*dep[i]).Import(*( (*(depSolutions.find(depProblems[i]))).second ), *Importer, Insert);

  xvec.Import(*xptr, *Importer, Insert);

  if( flag == NOX::Epetra::Interface::Required::FD_Res)
    // Overlap vector for solution received from FD coloring, so simply reorder
    // on processor
    u.Export(*soln, *ColumnToOverlapImporter, Insert);
  else // Communication to Overlap vector is needed
    u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int j,ierr;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;

  int row, column;
  double jac;
  double xx[2];
  double uu[2]; 
  double uuold[2];
  vector<double*> ddep(numDep);
  for( int i = 0; i < numDep; ++i)
    ddep[i] = new double[2];
  double *srcTerm = new double[2];
  Basis basis;

  // Bundle up the dependent variables in the way needed for computing
  // the source terms of each reaction
  map<string, double*> depVars;
  depVars.insert( pair< string, double*>(getName(), uu) );
  for( int i = 0; i < numDep; ++i )
    depVars.insert( pair<string, double*>(myManager->getName(depProblems[i]), ddep[i]) );

  // Do a check on this fill
//  map<string, double*>::iterator iter;
//  for( iter = depVars.begin(); iter != depVars.end(); iter++)
//    cout << "Inserted ... " << iter->first << "\t" << iter->second << endl;
//  cout << "--------------------------------------------------" << endl;
//  for( iter = depVars.begin(); iter != depVars.end(); iter++)
//	  cout << iter->first << "\t" << (iter->second)[0] << ", " 
//               << (iter->second)[1] << endl;
//  cout << "--------------------------------------------------" << endl;

  // Zero out the objects that will be filled
  if ( fillMatrix ) A->PutScalar(0.0);
  if ( fillF ) rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for( int ne = 0; ne < OverlapNumMyNodes-1; ++ne )
  {
    // Loop Over Gauss Points
    for( int gp = 0; gp < 2; ++gp ) 
    {
      // Get the solution and coordinates at the nodes 
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      uu[0] = u[ne];
      uu[1] = u[ne+1];
      uuold[0] = uold[ne];
      uuold[1] = uold[ne+1];
      for( int i = 0; i<numDep; i++ ) 
      {
        ddep[i][0] = (*dep[i])[ne];
        ddep[i][1] = (*dep[i])[ne+1];
      }
      // Calculate the basis function and variables at the gauss points
      basis.getBasis(gp, xx, uu, uuold, ddep);

      // Loop over Nodes in Element
      for (int i=0; i< 2; i++) 
      {
	row=OverlapMap->GID(ne+i);
	if (StandardMap->MyGID(row)) 
        {
	  if ( fillF ) 
          {
	    (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))] +=
	      + basis.wt * basis.dx
	      * ( peclet * (basis.duu / basis.dx) * basis.phi[i] 
              +   kappa * (1.0/(basis.dx*basis.dx)) * basis.duu * basis.dphide[i] );
	  }
	}
	// Loop over Trial Functions
	if ( fillMatrix ) 
        {
          // No-op for now - can be filled in for analytical Jacoibian terms
	}
      }
    }
  } 

  // Apply BCs

  // "Left" boundary
  if (MyPID==0) 
  {
    if ( fillF )
      (*rhs)[0]= (*soln)[0] - Tleft;
    if ( fillMatrix ) 
    {
      int column=0;
      double jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }
  // "Right" boundary
  if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 ) 
  {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    if ( fillF )
      (*rhs)[lastDof] = (*soln)[lastDof] - Tright;
    if ( fillMatrix ) 
    {
      int row=StandardMap->MaxAllGID();
      int column = row;
      double jac = 1.0;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      jac=0.0;
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->FillComplete();

#ifdef DEBUG
  A->Print(cout);

  if( fillF )
    cout << "For residual fill :" << endl << *rhs << endl;

  if( fillMatrix ) 
  {
    cout << "For jacobian fill :" << endl;
    A->Print(cout);
  }
#endif

  // Cleanup
  for( int i = 0; i < numDep; ++i)
  {
    delete [] ddep[i];
    delete     dep[i];
  }
  delete [] srcTerm;

  return true;
}

//-----------------------------------------------------------------------------

void 
ConvDiff_PDE::generateGraph()
{
  
  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();
  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) 
  {
    // Loop over Nodes in Element
    for (i=0; i<2; i++) 
    {
      // If this node is owned by current processor, add indices
      if (StandardMap->MyGID(OverlapMap->GID(ne+i))) 
      {
        // Loop over unknowns in Node
        row=OverlapMap->GID(ne+i);

        // Loop over supporting nodes
        for(j=0; j<2; j++) 
        {
          // Loop over unknowns at supporting nodes
          column=OverlapMap->GID(ne+j);
          //printf("\t\tWould like to insert -> (%d, %d)\n",row,column);
          AA->InsertGlobalIndices(row, 1, &column);
        }
      }
    }
  }
  AA->FillComplete();
  
  return;
}

//-----------------------------------------------------------------------------

