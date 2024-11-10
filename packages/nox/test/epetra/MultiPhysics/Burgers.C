// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"
#include <math.h>

#include "Problem_Manager.H"
#include "Burgers.H"

// Constructor - creates the Epetra objects (maps and vectors)
Burgers::Burgers(Epetra_Comm& comm, int numGlobalNodes, std::string name_) :
  GenericEpetraProblem(comm, numGlobalNodes, name_),
  xFactor(5.0),
  viscosity(0.100),
  xmin(0.0),
  xmax(1.0),
  dt(1.0e-3)
{

  // Create the nodal coordinates
  xptr = Teuchos::rcp( new Epetra_Vector(*StandardMap) );
  double Length= xmax - xmin;
  dx=Length/((double) NumGlobalNodes-1);
  for ( int i = 0; i < NumMyNodes; ++i )
    (*xptr)[i] = xmin + dx*((double) StandardMap->MinMyGID()+i);

  // Create extra vector needed for this transient problem
  oldSolution = new Epetra_Vector(*StandardMap);
  exactSolution = new Epetra_Vector(*StandardMap);

  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  initializeSoln();

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  AA = Teuchos::rcp( new Epetra_CrsGraph(Copy, *StandardMap, 5) );
  generateGraph();

  // Create a second matrix using graph of first matrix - this creates a
  // static graph so we can refill the new matirx after FillComplete()
  // is called.
  A = Teuchos::rcp(new Epetra_CrsMatrix (Copy, *AA));
  A->FillComplete();

  // Create the Importer needed for FD coloring
  ColumnToOverlapImporter = new Epetra_Import(A->ColMap(),*OverlapMap);

}

// Destructor
Burgers::~Burgers()
{
  delete oldSolution;
  delete exactSolution;
  delete ColumnToOverlapImporter;
}

// Reset function
void Burgers::reset(const Epetra_Vector& x)
{
  *oldSolution = x;
}

// Empty Reset function
void Burgers::reset()
{
  std::cout << "WARNING: reset called without passing any update vector !!"
       << std::endl;
}

// Set initialSolution to desired initial condition
bool Burgers::initializeSoln()
{
  Epetra_Vector& x = *xptr;

  double arg;
  for( int i = 0; i < NumMyNodes; ++i)
  {
    arg = ( 20.0 * x[i] - 10.0 ) / xFactor;
    (*initialSolution)[i] = (1.0 - ( exp(arg) - exp(-arg) ) /
                                   ( exp(arg) + exp(-arg) )) - 1.0;
//    (*initialSolution)[i] = 1.0*sin(2.0*pi*x[i]);
  }

  *oldSolution = *initialSolution;

  return true;
}

// Matrix and Residual Fills
bool Burgers::evaluate(
             NOX::Epetra::Interface::Required::FillType fillType,
             const Epetra_Vector* soln,
             Epetra_Vector* rhs)
{
  bool fillRes = false;
  bool fillJac = false;

  // We currently make no provision for combined Residual and Jacobian fills
  // i.e. it is either or

  if( NOX::Epetra::Interface::Required::Jac == fillType )
    fillJac = true;
   else
   {
    fillRes = true;
    assert( NULL != rhs );
   }

  int numDep = depProblems.size();

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  std::vector<Epetra_Vector*> dep(numDep);
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
  {
    dep[i]->Import(*( (*(depSolutions.find(depProblems[i]))).second ),
                   *Importer, Insert);
    //cout << "depSoln[" << i << "] :" << dep[i] << std::endl;
  }
  xvec.Import(*xptr, *Importer, Insert);
  if( NOX::Epetra::Interface::Required::FD_Res == fillType )
    // Overlap vector for solution received from FD coloring, so simply reorder
    // on processor
    u.Export(*soln, *ColumnToOverlapImporter, Insert);
  else // Communication to Overlap vector is needed
    u.Import(*soln, *Importer, Insert);

  int OverlapNumMyNodes = OverlapMap->NumMyElements();

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;

  int row, column;
  double jac;
  double xx[2];
  double uu[2];
  double uuold[2];
  std::vector<double*> ddep(numDep);
  for( int i = 0; i<numDep; i++)
    ddep[i] = new double[2];
  Basis basis;

  int id_temp; // Index for needed dependent Temperature vector

  //
  map<string, int>::iterator id_ptr = nameToMyIndex.find("Temperature");
  if( id_ptr == nameToMyIndex.end() ) {
    std::cout << "WARNING: Burgers (\"" << myName << "\") could not get "
         << "vector for problem \"Temperature\" !!" << std::endl;
    throw "Burgers ERROR";
  }
  else
    id_temp = (*id_ptr).second;
    //

  // Zero out the objects that will be filled
  if( fillJac )
    A->PutScalar(0.0);
  if( fillRes )
    rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {

    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      uu[0]=u[ne];
      uu[1]=u[ne+1];
      uuold[0]=uold[ne];
      uuold[1]=uold[ne+1];
      for( int i = 0; i<numDep; i++ )
      {
        ddep[i][0] = (*dep[i])[ne];
        ddep[i][1] = (*dep[i])[ne+1];
      }
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uuold, ddep);

      // Loop over Nodes in Element
      for ( int i = 0; i < 2; ++i )
      {
    row=OverlapMap->GID(ne+i);
    //printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
    //     MyPID, row, ne+i,StandardMap.MyGID(row));
    if (StandardMap->MyGID(row)) {
      if ( fillRes )
          {
        (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
              +basis.wt*basis.dx*(
                (basis.uu-basis.uuold)/dt*basis.phi[i] +
                (1.0/(basis.dx*basis.dx))*viscosity*pow(1.0*(basis.ddep[id_temp]-0.2),1.5)
                //(1.0/(basis.dx*basis.dx))*viscosity
                  *basis.duu*basis.dphide[i] -
                (1.0/basis.dx)*0.5*basis.uu*basis.uu*basis.dphide[i]);
      }
    }
    // Loop over Trial Functions
    if ( fillRes )
        {
      for(int j = 0; j < 2; ++j)
          {
        if (StandardMap->MyGID(row)) {
          column=OverlapMap->GID(ne+j);
              jac=basis.wt*basis.dx*(
                           basis.phi[j]/dt*basis.phi[i] +
                           (1.0/(basis.dx*basis.dx))*
                             basis.dphide[j]*basis.dphide[i] -
                            8.0/xFactor/xFactor*
                             (2.0*basis.uu-3.0*basis.uu*basis.uu)*
                             basis.phi[j]*basis.phi[i]);
          A->SumIntoGlobalValues(row, 1, &jac, &column);
        }
      }
    }
      }
    }
  }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(xmin)=1
  if( MyPID == 0 )
  {
    if( fillRes )
      (*rhs)[0]= (*soln)[0] - 1.0;
    if( fillJac )
    {
      column=0;
      jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }
  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(xmax)=0
  if( MyPID == NumProc-1 )
  {
    if ( fillRes )
      (*rhs)[NumMyNodes-1]= (*soln)[OverlapNumMyNodes-1] - (-1.0);
    if( fillJac )
    {
      row=NumGlobalNodes-1;
      column=row;
      jac=1.0;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      jac=0.0;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();

  if ( fillJac )
    A->FillComplete();

  // Cleanup
  for( int i = 0; i < numDep; ++i )
  {
    delete [] ddep[i];
    delete    dep[i];
  }

  return true;
}

Teuchos::RCP<Epetra_Vector> Burgers::getSolution()
{
  return initialSolution;
}

Epetra_Vector& Burgers::getExactSoln(double time)
{
  Epetra_Vector& x = *xptr;

  for(int i=0; i<NumMyNodes; i++)
    (*exactSolution)[i] = (1.0 - tanh( (x[i]-2.0*time/xFactor)/xFactor )) / 2.0;

  return *exactSolution;
}

Epetra_Vector& Burgers::getOldSoln()
{
  return *oldSolution;
}

void Burgers::generateGraph()
{

  // Declare required variables
  int i;
  int row, column;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();
  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {

    // Loop over Nodes in Element
    for (i=0; i< 2; i++) {
      row=OverlapMap->GID(ne+i);

      // Loop over Trial Functions
      for( int j = 0; j < 2; ++j )
      {
    // If this row is owned by current processor, add the index
    if (StandardMap->MyGID(row))
        {
      column=OverlapMap->GID(ne+j);
      AA->InsertGlobalIndices(row, 1, &column);
    }
      }
    }
  }
  AA->FillComplete();
}
