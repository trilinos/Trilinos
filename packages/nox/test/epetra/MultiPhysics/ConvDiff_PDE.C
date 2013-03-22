//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
            double bcWeight_             ,
            double xmin_                 ,
            double xmax_                 ,
            double Tleft_                ,
            double Tright_               ,
            int NumGlobalUnknowns_       ,
            std::string name_                   ) :
  GenericEpetraProblem(Comm_, NumGlobalUnknowns_, name_),
  xmin            ( xmin_      )   ,
  xmax            ( xmax_      )   ,
  Tleft           ( Tleft_     )   ,
  Tright          ( Tright_    )   ,
  peclet          ( peclet_    )   ,
  radiation       ( radiation_ )   ,
  kappa           ( kappa_     )   ,
  bcWeight        ( bcWeight_  )   ,
  expandJacobian  ( true      )    ,
  depProbPtr      ( NULL       )
{
  // Create mesh and solution vectors

  // We first initialize the mesh and then the solution since the latter can depend on the mesh.

  xptr = Teuchos::rcp( new Epetra_Vector(*StandardMap) );
  dx   = (xmax - xmin) / ( (double) NumGlobalNodes - 1 );

  for( int i = 0; i < NumMyNodes; ++i ) 
    (*xptr)[i]=xmin + dx*((double) StandardMap->MinMyGID()+i);

  // Create extra vector needed for transient problem interface
  oldSolution = Teuchos::rcp( new Epetra_Vector(*StandardMap) );

  // Create and initialize (using default provided) the solution vector
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  initializeSolution();

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  AA = Teuchos::rcp( new Epetra_CrsGraph(Copy, *StandardMap, 0) );
  generateGraph();

#ifdef DEBUG
  AA->Print(cout);
#endif

  // Create a matrix using the graph just created - this creates a
  // static graph so we can refill the new matirx after FillComplete()
  // is called.
  A = Teuchos::rcp(new Epetra_CrsMatrix (Copy, *AA));
  A->FillComplete();

  // Create the Importer needed for FD coloring
  ColumnToOverlapImporter = Teuchos::rcp( new Epetra_Import(A->ColMap(),*OverlapMap) );

}

//-----------------------------------------------------------------------------

// Destructor
ConvDiff_PDE::~ConvDiff_PDE()
{
}

//-----------------------------------------------------------------------------

// Initialize method
void 
ConvDiff_PDE::initialize()
{
  // Verify we have only one dependent problem and that it is of the appropriate type
  if( 1 < depProblems.size() )
  {
    std::string msg = "ERROR: ConvDiff_PDE::initialize : Problem \"" + myName
                    + "\" depends on more than one other problem.";
    throw msg;
  }

  GenericEpetraProblem & depProb = myManager->getProblem( depProblems[0] );
  depProbPtr = dynamic_cast<ConvDiff_PDE *>(&depProb);
  if( NULL == depProbPtr )
  {
    std::string msg = "ERROR: ConvDiff_PDE::initialize : Dependent problem \"" 
                    + depProb.getName() + "\" is not of type ConvDiff_PDE.";
    throw msg;
  }

  // Determine relative location of coupling interface
  if( (*xptr)[0] < depProbPtr->getMesh()[0] )
    myInterface = RIGHT;
  else
    myInterface = LEFT;

  // Set exact interface temperature accordingly
  if( LEFT == myInterface )
    T1_exact = Tleft ;
  else
    T1_exact = Tright;

  // Choose appropriate element and nodes for interfacial BC enforcement
  int OverlapNumMyNodes = OverlapMap->NumMyElements(); // this may break in parallel
  if( LEFT == myInterface )
  {
    interface_elem = 0                     ;
    local_node     = 0                     ;
    interface_node = 0                     ;
    opposite_node  = OverlapNumMyNodes - 1 ;
    dirScale       = 1.0                   ;
  }
  else
  {
    interface_elem = OverlapNumMyNodes - 2 ;
    local_node     = 1                     ;
    interface_node = OverlapNumMyNodes - 1 ;
    opposite_node  = 0                     ;
    dirScale       = -1.0                  ;
  }

  // Use this method to compute and output the analytic solution and its first derivative
  exactSolution     = Teuchos::rcp( new Epetra_Vector(*initialSolution) );
  dTdx              = Teuchos::rcp( new Epetra_Vector(*initialSolution) );
  Epetra_Vector & x = *xptr;

  if( 1.e-20 < fabs(peclet) )
  {
    for( int i = 0; i < NumMyNodes; ++i ) 
    {
      (*exactSolution)[i] = (T1_exact - Tleft*exp(peclet) + (Tleft - T1_exact)*exp(peclet*(x[i]-xmin))) /
                            ( 1.0 - exp(peclet) );
      (*dTdx         )[i] = peclet * (Tleft - T1_exact)*exp(peclet*(x[i]-xmin)) /
                            ( 1.0 - exp(peclet) );
    }
  }
  else
  {
    for( int i = 0; i < NumMyNodes; ++i ) 
    {
      (*exactSolution)[i] = (Tright - T1_exact)*(x[i] - xmax) + Tright  ;
      (*dTdx         )[i] =  Tright - T1_exact                          ;
    }
  }

  ostringstream sval; 
  sval << myId << flush; 
  std::string fileName1 = "analytic_" + sval.str(); 
  std::string fileName2 = "dTdx_"     + sval.str(); 

  std::ofstream outFile1(fileName1.c_str());
  std::ofstream outFile2(fileName2.c_str());
  if( !outFile1 || !outFile2 )
  {
    std::string msg = "ERROR: Could not open one of files \"" + fileName1 + "\"" + " or \""
                      + fileName2 + "\"";
    throw msg;
  }

  for( int i = 0; i < NumMyNodes; ++i ) 
  {
    outFile1 << i << "  " << x[i] << "  " << (*exactSolution)[i] << std::endl;
    outFile2 << i << "  " << x[i] << "  " << (*dTdx         )[i] << std::endl;
  }
}

//-----------------------------------------------------------------------------

// Set initialSolution to desired initial condition
void 
ConvDiff_PDE::initializeSolution(double val)
{
  // Aliases for convenience
  Epetra_Vector& soln = *initialSolution;

  soln.PutScalar(val);

  //Epetra_Vector perturb(soln);
  //perturb.Random();
  //perturb.Scale(0.20);

  //soln.Update(1.0, perturb, 1.0);
  
  *oldSolution = soln;
} 

//-----------------------------------------------------------------------------

// Matrix and Residual Fills
bool 
ConvDiff_PDE::evaluate(
                    NOX::Epetra::Interface::Required::FillType flag,
		    const Epetra_Vector * soln, 
		    Epetra_Vector * rhs)
{

  if( rhs == 0 ) 
  {
    std::string msg = "ERROR: ConvDiff_PDE::evaluate : callback appears to be other than a residual fill.  Others are not support for this type.";
    throw msg;
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
    (*dep[i]).Import(*( (*(depSolutions.find(depProblems[i]))).second ), *Importer, Insert);

  xvec.Import(*xptr, *Importer, Insert);

  if( flag == NOX::Epetra::Interface::Required::FD_Res)
    // Overlap vector for solution received from FD coloring, so simply reorder
    // on processor
    u.Export(*soln, *ColumnToOverlapImporter, Insert);
  else // Communication to Overlap vector is needed
    u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int OverlapNumMyNodes = OverlapMap->NumMyElements();

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;

  int row;

  double * xx    = new double[2];
  double * uu    = new double[2]; 
  double * uuold = new double[2];
  std::vector<double*> ddep(numDep);
  for( int i = 0; i < numDep; ++i)
    ddep[i] = new double[2];

  Basis basis;

  // Bundle up the dependent variables in the way needed for computing
  // the source terms of each reaction
  map<string, double*> depVars;
  depVars.insert( pair< std::string, double*>(getName(), uu) );
  for( int i = 0; i < numDep; ++i )
    depVars.insert( pair<string, double*>(myManager->getProblemName(depProblems[i]), ddep[i]) );

  // Zero out the objects that will be filled
  rhs->PutScalar(0.0);

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
      for( int i = 0; i < numDep; ++i ) 
      {
        ddep[i][0] = (*dep[i])[ne];
        ddep[i][1] = (*dep[i])[ne+1];
      }
      // Calculate the basis function and variables at the gauss points
      basis.getBasis(gp, xx, uu, uuold, ddep);

      // Loop over Nodes in Element
      for( int i = 0; i < 2; ++i )
      {
	row = OverlapMap->GID(ne+i);
	if( StandardMap->MyGID(row) ) 
        {
          (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))] +=
            + basis.wt * basis.dx
            * ( peclet * (basis.duu / basis.dx) * basis.phi[i] 
            +   kappa * (1.0/(basis.dx*basis.dx)) * basis.duu * basis.dphide[i] );
	}
      }
    }
  } 

  //if( NOX::Epetra::Interface::Required::Residual == flag )
  //{
  //  int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
  //  std::cout << "\t\"" << myName << "\" u[0] = " << (*soln)[0] 
  //       << "\tu[N] = " << (*soln)[lastDof] << std::endl;
  //  std::cout << "\t\"" << myName << "\" RHS[0] = " << (*rhs)[0] 
  //       << "\tRHS[N] = " << (*rhs)[lastDof] << std::endl << std::endl;
  //}


  // Apply BCs

  computeHeatFlux( soln );

  double bcResidual = bcWeight         * (myFlux - depProbPtr->getHeatFlux()                 ) -
                      (1.0 - bcWeight) * (u[interface_node] - depProbPtr->getInterfaceTemp() );

  int lastDof = StandardMap->LID(StandardMap->MaxAllGID());

  // "Left" boundary
  if( LEFT == myInterface ) // this may break in parallel
  {
    (*rhs)[0]       = bcResidual;
    (*rhs)[lastDof] = (*soln)[lastDof] - Tright;
  }
  // "Right" boundary
  else
  {
    (*rhs)[0]       = (*soln)[0] - Tleft;
    (*rhs)[lastDof] = bcResidual;
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->FillComplete();

#ifdef DEBUG
  std::cout << "For residual fill :" << std::endl << *rhs << std::endl;
#endif

  // Cleanup
  for( int i = 0; i < numDep; ++i)
  {
    delete [] ddep[i];
    delete     dep[i];
  }

  delete [] xx    ;
  delete [] uu    ;
  delete [] uuold ;

  return true;
}

//-----------------------------------------------------------------------------

// Initialize method
void 
ConvDiff_PDE::getOffBlockIndices( map<int, std::vector<int> > & indices )
{

  // Get the equation map from the other problem
  Epetra_Map & otherMap = *(depProbPtr->StandardMap);
  int otherMinEqID = otherMap.MinAllGID();
  int otherMaxEqID = otherMap.MaxAllGID();

  std::vector<int> offBlockColumns(2);

  if( LEFT == myInterface )
  {
    offBlockColumns[0]      = otherMaxEqID - 1 ;
    offBlockColumns[1]      = otherMaxEqID     ;
    indices[interface_node] = offBlockColumns  ;
  }
  else
  {
    offBlockColumns[0]      = otherMinEqID     ;
    offBlockColumns[1]      = otherMinEqID + 1 ;
    indices[interface_node] = offBlockColumns  ;
  }

  return;
}

//-----------------------------------------------------------------------------

void 
ConvDiff_PDE::prepare_data_for_transfer()
{
  // Redirect to compute heat fluxes
  return computeHeatFlux();
}

//-----------------------------------------------------------------------------

// A fill specialized to the single node at the coupling interface
void 
ConvDiff_PDE::computeHeatFlux( const Epetra_Vector * soln )
{
  
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
    (*dep[i]).Import(*( (*(depSolutions.find(depProblems[i]))).second ), *Importer, Insert);

  xvec.Import(*xptr, *Importer, Insert);

  if( NULL == soln )
    u.Import(*initialSolution, *Importer, Insert);
  else
    u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int row;
  double * xx = new double[2];
  double * uu = new double[2]; 
  double * uuold = new double[2];
  std::vector<double*> ddep(numDep);
  for( int i = 0; i < numDep; ++i)
    ddep[i] = new double[2];

  Basis basis;

  // Bundle up the dependent variables in the way needed for computing
  // the source terms of each reaction
  map<string, double*> depVars;
  depVars.insert( pair< std::string, double*>(getName(), uu) );
  for( int i = 0; i < numDep; ++i )
    depVars.insert( pair<string, double*>(myManager->getProblemName(depProblems[i]), ddep[i]) );

  myFlux = 0.0;

  // Loop Over Gauss Points
  for( int gp = 0; gp < 2; ++gp ) 
  {
    // Get the solution and coordinates at the nodes 
    xx[0]=xvec[interface_elem];
    xx[1]=xvec[interface_elem+1];
    uu[0] = u[interface_elem];
    uu[1] = u[interface_elem+1];
    uuold[0] = uold[interface_elem];
    uuold[1] = uold[interface_elem+1];
    for( int i = 0; i < numDep; ++i ) 
    {
      ddep[i][0] = (*dep[i])[interface_elem];
      ddep[i][1] = (*dep[i])[interface_elem+1];
    }

    // Calculate the basis function and variables at the gauss points
    basis.getBasis(gp, xx, uu, uuold, ddep);

    row = OverlapMap->GID( interface_elem + local_node );

    if( StandardMap->MyGID(row) ) 
    {
      myFlux += 
        + basis.wt * basis.dx
        * ( peclet * (basis.duu / basis.dx) * basis.phi[local_node] 
        +   kappa * (1.0/(basis.dx*basis.dx)) * basis.duu * basis.dphide[local_node] );
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  // Cleanup
  for( int i = 0; i < numDep; ++i)
  {
    delete [] ddep[i];
    delete     dep[i];
  }

  delete [] xx    ;
  delete [] uu    ;
  delete [] uuold ;

  //int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
  //cout << "\t\"" << myName << "\" u[0] = " << u[0] 
  //     << "\tu[N] = " << u[lastDof] << std::endl;
  //cout << u << std::endl;
  //cout << "\t\"" << myName << "\" myFlux = " << myFlux << std::endl << std::endl;

  // Scale domain integration according to interface position
  myFlux *= dirScale;

  // Now add radiation contribution to flux
  myFlux += radiation * ( pow(u[interface_node], 4) - pow(u[opposite_node], 4) );

  return;
}

//-----------------------------------------------------------------------------

double 
ConvDiff_PDE::getInterfaceTemp()
{
  
  if( LEFT == myInterface )
    return (*initialSolution)[0];
  else
  {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    return (*initialSolution)[lastDof];
  }
}

//-----------------------------------------------------------------------------

void 
ConvDiff_PDE::process_transferred_data()
{
  // Now that each problem has computed its fluxes, get these as well as boundary
  // temperature from problems on which we depend

  return;
}

//-----------------------------------------------------------------------------

void 
ConvDiff_PDE::outputStatus( std::ostream & os ) 
{
  
  std::string location = ( myInterface == LEFT ) ? "Left" : "Right";

  os << "\"" << myName << "\" couples at an interface with Problem \"" << depProbPtr->getName()
     << "\" on the " << location << std::endl;

  return;
}

//-----------------------------------------------------------------------------

void 
ConvDiff_PDE::generateGraph()
{
  
  // Declare required variables
  int row, column;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();
  
  // Loop Over # of Finite Elements on Processor
  for( int ne = 0; ne < OverlapNumMyNodes - 1; ++ne) 
  {
    // Loop over Nodes in Element
    for( int i = 0; i < 2; ++i ) 
    {
      // If this node is owned by current processor, add indices
      if (StandardMap->MyGID(OverlapMap->GID(ne+i))) 
      {
        // Loop over unknowns in Node
        row = OverlapMap->GID(ne+i);

        // Loop over supporting nodes
        for( int j = 0; j < 2; ++j) 
        {
          // Loop over unknowns at supporting nodes
          column = OverlapMap->GID(ne+j);
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

void 
ConvDiff_PDE::doTransfer()
{
  
  depProbPtr->computeHeatFlux();

  return;
}

//-----------------------------------------------------------------------------

double 
ConvDiff_PDE::computeAnalyticInterfaceTemp(
                  double radiation      ,
                  double T_left         ,
                  double T_right        ,
                  double kappa          ,
                  double peclet           )
{
  
  int    nIters         = 0;            // iteration counter
  double T_int          = 1.0;          // initial guess
  double tol            = 1.e-12;       // solve tolerance
  double residual       = radiation*( pow(T_int, 4) - pow(T_right, 4) ) 
                          + (T_left - T_int)*peclet*exp(peclet)/(1.0 - exp(peclet))
                          - kappa*(T_right - T_int);
  double dfdT           = 0.0;

  // Simple Newton loop
  while( tol < fabs(residual) )
  {
    ++nIters;

    dfdT = radiation*4.0*pow(T_int,3)
           - peclet*exp(peclet)/(1.0 - exp(peclet))
           + kappa;

    if( 1.e-15 > fabs(dfdT) )
    {
      std::cout << "ConvDiff_PDE::computeAnalyticInterfaceTemp:\n"
              "Warning: Obtained near-zero derivtative. Aborting calculation." << std::endl;
      return(T_int);
    }

    T_int = T_int - residual/dfdT;

    residual = radiation*( pow(T_int, 4) - pow(T_right, 4) ) 
               + (T_left - T_int)*peclet*exp(peclet)/(1.0 - exp(peclet))
               - kappa*(T_right - T_int);
  }

  std::cout << "Analytic interfacial temperature = " << setprecision(8) << T_int << std::endl;
  std::cout << "Residual = " << residual << " in " << nIters << " iterations." << std::endl;

  return T_int;
}

//-----------------------------------------------------------------------------
