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
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"

#include "../test/utils/NOX_Epetra_DebugTools.H"
//#include "NOX_Epetra_DebugTools.H"

#include "PelletTransport.H"
#include "MaterialProps.H"

// Constructor - creates the Epetra objects (maps and vectors)
PelletTransport::PelletTransport( int NumGlobalElementsUO2_  , double xminUO2_  , double xmaxUO2_  ,
                                  int NumGlobalElementsHe_   , double xminHe_   , double xmaxHe_   ,
                                  int NumGlobalElementsClad_ , double xminClad_ , double xmaxClad_ ,
                                  Epetra_Comm& Comm_, bool restart_ ) :
  NumGlobalElementsUO2(NumGlobalElementsUO2_),
  NumGlobalElementsHe(NumGlobalElementsHe_),
  NumGlobalElementsClad(NumGlobalElementsClad_),
  xminUO2(xminUO2_),
  xmaxUO2(xmaxUO2_),
  xminHe(xminHe_),
  xmaxHe(xmaxHe_),
  xminClad(xminClad_),
  xmaxClad(xmaxClad_),
  dt(1.0e+10),
  restart(restart_),
  Comm(&Comm_)
{

  NumGlobalNodes = NumGlobalElementsUO2 + NumGlobalElementsHe + NumGlobalElementsClad + 1;

  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Here we assume a 2-species PelletTransport model, ie 2 dofs per node
  // Note that this needs to be echoed in thew anonymous enum for NUMSPECIES.
  NumSpecies = 2;
  NumGlobalUnknowns = NumSpecies * NumGlobalNodes;

  // Construct a Source Map that puts approximately the same
  // number of equations on each processor

  // Begin by distributing nodes fairly equally
  StandardNodeMap = new Epetra_Map(NumGlobalNodes, 0, *Comm);

  // Get the number of nodes owned by this processor
  NumMyNodes = StandardNodeMap->NumMyElements();

  // Construct an overlap node map for the finite element fill
  // For single processor jobs, the overlap and standard maps are the same
  if (NumProc == 1)
    OverlapNodeMap = new Epetra_Map(*StandardNodeMap);
  else
  {
    // Distribute elements such that nodes are ghosted
    int OverlapNumMyNodes;
    int OverlapMinMyNodeGID;

    OverlapNumMyNodes = NumMyNodes + 1;
    if ( MyPID == NumProc - 1 )
      OverlapNumMyNodes --;

    OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();

    int* OverlapMyGlobalNodes = new int[OverlapNumMyNodes];

    for( int i = 0; i < OverlapNumMyNodes; ++i )
      OverlapMyGlobalNodes[i] = OverlapMinMyNodeGID + i;

    OverlapNodeMap = new Epetra_Map(-1, OverlapNumMyNodes,
                          OverlapMyGlobalNodes, 0, *Comm);

    delete [] OverlapMyGlobalNodes;

  } // End Overlap node map construction ********************************

  // Now create the unknowns maps from the node maps
  NumMyUnknowns = NumSpecies * NumMyNodes;
  int* StandardMyGlobalUnknowns = new int[NumMyUnknowns];
  for( int k = 0; k < NumSpecies; ++k )
    for( int i = 0; i < NumMyNodes; ++i )
    {
      // For now, we employ an interleave of unknowns
      StandardMyGlobalUnknowns[ NumSpecies * i + k ] = NumSpecies * StandardNodeMap->GID(i) + k;
    }

  StandardMap = new Epetra_Map(-1, NumMyUnknowns, StandardMyGlobalUnknowns, 0, *Comm);
  delete [] StandardMyGlobalUnknowns;

  assert(StandardMap->NumGlobalElements() == NumGlobalUnknowns);

  if (NumProc == 1)
    OverlapMap = new Epetra_Map(*StandardMap);
  else
  {
    int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();
    int OverlapNumMyUnknowns = NumSpecies * OverlapNumMyNodes;
    int* OverlapMyGlobalUnknowns = new int[OverlapNumMyUnknowns];
    for (int k=0; k<NumSpecies; k++)
      for( int i = 0; i < OverlapNumMyNodes; ++i )
        OverlapMyGlobalUnknowns[ NumSpecies * i + k ] = NumSpecies * OverlapNodeMap->GID(i) + k;

    OverlapMap = new Epetra_Map(-1, OverlapNumMyUnknowns, OverlapMyGlobalUnknowns, 0, *Comm);
    delete [] OverlapMyGlobalUnknowns;

  } // End Overlap unknowns map construction ***************************


#ifdef DEBUG_MAPS
  // Output to check progress so far
  printf("[%d] NumSpecies, NumMyNodes, NumGlobalNodes, NumMyUnknowns, NumGlobalUnknowns :\n %d\t%d\t%d\t%d\n",
    MyPID,NumSpecies, NumMyNodes, NumGlobalNodes, NumMyUnknowns, NumGlobalUnknowns);

  Comm->Barrier();
  Comm->Barrier();
  Comm->Barrier();
  printf("[%d] StandardNodeMap :\n", MyPID);
  Comm->Barrier();
  std::cout << *StandardNodeMap << std::endl;
  Comm->Barrier();
  printf("[%d] OverlapNodeMap :\n", MyPID);
  Comm->Barrier();
  std::cout << *OverlapNodeMap << std::endl;
  Comm->Barrier();
  printf("[%d] StandardMap :\n", MyPID);
  Comm->Barrier();
  std::cout << *StandardMap << std::endl;
  Comm->Barrier();
  printf("[%d] StandardMap :\n", MyPID);
  Comm->Barrier();
  std::cout << *OverlapMap << std::endl;
  Comm->Barrier();
#endif

  // Construct Linear Objects
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  nodeImporter = new Epetra_Import(*OverlapNodeMap, *StandardNodeMap);
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  oldSolution = new Epetra_Vector(*StandardMap);
  AA = Teuchos::rcp(new Epetra_CrsGraph(Copy, *StandardMap, 0));

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraphUsingNodes(*AA);

#ifdef DEBUG_GRAPH
  AA->Print(cout);
#endif
#ifdef DEBUG_IMPORTER
  Importer->Print(cout);
#endif

  // Create the element type int vector
  int numLocalElems = OverlapNodeMap->NumMyElements() - 1;
  Epetra_Map * elemMap = new Epetra_Map(-1, numLocalElems, 0, *Comm);
  elemTypes = Teuchos::rcp(new Epetra_IntVector(*elemMap));
  for( int i = 0; i < numLocalElems; ++i )
    if( i < NumGlobalElementsUO2 )
      (*elemTypes)[i] = UO2;
    else if( i < (NumGlobalElementsUO2 + NumGlobalElementsHe) )
      (*elemTypes)[i] = HELIUM;
      //(*elemTypes)[i] = UO2;
    else
      (*elemTypes)[i] = CLAD;
      //(*elemTypes)[i] = UO2;
  delete elemMap;

  // Create the nodal coordinates - only works in serial for now - RWH
  xptr = Teuchos::rcp(new Epetra_Vector(*StandardNodeMap));
  int inode = 0;
  double length = xmaxUO2 - xminUO2;
  double dx     = length/((double) NumGlobalElementsUO2 );
  (*xptr)[inode++] = xminUO2;
  for( int i = 1; i < NumGlobalElementsUO2; ++i )
  {
    (*xptr)[inode] = (*xptr)[inode-1] + dx;
    inode++;
  }
  length = xmaxHe - xminHe;
  dx     = length/((double) NumGlobalElementsHe );
  (*xptr)[inode++] = xminHe;
  for( int i = 1; i < NumGlobalElementsHe; ++i )
  {
    (*xptr)[inode] = (*xptr)[inode-1] + dx;
    inode++;
  }
  length = xmaxClad - xminClad;
  dx     = length/((double) NumGlobalElementsClad );
  (*xptr)[inode++] = xminClad;
  for( int i = 0; i < NumGlobalElementsClad; ++i )
  {
    (*xptr)[inode] = (*xptr)[inode-1] + dx;
    inode++;
  }

  initializeSoln();

}

// Destructor
PelletTransport::~PelletTransport()
{
  delete oldSolution;
  delete nodeImporter;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
  delete StandardNodeMap;
  delete OverlapNodeMap;
}

// Reset function
void
PelletTransport::reset(const Epetra_Vector& x)
{
  *oldSolution = x;
}

// Set initialSolution to desired initial condition
void
PelletTransport::initializeSoln()
{
  Epetra_Vector& soln = *initialSolution;
  Epetra_Vector& x = *xptr;

  // Here we do a sinusoidal perturbation of the unstable
  // steady state.

  if( restart )
  {
    std::string name = "restartVec";
    Epetra_Vector * tmpVec = NULL;
    NOX::Epetra::DebugTools::readVector(name, soln.Comm(), tmpVec );
    soln = *tmpVec;
  }
  else
  {
    for( int i = 0; i < x.MyLength(); ++i )
    {
      soln[2*i]   = 750.0 ;
      soln[2*i+1] = 0.02  ;
    }
  }
  *oldSolution = soln;
}

bool
PelletTransport::evaluate(NOX::Epetra::Interface::Required::FillType fType,
                    const Epetra_Vector* soln,
                    Epetra_Vector* tmp_rhs,
                    Epetra_RowMatrix* tmp_matrix)
{

  if( fType == NOX::Epetra::Interface::Required::Jac )
  {
    std::cout << "This problem only works for Finite-Difference or Matrix-Free Based Jacobians."
         << " No analytic Jacobian fill available !!" << std::endl;
    throw "Problem ERROR";
  }
  else
    rhs = new Epetra_Vector(*OverlapMap);

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  Epetra_Vector xvec(*OverlapNodeMap);

  // Export Solution to Overlap vector
  uold.Import(*oldSolution, *Importer, Insert);
  xvec.Import(*xptr, *nodeImporter, Insert);
  u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i;
  int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();

  MaterialPropBase::PropData materialProps;
  // Hard-code use of UO2 for now - RWH 5/14/2007
  std::string fixedType = "UO2";

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardNodeMap->MinMyGID()-1;

  int row1, row2;
  double term1, term2;
  double flux1, flux2;
  double xx[2];
  double uu[2*NUMSPECIES]; // Use of the anonymous enum is needed for SGI builds
  double uuold[2*NUMSPECIES];
  Basis basis(NumSpecies);


  // Zero out the objects that will be filled
  rhs->PutScalar(0.0);

  ACTIVE_REGIONS propType;

  // Loop Over # of Finite Elements on Processor
  for( int ne = 0; ne < OverlapNumMyNodes - 1; ++ne )
  {
    propType = (ACTIVE_REGIONS) (*elemTypes)[ne];

    // Loop Over Gauss Points
    for( int gp = 0; gp < 2; ++gp )
    {
      // Get the solution and coordinates at the nodes
      xx[0] = xvec[ne];
      xx[1] = xvec[ne+1];
      for (int k=0; k<NumSpecies; k++) {
        uu[NumSpecies * 0 + k] = u[NumSpecies * ne + k];
        uu[NumSpecies * 1 + k] = u[NumSpecies * (ne+1) + k];
        uuold[NumSpecies * 0 + k] = uold[NumSpecies * ne + k];
        uuold[NumSpecies * 1 + k] = uold[NumSpecies * (ne+1) + k];
      }
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uuold);
      MaterialPropFactory::computeProps( propType, basis.uu[0], basis.uu[1], materialProps );
      double & rho1    = materialProps.density  ;
      double & k1      = materialProps.k_thermal;
      double & Cp1     = materialProps.Cp       ;
      double & Qstar1  = materialProps.Qstar    ;
      double & Qdot1   = materialProps.Qdot     ;
      double & Ffunc   = materialProps.thermoF  ;
      double & D_diff1 = materialProps.D_diff   ;


      // Loop over Nodes in Element
      for( i = 0; i < 2; ++i )
      {
    row1=OverlapMap->GID(NumSpecies * (ne+i));
    row2=OverlapMap->GID(NumSpecies * (ne+i) + 1);
        flux1 = -k1*basis.duu[0]/basis.dx;
        flux2 = -0.5*D_diff1*(basis.duu[1]/basis.dx
              + basis.uu[1]*Qstar1/(Ffunc*8.3142*basis.uu[0]*basis.uu[0])*basis.duu[0]/basis.dx);
        term1 = basis.wt*basis.dx*basis.xx*(
                    rho1*Cp1*(basis.uu[0] - basis.uuold[0])/dt * basis.phi[i]
                  - flux1*basis.dphide[i]/basis.dx
                  - Qdot1*basis.phi[i]
                );
        term2 = basis.wt*basis.dx*basis.xx*(
                    0.5*(basis.uu[1] - basis.uuold[1])/dt * basis.phi[i]
                  - flux2*basis.dphide[i]/basis.dx
                );
        (*rhs)[NumSpecies*(ne+i)]   += term1;
        (*rhs)[NumSpecies*(ne+i)+1] += term2;
      }
    }
  }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // Dirichlet BCs at xminUO2
  const double xB = dynamic_cast<const MaterialProp_UO2 &>(MaterialPropFactory::factory().get_model(PelletTransport::UO2)).get_xB();
  if (MyPID==0)
  {
    // Use no-flux BCs
    //(*rhs)[0]= (*soln)[0] - 0.6;
    //(*rhs)[1]= (*soln)[1] - 10.0/3.0;
  }
  // Dirichlet BCs at xmaxUO2
  if( StandardMap->LID(NumSpecies*(NumGlobalElementsUO2)) >= 0 )
  {
    int lastUO2Dof = StandardMap->LID(NumSpecies*(NumGlobalElementsUO2) + 1);
    //(*rhs)[lastDof - 1] = (*soln)[lastDof - 1] - 840.0;
    (*rhs)[lastUO2Dof] = (*soln)[lastUO2Dof] - xB;
  }
  // Dirichlet BCs at all He and Clad species variables
  int lastUO2DofGID = NumSpecies*NumGlobalElementsUO2 + 1;
  int lastGID = StandardMap->MaxAllGID();
  for( int i = lastUO2DofGID; i < lastGID; i+=2 )
    if( StandardMap->LID(i) >= 0 )
      (*rhs)[StandardMap->LID(i)] = (*soln)[StandardMap->LID(i)] - xB;

  // Dirichlet BCs at xmaxClad
  if( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 )
  {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    (*rhs)[lastDof - 1] = (*soln)[lastDof - 1] - 750.0;
    (*rhs)[lastDof] = (*soln)[lastDof] - xB;
  }

  // Sync up processors to be safe
  Comm->Barrier();

  // Do an assemble for overlap nodes
  tmp_rhs->Export(*rhs, *Importer, Add);

  delete rhs;

  return true;
}

Teuchos::RCP<Epetra_Vector>
PelletTransport::getSolution()
{
  return initialSolution;
}

Teuchos::RCP<Epetra_CrsMatrix>
PelletTransport::getJacobian()
{
  if( Teuchos::is_null(A) ) return A;
  else {
    std::cout << "No valid Jacobian matrix for this problem.  This is likely the"
         << " result of overlapping NODES rather than ELEMENTS.\n" << std::endl;
    throw "PelletTransport Error";
  }
}

Teuchos::RCP<Epetra_Vector>
PelletTransport::getMesh()
{
  return xptr;
}

Epetra_Vector&
PelletTransport::getOldSoln()
{
  return *oldSolution;
}

double
PelletTransport::getdt()
{
  return dt;
}

Teuchos::RCP<Epetra_CrsGraph>
PelletTransport::getGraph()
{
  return AA;
}


Epetra_CrsGraph&
PelletTransport::generateGraphUsingNodes(Epetra_CrsGraph& AA)
{

  int row, column;

  int myMinNodeGID = StandardNodeMap->MinMyGID();
  int myMaxNodeGID = StandardNodeMap->MaxMyGID();

  int leftNodeGID, rightNodeGID;
  for( int myNode = myMinNodeGID; myNode <= myMaxNodeGID; myNode++ ) {

    leftNodeGID  = myNode - 1;
    rightNodeGID = myNode + 1;

    if( leftNodeGID < StandardNodeMap->MinAllGID() )
      leftNodeGID = StandardNodeMap->MinAllGID();

    if( rightNodeGID > StandardNodeMap->MaxAllGID() )
      rightNodeGID = StandardNodeMap->MaxAllGID();

    for( int dependNode = leftNodeGID; dependNode <= rightNodeGID; dependNode++ ) {

      // Loop over unknowns in Node
      for (int j = 0; j < NumSpecies; j++) {
        row = NumSpecies * myNode + j; // Interleave scheme

        // Loop over unknowns at supporting nodes
        for (int m = 0; m < NumSpecies; m++) {
          column = NumSpecies * dependNode + m;
          //printf("\t\tWould like to insert -> (%d, %d)\n",row,column);
          AA.InsertGlobalIndices(row, 1, &column);
        }
      }
    }
  }

  AA.FillComplete();

  return AA;
}

