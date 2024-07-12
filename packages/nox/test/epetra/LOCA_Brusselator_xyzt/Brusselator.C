// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
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

#include "Brusselator.H"

// Constructor - creates the Epetra objects (maps and vectors)
Brusselator::Brusselator(int numGlobalNodes, Epetra_Comm& comm,
                         OverlapType OType_) :
  xmin(0.0),
  xmax(1.0),
  Comm(&comm),
  overlapType(OType_),
  NumGlobalNodes(numGlobalNodes),
  ColumnToOverlapImporter(0),
  AA(0),
  A(0),
  alpha(0.6),
  beta(2.0)
{

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Here we assume a 2-species Brusselator model, ie 2 dofs per node
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
  else {

    if( overlapType == ELEMENTS ) {

      int OverlapNumMyNodes;
      int OverlapMinMyNodeGID;

      OverlapNumMyNodes = NumMyNodes + 2;
      if ((MyPID == 0) || (MyPID == NumProc - 1))
        OverlapNumMyNodes --;

      if (MyPID==0)
        OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();
      else
        OverlapMinMyNodeGID = StandardNodeMap->MinMyGID() - 1;

      int* OverlapMyGlobalNodes = new int[OverlapNumMyNodes];

      for (i = 0; i < OverlapNumMyNodes; i ++)
        OverlapMyGlobalNodes[i] = OverlapMinMyNodeGID + i;

      OverlapNodeMap = new Epetra_Map(-1, OverlapNumMyNodes,
                                OverlapMyGlobalNodes, 0, *Comm);

      delete [] OverlapMyGlobalNodes;
    }

    else { // overlapType == NODES

      // Here we distribute elements such that nodes are ghosted, i.e. no
      // overlapping elements but instead overlapping nodes
      int OverlapNumMyNodes;
      int OverlapMinMyNodeGID;

      OverlapNumMyNodes = NumMyNodes + 1;
      if ( MyPID == NumProc - 1 )
        OverlapNumMyNodes --;

      OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();

      int* OverlapMyGlobalNodes = new int[OverlapNumMyNodes];

        for (i = 0; i < OverlapNumMyNodes; i ++)
        OverlapMyGlobalNodes[i] = OverlapMinMyNodeGID + i;

      OverlapNodeMap = new Epetra_Map(-1, OverlapNumMyNodes,
                  OverlapMyGlobalNodes, 0, *Comm);

      delete [] OverlapMyGlobalNodes;
    }
  } // End Overlap node map construction ********************************

  // Now create the unknowns maps from the node maps
  NumMyUnknowns = NumSpecies * NumMyNodes;
  int* StandardMyGlobalUnknowns = new int[NumMyUnknowns];
  for (int k=0; k<NumSpecies; k++)
    for (i=0; i<NumMyNodes; i++)

      // For now, we employ an interleave of unknowns

      StandardMyGlobalUnknowns[ NumSpecies * i + k ] =
               NumSpecies * StandardNodeMap->GID(i) + k;

  StandardMap = new Epetra_Map(-1, NumMyUnknowns, StandardMyGlobalUnknowns,
                                 0, *Comm);
  delete [] StandardMyGlobalUnknowns;

  assert(StandardMap->NumGlobalElements() == NumGlobalUnknowns);

  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  }
  else {
    int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();
    int OverlapNumMyUnknowns = NumSpecies * OverlapNumMyNodes;
    int* OverlapMyGlobalUnknowns = new int[OverlapNumMyUnknowns];
    for (int k=0; k<NumSpecies; k++)
      for (i=0; i<OverlapNumMyNodes; i++)
        OverlapMyGlobalUnknowns[ NumSpecies * i + k ] =
                 NumSpecies * OverlapNodeMap->GID(i) + k;

    OverlapMap = new Epetra_Map(-1, OverlapNumMyUnknowns,
                                OverlapMyGlobalUnknowns, 0, *Comm);
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
  initialSolution = new Epetra_Vector(*StandardMap);
  solnDot = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 0);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  if( overlapType == NODES )
    generateGraphUsingNodes(*AA);
  else
    generateGraphUsingElements(*AA);

#ifdef DEBUG_GRAPH
  AA->Print(cout);
#endif
#ifdef DEBUG_IMPORTER
  Importer->Print(cout);
#endif

  // Create the Importer needed for FD coloring using element overlap
  // as well as the problem Jacobian matrix which may or may not be used
  // depending on the choice of Jacobian operator
  if( overlapType == ELEMENTS ) {
    ColumnToOverlapImporter = new Epetra_Import(AA->ColMap(),*OverlapMap);
    A = new Epetra_CrsMatrix (Copy, *AA);
    A->FillComplete();
  }

  // Create the nodal coordinates
  xptr = new Epetra_Vector(*StandardNodeMap);
  double Length= xmax - xmin;
  dx=Length/((double) NumGlobalNodes-1);
  for (i=0; i < NumMyNodes; i++) {
    (*xptr)[i]=xmin + dx*((double) StandardNodeMap->MinMyGID()+i);
  }

  initializeSoln();

}

// Destructor
Brusselator::~Brusselator()
{
  delete AA;
  delete xptr;
  delete initialSolution;
  delete solnDot;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

// Reset function
void Brusselator::setxDot(const Epetra_Vector& xdot)
{
  *solnDot = xdot;
}

// Set initialSolution to desired initial condition
void Brusselator::initializeSoln()
{
  Epetra_Vector& soln = *initialSolution;
  Epetra_Vector& x = *xptr;

  // Here we do a sinusoidal perturbation of the unstable
  // steady state.

  for (int i=0; i<x.MyLength(); i++) {
    soln[2*i] = alpha;
    soln[2*i+1] = beta/alpha;
  }

/*
  double pi = 4.*atan(1.0);

  for (int i=0; i<x.MyLength(); i++) {
    soln[2*i] = 0.6 + 1.e-1*sin(1.0*pi*x[i]);
    soln[2*i+1] = 10.0/3.0 + 1.e-1*sin(1.0*pi*x[i]);
  }
*/
  solnDot->PutScalar(0.0);
}

void Brusselator::setParameters(const double alpha_, const double beta_)
{
  alpha = alpha_;
  beta  = beta_;
  //cout<<" Brusselator::setParameters: alpha = "<<alpha<<"  beta = "<<beta<<endl;
}

// Matrix and Residual Fills
bool Brusselator::evaluate(NOX::Epetra::Interface::Required::FillType fType,
                    const Epetra_Vector* soln,
                    Epetra_Vector* tmp_rhs,
                    Epetra_RowMatrix* tmp_matrix,
                                    double jac_coeff, double mass_coeff)
{

  if( fType == NOX::Epetra::Interface::Required::Jac ) {
    if( overlapType == NODES ) {
      std::cout << "This problem only works for Finite-Difference Based Jacobians"
           << std::endl << "when overlapping nodes." << std::endl
           << "No analytic Jacobian fill available !!" << std::endl;
      exit(1);
    }
    flag = MATRIX_ONLY;
  }
  else {
    flag = F_ONLY;
    if( overlapType == NODES )
      rhs = new Epetra_Vector(*OverlapMap);
    else
      rhs = tmp_rhs;
  }

//  std::cout << "\nfTpye--> " << fType << std::endl
//       << "Incoming soln vector :\n" << std::endl << *soln << std::endl;

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector udot(*OverlapMap);
  Epetra_Vector xvec(*OverlapNodeMap);

  // Export Solution to Overlap vector
  // If the vector to be used in the fill is already in the Overlap form,
  // we simply need to map on-processor from column-space indices to
  // OverlapMap indices. Note that the old solution is simply fixed data that
  // needs to be sent to an OverlapMap (ghosted) vector.  The conditional
  // treatment for the current soution vector arises from use of
  // FD coloring in parallel.
  udot.Import(*solnDot, *Importer, Insert);
  xvec.Import(*xptr, *nodeImporter, Insert);
  if( (overlapType == ELEMENTS) &&
      (fType == NOX::Epetra::Interface::Required::FD_Res) )
    // Overlap vector for solution received from FD coloring, so simply reorder
    // on processor
    u.Export(*soln, *ColumnToOverlapImporter, Insert);
  else // Communication to Overlap vector is needed
    u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i,j;
  int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();

  int row1, row2, column1, column2;
  double term1, term2;
  double Dcoeff1 = 0.025;
  double Dcoeff2 = 0.025;
//  double alpha = 0.6;
//  double beta = 2.0;
  double jac11, jac12, jac21, jac22;
  double xx[2];
  double uu[2*NUMSPECIES] = {}; // Use of the anonymous enum is needed for SGI builds
  double uudot[2*NUMSPECIES] = {};
  Basis basis(NumSpecies);


  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) A->PutScalar(0.0);
  if ((flag == F_ONLY)      || (flag == ALL)) rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes - 1; ne++) {

    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      for (int k=0; k<NumSpecies; k++) {
        uu[NumSpecies * 0 + k] = u[NumSpecies * ne + k];
        uu[NumSpecies * 1 + k] = u[NumSpecies * (ne+1) + k];
        uudot[NumSpecies * 0 + k] = udot[NumSpecies * ne + k];
        uudot[NumSpecies * 1 + k] = udot[NumSpecies * (ne+1) + k];
      }
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uudot);

      // Loop over Nodes in Element
      for (i=0; i< 2; i++) {
    row1=OverlapMap->GID(NumSpecies * (ne+i));
    row2=OverlapMap->GID(NumSpecies * (ne+i) + 1);
        term1 = basis.wt*basis.dx * (
             basis.uudot[0] * basis.phi[i]
              +(1.0/(basis.dx*basis.dx))*Dcoeff1*basis.duu[0]*basis.dphide[i]
              + basis.phi[i] * ( -alpha + (beta+1.0)*basis.uu[0]
              - basis.uu[0]*basis.uu[0]*basis.uu[1]) );
        term2 = basis.wt*basis.dx * (
              basis.uudot[1] * basis.phi[i]
              +(1.0/(basis.dx*basis.dx))*Dcoeff2*basis.duu[1]*basis.dphide[i]
              + basis.phi[i] * ( -beta*basis.uu[0]
              + basis.uu[0]*basis.uu[0]*basis.uu[1]) );
    //printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
    //     MyPID, row, ne+i,StandardMap.MyGID(row));
        if ((flag == F_ONLY)    || (flag == ALL)) {
          if( overlapType == NODES ) {
            (*rhs)[NumSpecies*(ne+i)]   += term1;
        (*rhs)[NumSpecies*(ne+i)+1] += term2;
          }
      else
        if (StandardMap->MyGID(row1)) {
              (*rhs)[StandardMap->LID(OverlapMap->GID(NumSpecies*(ne+i)))]+=
                term1;
              (*rhs)[StandardMap->LID(OverlapMap->GID(NumSpecies*(ne+i)+1))]+=
                term2;
      }
    }
        // Loop over Trial Functions
        if ((flag == MATRIX_ONLY) || (flag == ALL)) {
          for(j=0;j < 2; j++) {
            if (StandardMap->MyGID(row1)) {
              column1=OverlapMap->GID(NumSpecies * (ne+j));
              column2=OverlapMap->GID(NumSpecies * (ne+j) + 1);
              jac11=basis.wt*basis.dx*(
                      mass_coeff*basis.phi[j]*basis.phi[i]
                      +jac_coeff*
                      +((1.0/(basis.dx*basis.dx))*Dcoeff1*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * ( (beta+1.0)*basis.phi[j]
                      - 2.0*basis.uu[0]*basis.phi[j]*basis.uu[1])) );
              jac12=basis.wt*basis.dx*(
                      jac_coeff*basis.phi[i] * ( -basis.uu[0]*basis.uu[0]*basis.phi[j]) );
              jac21=basis.wt*basis.dx*(
                      jac_coeff*basis.phi[i] * ( -beta*basis.phi[j]
                      + 2.0*basis.uu[0]*basis.phi[j]*basis.uu[1]) );
              jac22=basis.wt*basis.dx*(
                      mass_coeff*basis.phi[j]*basis.phi[i]
                      +jac_coeff*
                      +((1.0/(basis.dx*basis.dx))*Dcoeff2*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * basis.uu[0]*basis.uu[0]*basis.phi[j] ));
              A->SumIntoGlobalValues(row1, 1, &jac11, &column1);
              column1++;
              A->SumIntoGlobalValues(row1, 1, &jac12, &column1);
              A->SumIntoGlobalValues(row2, 1, &jac22, &column2);
              column2--;
              A->SumIntoGlobalValues(row2, 1, &jac21, &column2);
            }
          }
        }
      }
    }
  }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(0)=1
  if (MyPID==0) {
    if ((flag == F_ONLY)    || (flag == ALL)) {
      (*rhs)[0]= (*soln)[0] - alpha;
      (*rhs)[1]= (*soln)[1] - beta/alpha;
    }
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int column=0;
      double jac=1.0*jac_coeff;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      jac=0.0;
      column=0;
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      A->ReplaceGlobalValues(1, 1, &jac, &column);
    }
  }
  // U(1)=1
  if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 ) {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    if ((flag == F_ONLY)    || (flag == ALL)) {
      (*rhs)[lastDof - 1] = (*soln)[lastDof - 1] - alpha;
      (*rhs)[lastDof] = (*soln)[lastDof] - beta/alpha;
    }
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int row=StandardMap->MaxAllGID() - 1;
      int column = row;
      double jac = 1.0*jac_coeff;
      A->ReplaceGlobalValues(row++, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      jac=0.0;
      row = column - 1;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();

  // Do an assemble for overlap nodes
  if( overlapType == NODES )
    tmp_rhs->Export(*rhs, *Importer, Add);

//  Comm->Barrier();
//  std::cout << "Returned tmp_rhs residual vector :\n" << std::endl << *tmp_rhs << std::endl;

  return true;
}

Epetra_Vector& Brusselator::getSolution()
{
  return *initialSolution;
}

Epetra_CrsMatrix& Brusselator::getJacobian()
{
  if( A ) return *A;
  else {
    std::cout << "No valid Jacobian matrix for this problem.  This is likely the"
         << " result of overlapping NODES rather than ELEMENTS.\n" << std::endl;
    throw "Brusselator Error";
  }
}

Epetra_Vector& Brusselator::getMesh()
{
  return *xptr;
}

Epetra_Vector& Brusselator::getSolnDot()
{
  return *solnDot;
}

Epetra_CrsGraph& Brusselator::getGraph()
{
  return *AA;
}


Epetra_CrsGraph& Brusselator::generateGraphUsingNodes(Epetra_CrsGraph& AA)
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
//   AA.SortIndices();
//   AA.RemoveRedundantIndices();
  return AA;
}

Epetra_CrsGraph& Brusselator::generateGraphUsingElements(Epetra_CrsGraph& AA)
{

  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {

    // Loop over Nodes in Element
    for (i=0; i<2; i++) {

      // If this node is owned by current processor, add indices
      if (StandardNodeMap->MyGID(OverlapNodeMap->GID(ne+i))) {

        // Loop over unknowns in Node
        for (int k=0; k<NumSpecies; k++) {
          row=OverlapMap->GID( NumSpecies*(ne+i) + k); // Interleave scheme

          // Loop over supporting nodes
          for(j=0; j<2; j++) {

            // Loop over unknowns at supporting nodes
            for (int m=0; m<NumSpecies; m++) {
              column=OverlapMap->GID( NumSpecies*(ne+j) + m);
              //printf("\t\tWould like to insert -> (%d, %d)\n",row,column);
              AA.InsertGlobalIndices(row, 1, &column);
            }
          }
        }
      }
    }
  }
  AA.FillComplete();
//   AA.SortIndices();
//   AA.RemoveRedundantIndices();
  return AA;
}

