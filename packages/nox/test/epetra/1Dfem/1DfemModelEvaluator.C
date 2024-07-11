// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// NOX include (for iostream, cmath, etc...)
#include "NOX_Common.H"

#include "Teuchos_ParameterList.hpp"

// Class Definition
#include "1DfemModelEvaluator.H"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

ModelEvaluatorInterface::ModelEvaluatorInterface(int numGlobalElements,
                         Epetra_Comm& comm,
                         double xmin_,
                         double xmax_) :
  NumGlobalElements(numGlobalElements),
  NumMyElements(0),  // gets set after map creation
  MyPID(comm.MyPID()),
  NumProc(comm.NumProc()),
  xmin(xmin_),
  xmax(xmax_),
  factor(1.0),
  Comm(&comm),
  Importer(0),
  rhs(0),
  Graph(0)
{

  // Construct a Source Map that puts approximately the same
  // Number of equations on each processor in uniform global ordering
  StandardMap = Teuchos::rcp(new Epetra_Map(NumGlobalElements, 0, *Comm));

  // Get the number of elements owned by this processor
  NumMyElements = StandardMap->NumMyElements();

  // Construct an overlaped map for the finite element fill *************
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) {
    OverlapMap = Teuchos::rcp(new Epetra_Map(*StandardMap));
  } else {

    int OverlapNumMyElements;
    int OverlapMinMyGID;
    OverlapNumMyElements = NumMyElements + 2;
    if ((MyPID == 0) || (MyPID == NumProc - 1))
      OverlapNumMyElements --;

    if (MyPID==0)
      OverlapMinMyGID = StandardMap->MinMyGID();
    else
      OverlapMinMyGID = StandardMap->MinMyGID() - 1;

    int* OverlapMyGlobalElements = new int[OverlapNumMyElements];

    for (int i = 0; i < OverlapNumMyElements; i ++)
      OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

    OverlapMap = Teuchos::rcp(new Epetra_Map(-1, OverlapNumMyElements,
                         OverlapMyGlobalElements, 0,
                         *Comm));

    delete [] OverlapMyGlobalElements;

  } // End Overlap map construction *************************************

  // Construct Linear Objects
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));

  // Assign non-zero entries in the graph
  createGraph();

  // Construct a matrix
  jacobian = Teuchos::rcp(new Epetra_CrsMatrix (Copy, *Graph));

  // Clean-up
  jacobian->FillComplete();

  // Create the nodal coordinates
  xptr = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  double Length = xmax - xmin;
  double dx = Length/((double) NumGlobalElements-1);
  for (int i=0; i < NumMyElements; i++) {
    (*xptr)[i] = xmin + dx*((double) StandardMap->MinMyGID()+i);
  }

  initializeSoln();

}

ModelEvaluatorInterface::~ModelEvaluatorInterface()
{
  delete Graph;
  delete Importer;
}

void ModelEvaluatorInterface::evalModel(const InArgs& inArgs,
                    const OutArgs& outArgs) const
{
  Epetra_RowMatrix* tmp = dynamic_cast<Epetra_RowMatrix*>
    (outArgs.get_W().get());

  this->evaluate( inArgs.get_x().get(), outArgs.get_f().get(), tmp );

  return;
}

bool ModelEvaluatorInterface::evaluate(const Epetra_Vector* soln,
                       Epetra_Vector* tmp_rhs,
                       Epetra_RowMatrix* tmp_matrix) const
{
  //Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  if (tmp_rhs != 0) {
    fillF = true;
    rhs = tmp_rhs;
  }
  if (tmp_matrix != 0) {
    fillMatrix = true;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector x(*OverlapMap);

  // Export Solution to Overlap vector
  u.Import(*soln, *Importer, Insert);
  x.Import(*xptr, *Importer, Insert);

  // Declare required variables
  int OverlapNumMyElements = OverlapMap->NumMyElements();

  int row, column;
  double jac;
  double xx[2];
  double uu[2];
  Basis basis;

  // Zero out the objects that will be filled
  if (fillF)
    rhs->PutScalar(0.0);
  if (fillMatrix)
    jacobian->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {

    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes
      xx[0]=x[ne];
      xx[1]=x[ne+1];
      uu[0]=u[ne];
      uu[1]=u[ne+1];
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu);

      // Loop over Nodes in Element
      for (int i=0; i< 2; i++) {
    row=OverlapMap->GID(ne+i);
    //printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
    //     MyPID, row, ne+i,StandardMap.MyGID(row));
    if (StandardMap->MyGID(row)) {
      if (fillF) {
        (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
          +basis.wt*basis.dx
          *((1.0/(basis.dx*basis.dx))*basis.duu*
        basis.dphide[i]+factor*basis.uu*basis.uu*basis.phi[i]);
      }
    }
    // Loop over Trial Functions
    if (fillMatrix) {
      for(int j=0;j < 2; j++) {
        if (StandardMap->MyGID(row)) {
          column=OverlapMap->GID(ne+j);
          jac=basis.wt*basis.dx*((1.0/(basis.dx*basis.dx))*
                     basis.dphide[j]*basis.dphide[i]
                     +2.0*factor*basis.uu*basis.phi[j]*
                     basis.phi[i]);
          jacobian->SumIntoGlobalValues(row, 1, &jac, &column);
        }
      }
    }
      }
    }
  }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(0)=1
  if (MyPID==0) {
    if (fillF)
      (*rhs)[0]= (*soln)[0] - 1.0;
    if (fillMatrix) {
      column=0;
      jac=1.0;
      jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();

  jacobian->FillComplete();

  return true;
}

Teuchos::RCP<const Epetra_Map> ModelEvaluatorInterface::
get_x_map() const
{
  return StandardMap;
}

Teuchos::RCP<const Epetra_Map> ModelEvaluatorInterface::
get_f_map() const
{
  return StandardMap;
}

Teuchos::RCP<const Epetra_Vector> ModelEvaluatorInterface::
get_x_init() const
{
  return initialSolution;
}

Teuchos::RCP<Epetra_Operator> ModelEvaluatorInterface::
create_W() const
{
  return jacobian;
}

EpetraExt::ModelEvaluator::InArgs ModelEvaluatorInterface::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription("NOX 1D Finite Element Test Problem");
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs ModelEvaluatorInterface::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("NOX 1D Finite Element Test Problem");
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,false // supportsAdjoint
      )
    );
  return outArgs;
}

bool ModelEvaluatorInterface::createGraph()
{
  if (Graph != 0) {
    delete Graph;
    Graph = 0;
  }

  // Create the shell for the
  Graph = new Epetra_CrsGraph(Copy, *StandardMap, 5);

  // Declare required variables
  int row, column;
  int OverlapNumMyElements = OverlapMap->NumMyElements();

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {

    // Loop over Nodes in Element
    for (int i=0; i< 2; i++) {
      row=OverlapMap->GID(ne+i);

      // Loop over Trial Functions
      for(int j=0;j < 2; j++) {

    // If this row is owned by current processor, add the index
    if (StandardMap->MyGID(row)) {
      column=OverlapMap->GID(ne+j);
      Graph->InsertGlobalIndices(row, 1, &column);
    }
      }
    }
  }
  Graph->FillComplete();
  return true;
}

// Set initialSolution to desired initial condition
bool ModelEvaluatorInterface::initializeSoln()
{
  initialSolution->PutScalar(1.0); // Default initialization
  return true;
}

//====================================================================
// Basis vector

// Constructor
Basis::Basis() {
  phi = new double[2];
  dphide = new double[2];
}

// Destructor
Basis::~Basis() {
  delete [] phi;
  delete [] dphide;
}

// Calculates a linear 1D basis
void Basis::computeBasis(int gp, double *x, double *u, double *uold) {
  int N = 2;
  if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;

  // Caculate basis function and derivative at GP.
  dx=0.5*(x[1]-x[0]);
  xx=0.0;
  uu=0.0;
  duu=0.0;
  uuold=0.0;
  duuold=0.0;
  for (int i=0; i < N; i++) {
    xx += x[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
    if (uold) {
      uuold += uold[i] * phi[i];
      duuold += uold[i] * dphide[i];
    }
  }

  return;
}
