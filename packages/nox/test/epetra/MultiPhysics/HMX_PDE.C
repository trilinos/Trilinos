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

#include "Problem_Manager.H"
#include "HMX_PDE.H"

// Constructor - creates the Epetra objects (maps and vectors)
HMX_PDE::HMX_PDE(Epetra_Comm& comm,
          double diffCoef_,
          double Const_R_,
          double Steric_C_,
      double PreExp_A_,
      double ActEnergy_,
      map<string, double> SrcTermExponent_,
      map<string, double> SrcTermWeight_,
          int numGlobalNodes,
          std::string name_) :
  GenericEpetraProblem(comm, numGlobalNodes, name_),
  xmin(0.0),
  xmax(1.0),
  dt(2.0e-1),
  diffCoef(diffCoef_),
  Const_R(Const_R_),
  StericCoef(Steric_C_),
  PreExp_A(PreExp_A_),
  ActEnergy(ActEnergy_),
  SrcTermExponent(SrcTermExponent_),
  SrcTermWeight(SrcTermWeight_)
{
  // Create mesh and solution vectors

  // We first initialize the mesh and then the solution since the latter
  // can depend on the mesh.
  xptr = Teuchos::rcp( new Epetra_Vector(*StandardMap) );
  double Length= xmax - xmin;
  dx=Length/((double) NumGlobalNodes-1);
  for (int i=0; i < NumMyNodes; i++) {
    (*xptr)[i]=xmin + dx*((double) StandardMap->MinMyGID()+i);
  }

  // Create extra vector needed for this transient problem
  oldSolution = new Epetra_Vector(*StandardMap);

  // Next we create and initialize (using default provided) the solution vector
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
  ColumnToOverlapImporter = new Epetra_Import(A->ColMap(),*OverlapMap);
}

//-----------------------------------------------------------------------------

// Destructor
HMX_PDE::~HMX_PDE()
{
  delete oldSolution; oldSolution = 0;
  delete ColumnToOverlapImporter; ColumnToOverlapImporter = 0;
}

//-----------------------------------------------------------------------------

// Reset function
void
HMX_PDE::reset(const Epetra_Vector& x)
{
  *oldSolution = x;
}

//-----------------------------------------------------------------------------

// Empty Reset function
void
HMX_PDE::reset()
{
  std::cout << "WARNING: reset called without passing any update vector !!"
       << std::endl;
}

//-----------------------------------------------------------------------------

// Set initialSolution to desired initial condition
void
HMX_PDE::initializeSolution(double val)
{
  // Aliases for convenience
  Epetra_Vector& soln = *initialSolution;

  soln.PutScalar(val);

  *oldSolution = soln;
}

//-----------------------------------------------------------------------------

// Matrix and Residual Fills
bool
HMX_PDE::evaluate(
                    NOX::Epetra::Interface::Required::FillType flag,
            const Epetra_Vector* soln,
            Epetra_Vector* tmp_rhs)
{
  // Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  if (tmp_rhs != 0) {
    fillF = true;
    rhs = tmp_rhs;
  }
  else {
    fillMatrix = true;
  }

  // "flag" can be used to determine how accurate your fill of F should be
  // depending on why we are calling evaluate (Could be using computeF to
  // populate a Jacobian or Preconditioner).
  if (flag == NOX::Epetra::Interface::Required::Residual) {
    // Do nothing for now
  }
  else if (flag == NOX::Epetra::Interface::Required::Jac) {
    // Do nothing for now
  }
  else if (flag == NOX::Epetra::Interface::Required::Prec) {
    // Do nothing for now
  }
  else if (flag == NOX::Epetra::Interface::Required::User) {
    // Do nothing for now
  }

  int numDep = depProblems.size();

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  std::vector<Epetra_Vector*> dep(numDep);
  for( int i = 0; i<numDep; i++)
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
  for( int i = 0; i<numDep; i++ )
    (*dep[i]).Import(*( (*(depSolutions.find(depProblems[i]))).second ),
                   *Importer, Insert);
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

  // Setup iterators for looping over each problem source term contribution
  // to this one's PDE
  map<string, double>::iterator srcTermIter;
  map<string, double>::iterator srcTermEnd = SrcTermWeight.end();

  // Bundle up the dependent variables in the way needed for computing
  // the source terms of each reaction
  /*
  Epetra_Vector debugSrcTerm(*OverlapMap);
  map<string, Epetra_Vector*> debugDepVars;
  debugDepVars.insert( pair<string, Epetra_Vector*>(getName(), &u) );
  for( int i = 0; i<numDep; i++ )
    debugDepVars.insert( pair<string, Epetra_Vector*>
                    (myManager->getName(depProblems[i]), &dep[i]) );

  for( srcTermIter = SrcTermWeight.begin();
    srcTermIter != srcTermEnd; srcTermIter++) {
                     HMX_PDE &srcTermProb =
                     dynamic_cast<HMX_PDE&>(
                     myManager->getProblem(srcTermIter->first) );
    std::cout << "Inside problem: \"" << getName() << "\" calling to get source term "
         << "from problem: \"" << srcTermIter->first << "\" :" << std::endl;
    srcTermProb.computeSourceTerm(debugDepVars, debugSrcTerm);
    std::cout << "Resulting source term :" << debugSrcTerm << std::endl;
  }
  */

  int row;
  double alpha = 500.0;
  double xx[2];
  double uu[2];
  double uuold[2];
  std::vector<double*> ddep(numDep);
  for( int i = 0; i<numDep; i++)
    ddep[i] = new double[2];
  double *srcTerm = new double[2];
  Basis basis;

  // Bundle up the dependent variables in the way needed for computing
  // the source terms of each reaction
  map<string, double*> depVars;
  depVars.insert( pair<string, double*>(getName(), uu) );
  for( int i = 0; i<numDep; i++ )
    depVars.insert( pair<string, double*>
                    (myManager->getProblemName(depProblems[i]), ddep[i]) );
  // Do a check on this fill
//  map<string, double*>::iterator iter;
//  for( iter = depVars.begin(); iter != depVars.end(); iter++)
//    std::cout << "Inserted ... " << iter->first << "\t" << iter->second << std::endl;
//  std::cout << "--------------------------------------------------" << std::endl;
//  for( iter = depVars.begin(); iter != depVars.end(); iter++)
//      std::cout << iter->first << "\t" << (iter->second)[0] << ", "
//               << (iter->second)[1] << std::endl;
//  std::cout << "--------------------------------------------------" << std::endl;

  // Zero out the objects that will be filled
  if ( fillMatrix ) A->PutScalar(0.0);
  if ( fillF ) rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++)
  {
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++)
    {
      // Get the solution and coordinates at the nodes
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      uu[0] = u[ne];
      uu[1] = u[ne+1];
      uuold[0] = uold[ne];
      uuold[1] = uold[ne+1];
      for( int i = 0; i<numDep; i++ ) {
        ddep[i][0] = (*dep[i])[ne];
        ddep[i][1] = (*dep[i])[ne+1];
      }
      // Calculate the basis function and variables at the gauss points
      basis.getBasis(gp, xx, uu, uuold, ddep);

      // Loop over Nodes in Element
      for (int i=0; i< 2; i++) {
    row=OverlapMap->GID(ne+i);
    if (StandardMap->MyGID(row)) {
      if ( fillF ) {

            // First do time derivative and diffusion operator
        (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
          +basis.wt*basis.dx
          *((basis.uu - basis.uuold)/dt * basis.phi[i]
              +(1.0/(basis.dx*basis.dx))*diffCoef*basis.duu*basis.dphide[i]);

            // Then do source term contributions
        //
            for( srcTermIter = SrcTermWeight.begin();
                          srcTermIter != srcTermEnd; srcTermIter++) {
              HMX_PDE &srcTermProb =
                       dynamic_cast<HMX_PDE&>(
                   myManager->getProblem((*srcTermIter).first) );
              srcTermProb.computeSourceTerm(2, depVars, srcTerm);
              (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
                +basis.wt*basis.dx
                *( basis.phi[i] * ( - (*srcTermIter).second * srcTerm[i] ));
            }
        //
      }
    }
    // Loop over Trial Functions
    if ( fillMatrix ) {
        /*
      for(j=0;j < 2; j++) {
        if (StandardMap->MyGID(row)) {
          column=OverlapMap->GID(ne+j);
          jac=basis.wt*basis.dx*(
                      basis.phi[j]/dt*basis.phi[i]
                      +(1.0/(basis.dx*basis.dx))*diffCoef*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * ( (beta+1.0)*basis.phi[j]
                      - 2.0*basis.uu*basis.phi[j]*basis.ddep[id_spec]) );
          ierr=A->SumIntoGlobalValues(row, 1, &jac, &column);
        }
      }
   */
    }
      }
    }
  }

  // Apply Dirichlet BC for Temperature problem only (for now); this implies
  // no-flux across domain boundary for all species.
  if( getName() == tempFieldName )
  {
    // Insert Boundary Conditions and modify Jacobian and function (F)
    // U(0)=1
    if (MyPID==0)
    {
      if ( fillF )
        (*rhs)[0]= (*soln)[0] - alpha;
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
    // U(1)=1
    if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 )
    {
      int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
      if ( fillF )
        (*rhs)[lastDof] = (*soln)[lastDof] - alpha;
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
  }

  // Sync up processors to be safe
  Comm->Barrier();

  A->FillComplete();

#ifdef DEBUG
  A->Print(cout);

  if( fillF )
    std::cout << "For residual fill :" << std::endl << *rhs << std::endl;

  if( fillMatrix )
  {
    std::cout << "For jacobian fill :" << std::endl;
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

Epetra_Vector &
HMX_PDE::getOldSoln()
{
  return *oldSolution;
}

//-----------------------------------------------------------------------------

double
HMX_PDE::getdt() const
{
  return dt;
}

//-----------------------------------------------------------------------------

void
HMX_PDE::generateGraph()
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

void
HMX_PDE::computeSourceTerm(map<string, Epetra_Vector*> fields,
                                Epetra_Vector& result)
{
  // All dependent variables should be in place before now

  // We assume all other registered dependent problems are species which
  // affect the reaction rate of this specie

  Epetra_Vector *TvecPtr = 0;
  map<string, Epetra_Vector*>::iterator iter = fields.find(tempFieldName);
  if( iter == fields.end() )
  {
    std::cout << "ERROR: Cannot find Temperature field \"" << tempFieldName
         << "\" for use in computeSourceTerm for problem \""
         << getName() << std::endl;
    throw "HMX_PDE ERROR";
  }
  else
    TvecPtr = (*iter).second;

  Epetra_Vector &T = *TvecPtr;

  // If this problem is the temperature equation, don't compute a source
  // term.  This would be where a volumetric heating term would go.
  if( getName() == tempFieldName )
  {
    result.PutScalar(0.0);
    return;
  }
  else
  {
    double rateK;

    map<string, double>::iterator requiredFieldIter;
    map<string, double>::iterator requiredFieldEnd = SrcTermExponent.end();

    for( int i = 0; i<result.MyLength(); i++ ) {

      rateK = pow(T[i],StericCoef) * PreExp_A *
              exp( -ActEnergy / (Const_R * T[i]) );

      result[i] = 1.0;  // Start point for product

      // Loop over all required fields and contribute to product
      for( requiredFieldIter = SrcTermExponent.begin();
           requiredFieldIter != requiredFieldEnd; requiredFieldIter++)
      {
        iter = fields.find( (*requiredFieldIter).first );
        if( iter == fields.end() )
        {
          std::cout << "ERROR: Cannot find required field \""
               << (*requiredFieldIter).first
               << "\" for use in computeSourceTerm for problem \""
               << getName() << std::endl;
          throw "HMX_PDE ERROR";
        }
    Epetra_Vector &reqFieldVec = *((*iter).second);

        result[i] *= pow( reqFieldVec[i], (*requiredFieldIter).second );
      }
      result[i] *= rateK;
    }
  }
}

//-----------------------------------------------------------------------------

void
HMX_PDE::computeSourceTerm(int arraySize,
                    map<string, double*> fields,
                    double*& result)
{
  // All dependent variables should be in place before now

  // We assume all other registered dependent problems are species which
  // affect the reaction rate of this specie

  double *T = 0;
  map<string, double*>::iterator iter = fields.find(tempFieldName);
  if( iter == fields.end() ) {
    std::cout << "ERROR: Cannot find Temperature field \"" << tempFieldName
         << "\" for use in computeSourceTerm for problem \""
         << getName() << std::endl;
    throw "HMX_PDE ERROR";
  }
  else
    T = (*iter).second;

  // If this problem is the temperature equation, don't compute a source
  // term.  This would be where a volumetric heating term would go.
  if( getName() == tempFieldName ) {
    for( int i = 0; i<arraySize; i++)
      result[i] = 0.0;

    return;
  }
  else {

    double rateK;

    map<string, double>::iterator requiredFieldIter;
    map<string, double>::iterator requiredFieldEnd = SrcTermExponent.end();

    for( int i = 0; i<arraySize; i++ ) {

      rateK = pow(T[i],StericCoef) * PreExp_A *
              exp( -ActEnergy / (Const_R * T[i]) );

      result[i] = 1.0;  // Start point for product

      // Loop over all required fields and contribute to product
      for( requiredFieldIter = SrcTermExponent.begin();
           requiredFieldIter != requiredFieldEnd; requiredFieldIter++) {

        iter = fields.find( (*requiredFieldIter).first );
        if( iter == fields.end() ) {
          std::cout << "ERROR: Cannot find required field \""
               << (*requiredFieldIter).first
               << "\" for use in computeSourceTerm for problem \""
               << getName() << "\"" << std::endl;
          throw "HMX_PDE ERROR";
        }
    double *reqFieldVec = (*iter).second;

        result[i] *= pow( reqFieldVec[i], (*requiredFieldIter).second );
      }
      result[i] *= rateK;
    }
  }
}

//-----------------------------------------------------------------------------

void
HMX_PDE::compute_dRdT()
{
  // Compute dependence of Rxn rate on Temperature
}

//-----------------------------------------------------------------------------

void
HMX_PDE::compute_dRdN()
{
  // Compute dependence of Rxn rate on Self-Species
}

//-----------------------------------------------------------------------------

