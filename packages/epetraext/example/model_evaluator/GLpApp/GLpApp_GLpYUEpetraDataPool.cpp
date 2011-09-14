/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "GLpApp_GLpYUEpetraDataPool.hpp"
#include "rect2DMeshGenerator.hpp"

#include "Amesos_Klu.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_Reindex_LinearProblem.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Epetra_BLAS.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif


// 2008/09/04: Added to address failed tests (see bug 4040)
using namespace std;


// Define to see all debug output for mesh generation
//#define GLPYUEPETRA_DATAPOOL_DUMP_ALL_MESH


namespace GLpApp {

//
// Declarations
//

const double GLp_pi = 3.14159265358979323846;

ostream& operator <<(ostream &, const Usr_Par &);

bool CrsMatrix2MATLAB(const Epetra_CrsMatrix &, ostream &);

bool Vector2MATLAB( const Epetra_Vector &, ostream &);

bool FEVector2MATLAB( const Epetra_FEVector &, ostream &);

int quadrature(const int, const int, Epetra_SerialDenseMatrix &,
               Epetra_SerialDenseVector &);

int meshreader(
  const Epetra_Comm &,
  Epetra_IntSerialDenseVector &,
  Epetra_SerialDenseMatrix &,
  Epetra_IntSerialDenseVector &,
  Epetra_SerialDenseMatrix &,
  Epetra_IntSerialDenseMatrix &,
  Epetra_IntSerialDenseMatrix &,
  const char geomFileBase[],
  const bool trace = true,
  const bool dumpAll = false
  );

int lassembly(const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseVector &,
              const Usr_Par &,
              Epetra_SerialDenseMatrix &,
              Epetra_SerialDenseVector &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &,
                  bool);

int assemble(const Epetra_Comm &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemble_bdry(
  const Epetra_Comm                                &Comm
  ,const Epetra_IntSerialDenseVector               &ipindx
  ,const Epetra_SerialDenseMatrix                  &ipcoords
  ,const Epetra_IntSerialDenseVector               &pindx
  ,const Epetra_SerialDenseMatrix                  &pcoords
  ,const Epetra_IntSerialDenseMatrix               &t
  ,const Epetra_IntSerialDenseMatrix               &e
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *B
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *R
  );

int nonlinvec(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FEVector> &);

int nonlinjac(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FECrsMatrix> &);

int nonlinhessvec(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);

int compproduct(Epetra_SerialDenseVector &, double *, double *);

int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);

double determinant(const Epetra_SerialDenseMatrix &);

int inverse(const Epetra_SerialDenseMatrix &, Epetra_SerialDenseMatrix &);

int quadrature(
  const int, const int, Epetra_SerialDenseMatrix &,
  Epetra_SerialDenseVector &);

void gpfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv);

void g2pfctn(const Epetra_SerialDenseVector & , Epetra_SerialDenseVector & );

void gfctn(const Epetra_SerialDenseVector & , Epetra_SerialDenseVector & );

//
// GLpYUEpetraDataPool
//

GLpYUEpetraDataPool::GLpYUEpetraDataPool(
  Teuchos::RefCountPtr<const Epetra_Comm>    const& commptr
  ,const double                              beta
  ,const double                              len_x
  ,const double                              len_y
  ,const int                                 local_nx
  ,const int                                 local_ny
  ,const char                                myfile[]
  ,const bool                                trace
  )
  :commptr_(commptr)
  ,beta_(beta)
{
  ipcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  ipindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  pcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  pindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  t_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );
  e_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );

  if( myfile && myfile[0] != '\0' ) {
    meshreader(*commptr_,*ipindx_,*ipcoords_,*pindx_,*pcoords_,*t_,*e_,myfile,trace);
  }
  else {
    rect2DMeshGenerator(
      commptr_->NumProc(),commptr_->MyPID()
      ,len_x,len_y,local_nx,local_ny,2 // 2==Neuman boundary conditions!
      ,&*ipindx_,&*ipcoords_,&*pindx_,&*pcoords_,&*t_,&*e_
#ifdef GLPYUEPETRA_DATAPOOL_DUMP_ALL_MESH
      ,&*Teuchos::VerboseObjectBase::getDefaultOStream(),true
#else
      ,NULL,false
#endif
      );
  }

  // Assemble volume and boundary mass and stiffness matrices, and the right-hand side of the PDE.
  assemble( *commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, A_, H_, b_ );
  assemble_bdry( *commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, &B_, &R_ );

  // Set desired state q.
  Epetra_Map standardmap(A_->DomainMap());
  q_ = Teuchos::rcp(new Epetra_FEVector(standardmap,1));
  int * qintvalues = new int[standardmap.NumMyElements()];
  double * qdvalues = new double[standardmap.NumMyElements()];
  standardmap.MyGlobalElements(qintvalues);
  for (int i = 0; i < standardmap.NumMyElements(); i++)
      qdvalues[i]=cos( GLp_pi* ((*ipcoords_)(i,0)) ) * cos( GLp_pi* ((*ipcoords_)(i,1)) );
  q_->ReplaceGlobalValues(standardmap.NumMyElements(), qintvalues, qdvalues);
  q_->GlobalAssemble();
}

void GLpYUEpetraDataPool::computeAll( const GenSQP::Vector &x )
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
        (Teuchos::dyn_cast<GenSQP::YUEpetraVector>(const_cast<GenSQP::Vector&>(x))).getYVector();

  computeNy(ey);

  computeNpy(ey);

  computeAugmat();
  
}

int GLpYUEpetraDataPool::solveAugsys(
  const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
  const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
  const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
  const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
  const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
  const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
  double tol )
{
/*
  int systemChoice = 1;   // 1 for full KKT system solve, 2 for Schur complement solve
  int solverChoice = 12;  // These options are for the full KKT system solve.
                          // 11 for AztecOO with built-in Schwarz DD preconditioning and ILU on subdomains
                          // 12 for AztecOO with IFPACK Schwarz DD preconditioning and Umfpack on subdomains
                          // 13 for a direct sparse solver (Umfpack, KLU)
  
  if (systemChoice == 1) {
    // We're using the full KKT system formulation to solve the augmented system.
   
    Epetra_Map standardmap(A_->DomainMap());
    int numstates = standardmap.NumGlobalElements();
    Epetra_Map bdryctrlmap(B_->DomainMap());
    int numcontrols = bdryctrlmap.NumGlobalElements();
    Epetra_Vector rhs( (Epetra_BlockMap&)Augmat_->RangeMap() );
    Epetra_Vector soln( (Epetra_BlockMap&)Augmat_->RangeMap() );
    soln.Random();  

    std::vector<double> values(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    std::vector<int> indices(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    ((Epetra_BlockMap&)Augmat_->RangeMap()).MyGlobalElements(&indices[0]);

    for (int i=0; i<rhsy->MyLength(); i++) {
      values[i] = (*((*rhsy)(0)))[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      values[i+rhsy->MyLength()] = (*((*rhsu)(0)))[i];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      values[i+rhsy->MyLength()+rhsu->MyLength()] = (*((*rhsp)(0)))[i];
    }

    rhs.ReplaceGlobalValues(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength(), &values[0], &indices[0]);

    if (solverChoice == 11) {
      int Overlap = 3;
      int ival = 4;

      AztecOO::AztecOO kktsolver(&(*Augmat_), &soln, &rhs);
      kktsolver.SetAztecOption( AZ_solver, AZ_gmres );
      kktsolver.SetAztecOption( AZ_precond, AZ_dom_decomp );
      kktsolver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      //kktsolver.SetAztecOption( AZ_kspace, 2*numstates+numcontrols );
      kktsolver.SetAztecOption( AZ_kspace, 9000 );
      kktsolver.SetAztecOption(AZ_overlap,Overlap);
      kktsolver.SetAztecOption(AZ_graph_fill,ival);
      //kktsolver.SetAztecOption(AZ_poly_ord, ival);
      //kktsolver.SetAztecParam(AZ_drop, 1e-9);
      kktsolver.SetAztecParam(AZ_athresh, 1e-5);
      //kktsolver.SetAztecParam(AZ_rthresh, 0.0);
      kktsolver.SetAztecOption( AZ_reorder, 0 );
      //kktsolver.SetAztecParam44( AZ_ilut_fill, 1.5 );
      kktsolver.SetAztecOption( AZ_output, AZ_last );
      //kktsolver.Iterate(2*numstates+numcontrols,1e-12);
      kktsolver.Iterate(9000,1e-11);
      //cout << soln;
    }
    else if (solverChoice == 12) {
      // =============================================================== //
      // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
      // =============================================================== //

      Teuchos::ParameterList List;

      // allocates an IFPACK factory. No data is associated
      // to this object (only method Create()).
      Ifpack Factory;

      // create the preconditioner. For valid PrecType values,
      // please check the documentation
      string PrecType = "Amesos";
      int OverlapLevel = 2; // must be >= 0. If Comm.NumProc() == 1,
                            // it is ignored.
  
      Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &(*Augmat_), OverlapLevel);
      assert(Prec != 0);

      // specify the Amesos solver to be used.
      // If the selected solver is not available,
      // IFPACK will try to use Amesos' KLU (which is usually always
      // compiled). Amesos' serial solvers are:
      // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
      List.set("amesos: solver type", "Amesos_Umfpack");

      // sets the parameters
      IFPACK_CHK_ERR(Prec->SetParameters(List));

      // initialize the preconditioner. At this point the matrix must
      // have been FillComplete()'d, but actual values are ignored.
      // At this call, Amesos will perform the symbolic factorization.
      IFPACK_CHK_ERR(Prec->Initialize());

      // Builds the preconditioners, by looking for the values of
      // the matrix. At this call, Amesos will perform the
      // numeric factorization.
      IFPACK_CHK_ERR(Prec->Compute());

      // =================================================== //
      // E N D   O F   I F P A C K   C O N S T R U C T I O N //
      // =================================================== //

      // need an Epetra_LinearProblem to define AztecOO solver
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);

      // now we can allocate the AztecOO solver
      AztecOO kktsolver(Problem);

      // specify solver
      kktsolver.SetAztecOption(AZ_solver,AZ_gmres);
      kktsolver.SetAztecOption(AZ_kspace, 300 );
      kktsolver.SetAztecOption(AZ_output,AZ_last);

      // HERE WE SET THE IFPACK PRECONDITIONER
      kktsolver.SetPrecOperator(Prec);

      // .. and here we solve
      kktsolver.Iterate(300,1e-12);

      // delete the preconditioner
      delete Prec;
    }
    else if (solverChoice == 13) {
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);
      
      EpetraExt::LinearProblem_Reindex reindex(NULL);
      Epetra_LinearProblem newProblem = reindex(Problem);
      
      Amesos_Klu kktsolver(newProblem);
   
      AMESOS_CHK_ERR(kktsolver.SymbolicFactorization());
      AMESOS_CHK_ERR(kktsolver.NumericFactorization());
      AMESOS_CHK_ERR(kktsolver.Solve());
      kktsolver.PrintTiming();
    }
    
    
    for (int i=0; i<rhsy->MyLength(); i++) {
      (*((*y)(0)))[i] = soln[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      (*((*u)(0)))[i] = soln[i+rhsy->MyLength()];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      (*((*p)(0)))[i] = soln[i+rhsy->MyLength()+rhsu->MyLength()];
    }
    
  }
  else if (systemChoice == 2) {
    // We're using the Schur complement formulation to solve the augmented system.
  
    // Form linear operator.
    GLpApp::SchurOp schurop(A_, B_, Npy_);
  
    // Form Schur complement right-hand side.
    Epetra_MultiVector ny( (Epetra_BlockMap&)Npy_->RangeMap(), 1);
    Epetra_MultiVector ay( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector schurrhs( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector bu( (Epetra_BlockMap&)B_->RangeMap(), 1);
    A_->Multiply(false, *rhsy, ay);
    Npy_->Multiply(false, *rhsy, ny);
    B_->Multiply(false, *rhsu, bu);
    schurrhs.Update(1.0, ny, 1.0, ay, 0.0);
    schurrhs.Update(-1.0, *rhsp, 1.0, bu, 1.0);
  
    p->PutScalar(0.0);
    Epetra_LinearProblem linprob(&schurop, &(*p), &schurrhs);
    AztecOO::AztecOO Solver(linprob);
    Solver.SetAztecOption( AZ_solver, AZ_cg );
    Solver.SetAztecOption( AZ_precond, AZ_none );
    Solver.SetAztecOption( AZ_output, AZ_none );
    Solver.Iterate(8000,tol);
  
    Epetra_MultiVector bp( (Epetra_BlockMap&)B_->DomainMap(), 1);
    B_->Multiply(true, *p, bp);
    u->Update(1.0, *rhsu, -1.0, bp, 0.0);

    Epetra_MultiVector ap( (Epetra_BlockMap&)A_->DomainMap(), 1);
    Epetra_MultiVector np( (Epetra_BlockMap&)A_->DomainMap(), 1);
    A_->Multiply(true, *p, ap);
    Npy_->Multiply(true, *p, np);
    y->Update(1.0, *rhsy,0.0);
    y->Update(-1.0, ap, -1.0, np, 1.0);
  }
*/
  TEST_FOR_EXCEPT(true);
  return 0;
}

Teuchos::RefCountPtr<const Epetra_Comm> GLpYUEpetraDataPool::getCommPtr()   { return commptr_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getA()  { return A_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getB()  { return B_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getH()  { return H_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getR()  { return R_; }

Teuchos::RefCountPtr<Epetra_CrsMatrix> GLpYUEpetraDataPool::getAugmat()  { return Augmat_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getNpy()  { return Npy_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getb()  { return b_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getq()  { return q_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getNy()  { return Ny_; }

double GLpYUEpetraDataPool::getbeta()  { return beta_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLpYUEpetraDataPool::getipcoords()  { return ipcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLpYUEpetraDataPool::getipindx()  { return ipindx_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLpYUEpetraDataPool::getpcoords()  { return pcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLpYUEpetraDataPool::getpindx()  { return pindx_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLpYUEpetraDataPool::gett()  { return t_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLpYUEpetraDataPool::gete()  { return e_; }


void GLpYUEpetraDataPool::computeNy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinvec(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Ny_);
}


void GLpYUEpetraDataPool::computeNpy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinjac(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Npy_);
}


void GLpYUEpetraDataPool::computeAugmat()
{
  Epetra_Map standardmap(A_->DomainMap());
  Epetra_Map bdryctrlmap(B_->DomainMap());

  int indexBase = 1;

  int numstates = standardmap.NumGlobalElements();
  //int numcontrols = bdryctrlmap.NumGlobalElements();
  int nummystates = standardmap.NumMyElements();
  int nummycontrols = bdryctrlmap.NumMyElements();

  Epetra_IntSerialDenseVector KKTmapindx(2*nummystates+nummycontrols);
  
  // Build KKT map.
  Epetra_IntSerialDenseVector states(nummystates);
  Epetra_IntSerialDenseVector controls(nummycontrols);
  standardmap.MyGlobalElements(states.Values());
  bdryctrlmap.MyGlobalElements(controls.Values());
  for (int i=0; i<nummystates; i++) {
    KKTmapindx(i) = states(i);
    KKTmapindx(nummystates+nummycontrols+i) = 2*numstates+states(i);
  }
  for (int i=0; i<nummycontrols; i++) {
    KKTmapindx(nummystates+i) = numstates+controls(i);
  }
  Epetra_Map KKTmap(-1, 2*nummystates+nummycontrols, KKTmapindx.Values(), indexBase, *(commptr_));
  
  
  // Start building the KKT matrix.
  
  Augmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, KKTmap, 0));

  double one[1];
  one[0]=1.0;
  for (int i=0; i<nummystates+nummycontrols; i++) {
    Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], 1, one, KKTmapindx.Values()+i);
  }
  
  int checkentries=0;
  int nummyentries=0;
  Epetra_SerialDenseVector values(nummyentries);
  Epetra_IntSerialDenseVector indices(nummyentries);
  // Insert A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = A_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    A_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = Npy_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    Npy_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = B_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    B_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  
  bool MakeDataContiguous = false;
  EpetraExt::RowMatrix_Transpose transposer( MakeDataContiguous );
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A_));
  Epetra_CrsMatrix & transB = dynamic_cast<Epetra_CrsMatrix&>(transposer(*B_));
  Epetra_CrsMatrix & transNpy = dynamic_cast<Epetra_CrsMatrix&>(transposer(*Npy_));
  // Insert transpose of A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transA.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transA.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = transNpy.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transNpy.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                  indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert transpose of B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transB.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transB.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    int err = 0;
    if (nummyentries > 0)
      err = Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+numstates, nummyentries,
                                        values.Values(), indices.Values());
    // This will give a nasty message if something goes wrong with the insertion of B transpose.
    if (err < 0) {
      cout << "Insertion of entries failed:\n";
      cout << indices;
      cout << nummyentries << endl;
      cout << "at row: " << KKTmapindx.Values()[i]+numstates << endl << endl;
    }
  }

  Augmat_->FillComplete(KKTmap, KKTmap);
  // End building the KKT matrix.

}

void GLpYUEpetraDataPool::PrintVec( const Teuchos::RefCountPtr<const Epetra_Vector> & x )
{
  Vector2MATLAB(*x, cout);
}

Usr_Par::Usr_Par() {
  Epetra_BLAS blas;
  Epetra_SerialDenseVector product(4);

  // get nodes and weights
  quadrature(2,3,Nodes,Weights);
    
  // Evaluate nodal basis functions and their derivatives at quadrature
  // pts N(i,j) = value of the j-th basis function at quadrature node i.
  N.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }
  // Nx1(i,j) partial derrivative of basis function j wrt x1 at node i
  Nx1.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    Nx1(i,0) = -1.0;
    Nx1(i,1) = 1.0;
    Nx1(i,2) = 0.0;
  }
  // Nx2(i,j) partial derrivative of basis function j wrt x2 at node i
  Nx2.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    Nx2(i,0) = -1.0;
    Nx2(i,1) = 0.0;
    Nx2(i,2) = 1.0;
  }

  // Evaluate integrals of various combinations of partial derivatives
  // of the local basis functions (they're constant).
  S1.Shape(3,3);
  S1(0,0)= 1.0; S1(0,1)=-1.0; S1(0,2)= 0.0;
  S1(1,0)=-1.0; S1(1,1)= 1.0; S1(1,2)= 0.0;
  S1(2,0)= 0.0; S1(2,1)= 0.0; S1(2,2)= 0.0;
  S2.Shape(3,3);
  S2(0,0)= 2.0; S2(0,1)=-1.0; S2(0,2)=-1.0;
  S2(1,0)=-1.0; S2(1,1)= 0.0; S2(1,2)= 1.0;
  S2(2,0)=-1.0; S2(2,1)= 1.0; S2(2,2)= 0.0;
  S3.Shape(3,3);
  S3(0,0)= 1.0; S3(0,1)= 0.0; S3(0,2)=-1.0;
  S3(1,0)= 0.0; S3(1,1)= 0.0; S3(1,2)= 0.0;
  S3(2,0)=-1.0; S3(2,1)= 0.0; S3(2,2)= 1.0;
    
  // Evaluate integrals of basis functions N(i).
  Nw.Size(3);
  Nw.Multiply('T', 'N', 1.0, N, Weights, 0.0);

  // Evaluate integrals of 9 products N(i)*N(j).
  NNw.Shape(3,3);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      compproduct(product, N[i], N[j]);
      NNw(i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
    }
  }

  // Evaluate integrals of 27 products N(i)*N(j)*N(k).
  NNNw = new Epetra_SerialDenseMatrix[3];
  NNNw[0].Shape(3,3); NNNw[1].Shape(3,3); NNNw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], N[j], N[k]);
        NNNw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

  // Evaluate integrals of 27 products N(i)*(dN(j)/dx1)*N(k).
  NdNdx1Nw = new Epetra_SerialDenseMatrix[3];
  NdNdx1Nw[0].Shape(3,3); NdNdx1Nw[1].Shape(3,3); NdNdx1Nw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], Nx1[j], N[k]);
        NdNdx1Nw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

  // Evaluate integrals of 27 products N(i)*(dN(j)/dx2)*N(k).
  NdNdx2Nw = new Epetra_SerialDenseMatrix[3];
  NdNdx2Nw[0].Shape(3,3); NdNdx2Nw[1].Shape(3,3); NdNdx2Nw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], Nx2[j], N[k]);
        NdNdx2Nw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

}

void Usr_Par::Print(ostream& os) const {
  os << endl;
  os << "\n\nQuadrature nodes:";
  os << "\n-----------------";
  Nodes.Print(os);
  os << "\n\nQuadrature weights:";
  os << "\n-------------------\n";
  Weights.Print(os);
  os << "\n\nNodal basis functions:";
  os << "\n----------------------";
  N.Print(os);
  os << "\n\nIntegrals of combinations of partial derivatives:";
  os << "\n-------------------------------------------------";
  S1.Print(os);
  S2.Print(os);
  S3.Print(os);
  os << "\n\nIntegrals of basis functions:";
  os << "\n-----------------------------\n";
  Nw.Print(os);
  os << "\n\nIntegrals of products N(i)*N(j):";
  os << "\n--------------------------------\n";
  NNw.Print(os);
  os << "\n\nIntegrals of products N(i)*N(j)*N(k):";
  os << "\n-------------------------------------\n";
  NNNw[0].Print(os); NNNw[1].Print(os); NNNw[2].Print(os);
  os << "\n\nIntegrals of products N(i)*(dN(j)/dx1)*N(k):";
  os << "\n--------------------------------------------\n";
  NdNdx1Nw[0].Print(os); NdNdx1Nw[1].Print(os); NdNdx1Nw[2].Print(os);
  os << "\n\nIntegrals of products N(i)*(dN(j)/dx2)*N(k):";
  os << "\n--------------------------------------------\n";
  NdNdx2Nw[0].Print(os); NdNdx2Nw[1].Print(os); NdNdx2Nw[2].Print(os);
}

ostream& operator <<(ostream & out, const Usr_Par & usr_par)
{
  usr_par.Print(out);
  return out;
}

} // namespace GLpApp

//
// Implementations
//

int GLpApp::compproduct( Epetra_SerialDenseVector & product,
  double *first, double *second)
{
  for (int i=0; i<product.M(); i++) {
    product[i] = first[i]*second[i];
  }
  return(0);
}

int GLpApp::compproduct(Epetra_SerialDenseVector & product,
                double *first, double *second, double *third)
{
  for (int i=0; i<product.M(); i++) {
    product[i] = first[i]*second[i]*third[i];
  }
  return(0);
}

//#define GLPAPP_SHOW_BOUNDARY_ASSEMBLY

/*  \brief Performs finite-element assembly of face mass matrices.

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the 0-based j-th vertex in triangle i, where i = 0, ..., numElements-1
  \param  e         [in]  - Matrix (EDGE x 3) of edges. \n
                            e(i,1:2) contains the indices of the endpoints of edge i, where i = 0, ..., numEdges-1 \n
                            e(i,3) contains the boundary marker
  \param  B_out     [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
                            state/control face mass matrix.
  \param  R_out     [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
                            control/control volume mass matrix.
  \return 0                 if successful.

  \par Detailed Description:

  -# Assembles the FE state/control mass matrix \e B, given by
     \f[
       \mathbf{B}_{jk} = b(\mu_k,\phi_j) = - \int_{\partial \Omega}  \mu_k(x) \phi_j(x) dx,
     \f]
     where \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the piecewise linear nodal basis for the finite-dimensional
     state space, and \f$\{ \mu_j \}_{j = 1}^{n}\f$ is the piecewise linear nodal basis for the
     finite-dimensional control space.
  -# Assembles the FE control/control mass matrix \e R, given by
     \f[
       \mathbf{R}_{jk} = b(\mu_k,\mu_j) = - \int_{\partial \Omega}  \mu_k(x) \mu_j(x) dx,
     \f]
     where \f$\{ \mu_j \}_{j = 1}^{n}\f$ is the piecewise linear nodal basis for the finite-dimensional
     control space.
*/
int GLpApp::assemble_bdry(
  const Epetra_Comm                                &Comm
  ,const Epetra_IntSerialDenseVector               &ipindx
  ,const Epetra_SerialDenseMatrix                  &ipcoords
  ,const Epetra_IntSerialDenseVector               &pindx
  ,const Epetra_SerialDenseMatrix                  &pcoords
  ,const Epetra_IntSerialDenseMatrix               &t
  ,const Epetra_IntSerialDenseMatrix               &e
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *B_out
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *R_out
  )
{

  using Teuchos::rcp;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out << "\nEntering assemble_bdry(...) ...\n";
#endif

  int numLocElems = t.M();
  int numLocEdges = e.M();

  int indexBase = 1;

  //
  // Create a sorted (by global ID) list of boundry nodes
  //
  int * lastindx = 0;
  Epetra_IntSerialDenseVector BdryNode(2*numLocEdges);
  for (int j=0; j<numLocEdges; j++) {
    BdryNode[j] = e(j,0);
    BdryNode[j+numLocEdges] = e(j,1);
  }
  std::sort(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  lastindx  = std::unique(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  const int numBdryNodes = lastindx - BdryNode.Values();
  BdryNode.Resize(numBdryNodes); // RAB: This does not overwrite?

  //
  // Determine the boundary vertices that belong to this processor.
  //
  Epetra_IntSerialDenseVector MyBdryNode(numBdryNodes);
  lastindx = std::set_intersection(
    BdryNode.Values(), BdryNode.Values()+numBdryNodes,  // (Sorted) set A
    ipindx.Values(), ipindx.Values()+ipindx.M(),        // (Sorted) set B
    MyBdryNode.Values()                                 // Intersection
    );
  const int numMyBdryNodes = lastindx - MyBdryNode.Values();
  MyBdryNode.Resize(numMyBdryNodes); // RAB: This does not overwrite?
  
  //
  // Define the maps for the various lists
  //
  Epetra_Map standardmap(-1, ipindx.M(), const_cast<int*>(ipindx.A()), indexBase, Comm);
  Epetra_Map overlapmap(-1, pindx.M(), const_cast<int*>(pindx.A()), indexBase, Comm);
  Epetra_Map mybdryctrlmap(-1, numMyBdryNodes, const_cast<int*>(MyBdryNode.A()), indexBase, Comm);
  // Above, it is important to note what mybndyctrlmap represents.  It is the
  // a sorted list of global vertex node IDS for nodes on a boundary that are
  // uniquely owned by the local process.

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  *out << "\nstandardmap:\n";
  standardmap.Print(Teuchos::OSTab(out).o());
  *out << "\nmybdryctrlmap:\n";
  mybdryctrlmap.Print(Teuchos::OSTab(out).o());
#endif

  //
  // Allocate matrices to fill
  //
  Teuchos::RefCountPtr<Epetra_FECrsMatrix>
    B = rcp(new Epetra_FECrsMatrix(Copy,standardmap,0)),
    R = rcp(new Epetra_FECrsMatrix(Copy,standardmap,0));
  // NOTE: The data map is the same as for the volume matrices. Later, when
  // FillComplete is called, we will fix the proper domain and range maps. 
  // Declare quantities needed for the call to the local assembly routine.
  const int numNodesPerEdge = 2;
  Epetra_IntSerialDenseVector epetra_nodes(numNodesPerEdge);

  //
  // Load B and R by looping through the edges
  //

  Epetra_SerialDenseMatrix Bt(2,2);
  int err=0;

  for( int i=0; i < numLocEdges; i++ ) {

    const int
      global_id_0 = e(i,0),
      global_id_1 = e(i,1),
      local_id_0  = overlapmap.LID(global_id_0), // O(log(numip)) bindary search
      local_id_1  = overlapmap.LID(global_id_1); // O(log(numip)) bindary search

    epetra_nodes(0) = global_id_0;
    epetra_nodes(1) = global_id_1;

    const double
      x0 = pcoords(local_id_0,0),
      y0 = pcoords(local_id_0,1),
      x1 = pcoords(local_id_1,0),
      y1 = pcoords(local_id_1,1);
    
    const double l = sqrt(pow(x0-x1,2) + pow(y0-y1,2));  // Length of this edge
    
    // We have an explicit formula for Bt.
    const double l_sixth = l * (1.0/6.0);
    Bt(0,0) = l_sixth * 2.0;
    Bt(0,1) = l_sixth * 1.0;
    Bt(1,0) = l_sixth * 1.0;
    Bt(1,1) = l_sixth * 2.0;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
    *out
      << "\nEdge global nodes = ("<<global_id_0<<","<<global_id_1<<"),"
      << " local nodes = ("<<local_id_0<<","<<local_id_1<<"),"
      << " (x0,y0)(x1,y1)=("<<x0<<","<<y0<<")("<<x1<<","<<y1<<"),"
      << " Bt = ["<<Bt(0,0)<<","<<Bt(0,1)<<";"<<Bt(1,0)<<","<<Bt(1,1)<<"]\n";
#endif

    const int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
    err = B->InsertGlobalValues(epetra_nodes,Bt,format);
    if (err<0) return(err);
    err = R->InsertGlobalValues(epetra_nodes,Bt,format);
    if (err<0) return(err);
    
  }

/*

  err = B->GlobalAssemble(false);
  if (err<0) return(err);
  err = R->GlobalAssemble(false);
  if (err<0) return(err);

  err = B->FillComplete(mybdryctrlmap,standardmap);
  if (err<0) return(err);
  err = R->FillComplete(mybdryctrlmap,mybdryctrlmap);
  if (err<0) return(err);

*/

  err = B->GlobalAssemble(mybdryctrlmap,standardmap,true);
  if (err<0) return(err);
  err = R->GlobalAssemble(mybdryctrlmap,mybdryctrlmap,true);
  if (err<0) return(err);

  if(B_out) *B_out = B;
  if(R_out) *R_out = R;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  *out << "\nB =\n";
  B->Print(Teuchos::OSTab(out).o());
  *out << "\nLeaving assemble_bdry(...) ...\n";
#endif

  return(0);

}

/*  \brief Performs finite-element assembly of volume stiffness and mass matrices,
            and the right-hand side (forcing term).

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the j-th vertex in triangle i, where i = 1, ..., ELE
  \param  e         [in]  - Matrix (EDGE x 3) of edges. \n
                            e(i,1:2) contains the indices of the endpoints of edge i, where i = 1, ..., EDGE \n
                            e(i,3) contains the boundary marker \n
                            e(i,3) = 1  Dirichlet boundary conditions on edge i \n
                            e(i,3) = 2  Neumann boundary conditions on edge i \n
                            e(i,3) = 3  Robin boundary conditions on edge i
  \param  A         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE volume
                            stiffness matrix for the advection-diffusion equation. Includes advection,
                            diffusion, and reaction terms, and modifications that account for the boundary
                            conditions.
  \param  M         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE volume
                            mass matrix.
  \param  b         [out] - Reference-counting pointer to the Epetra_FEVector describing the discretized
                            forcing term in the advection-diffusion equation. Includes modifications that
                            account for the boundary conditions.
  \return 0                 if successful.

  \par Detailed Description:

  -# Assembles the FE volume stiffness matrix \e A and the right-hand side \e b for the
  solution of an advection-diffusion equation using piecewise linear finite elements.
  The advection-diffusion equation is 
  \f{align*}
       - \nabla \cdot \left( K(x) \nabla y(x) \right) + \mathbf{c}(x) \cdot \nabla y(x) + r(x) y(x)
       &= f(x), & x &\in \Omega,\;  \\
       y(x) &= d(x),    & x &\in {\partial \Omega}_d, \\
       K(x) \frac{\partial}{\partial \mathbf{n}}  y(x) &= g(x), & x &\in {\partial \Omega}_n, \\
       \sigma_0 K(x) \frac{\partial}{\partial \mathbf{n}} y(x)
       + \sigma_1 y(x) &= g(x), & x &\in {\partial \Omega}_r,
  \f}
  where \f$ K \f$ represents scalar diffusion, \f$ \mathbf{c} \f$ is the advection vector field,
  \f$ r \f$ is the reaction term, \f$ d \f$ and  \f$ g \f$ are given functions, \f$sigma_0\f$ and
  \f$ sigma_1 \f$ are given quantities, \f$ {\partial \Omega}_d \f$ is the Dirichlet boundary,
  \f$ {\partial \Omega}_n \f$ is the Neumann boundary, and \f$ {\partial \Omega}_r \f$ is the Robin
  boundary. The quantities \f$ K \f$, \f$ \mathbf{c} \f$, \f$ r \f$, \f$ d \f$, and \f$ g \f$ are
  assumed to be piecewise linear. Currently, they are to be hard-coded inside this function.
  -# Assembles the FE volume mass matrix \e M.
*/
int GLpApp::assemble(const Epetra_Comm & Comm,
             const Epetra_IntSerialDenseVector & ipindx,
             const Epetra_SerialDenseMatrix & ipcoords,
             const Epetra_IntSerialDenseVector & pindx,
             const Epetra_SerialDenseMatrix & pcoords,
             const Epetra_IntSerialDenseMatrix & t,
             const Epetra_IntSerialDenseMatrix & e,
             RefCountPtr<Epetra_FECrsMatrix> & A,
             RefCountPtr<Epetra_FECrsMatrix> & M,
             RefCountPtr<Epetra_FEVector> & b)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();
  Usr_Par usr_par;

  int numLocElems = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, ipindx.M(), (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, pindx.M(), (int*)pindx.A(), indexBase, Comm);

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;

  A = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  M = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  b = rcp(new Epetra_FEVector(standardmap,1));

  // Declare quantities needed for the call to the local assembly routine.
  int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);


  /* ************************  First, we build A and b.  ************************ */
  Epetra_SerialDenseMatrix At;
  Epetra_SerialDenseVector bt;
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());
  
  Epetra_SerialDenseVector k(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) k(i)=1.0;
  Epetra_SerialDenseMatrix c(numNodesPerElem,2);
  for (i=0; i< numNodesPerElem; i++) { c(i,0)=0.0; c(i,1)=0.0; }
  Epetra_SerialDenseVector r(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) r(i)=0.0;
  Epetra_SerialDenseVector f(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) f(i)=0.0;
  Epetra_SerialDenseVector g(2);
  g(0)=0.0; g(1)=0.0;
  Epetra_SerialDenseVector sigma(2);
  sigma(0)=0.0; sigma(1)=0.0;
  for(i=0; i<numLocElems; i++) {
    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      // get vertex
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
      // set rhs function
      f(j) = cos(GLp_pi*vertices(j,0))*cos(GLp_pi*vertices(j,1)) * 
               (2*GLp_pi*GLp_pi + pow(cos(GLp_pi*vertices(j,0)),2)*pow(cos(GLp_pi*vertices(j,1)),2) - 1.0);
    }
    lassembly(vertices, k, c, r, f, usr_par, At, bt);
    err = A->InsertGlobalValues(epetra_nodes, At, format);
    if (err<0) return(err);
    err = b->SumIntoGlobalValues(epetra_nodes, bt);
    if (err<0) return(err);
  }

  // MAKE ADJUSTMENTS TO A and b TO ACCOUNT FOR BOUNDARY CONDITIONS.

  // Get Neumann boundary edges.
  Epetra_IntSerialDenseMatrix neumann(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==2) {
      neumann(j,0) = e(i,0);  neumann(j,1) = e(i,1);
      j++;
    }
  }
  neumann.Reshape(j,2);
  // Adjust for Neumann BC's.
  if (neumann.M() != 0) {
    Epetra_BLAS blas;
    Epetra_SerialDenseMatrix quadnodes;
    Epetra_SerialDenseVector quadweights;
    Epetra_SerialDenseMatrix N;
    Epetra_SerialDenseMatrix NN;
    Epetra_SerialDenseVector product(2);

    // Get quadrature weights and points.
    quadrature(1, 2, quadnodes, quadweights);
    
    // Evaluate nodal basis functions at quadrature points
    // N(i,j) value of basis function j at quadrature node i
    N.Shape(quadnodes.M(),2);
    for (i=0; i<quadnodes.M(); i++) {
      N(i,0) = 1.0 - quadnodes(i,0);
      N(i,1) = quadnodes(i,0);
    }

    // Evaluate integrals of 4 products N(i)*N(j).
    NN.Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        compproduct(product, N[i], N[j]);
        NN(i,j) = blas.DOT(quadweights.M(), quadweights.A(), product.A());
      }
    }

    Epetra_IntSerialDenseVector neumann_nodes(2);
    Epetra_SerialDenseVector neumann_add(2);
    double l;
    for (i=0; i<neumann.M(); i++) {
      neumann_nodes(0) = neumann(i,0); neumann_nodes(1) = neumann(i,1);
      neumann_add(0) = pcoords(overlapmap.LID(neumann_nodes(0)),0)
                      -pcoords(overlapmap.LID(neumann_nodes(1)),0);
      neumann_add(1) = pcoords(overlapmap.LID(neumann_nodes(0)),1)
                      -pcoords(overlapmap.LID(neumann_nodes(1)),1);
      l = blas.NRM2(neumann_add.M(), neumann_add.A());
      neumann_add.Multiply('N', 'N', 1.0, NN, g, 0.0);
      neumann_add.Scale(l);
      err = b->SumIntoGlobalValues(neumann_nodes, neumann_add); 
      if (err<0) return(err);
    }
  }

  // Get Robin boundary edges.
  Epetra_IntSerialDenseMatrix robin(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==3) {
      robin(j,0) = e(i,0);  robin(j,1) = e(i,1);
      j++;
    }
  }
  robin.Reshape(j,2);
  // Adjust for Robin BC's.
  if (robin.M() != 0) {
    Epetra_BLAS blas;
    Epetra_SerialDenseMatrix quadnodes;
    Epetra_SerialDenseVector quadweights;
    Epetra_SerialDenseMatrix N;
    Epetra_SerialDenseMatrix NN;
    Epetra_SerialDenseMatrix * NNN;
    Epetra_SerialDenseVector product(2);

    // Get quadrature weights and points.
    quadrature(1, 2, quadnodes, quadweights);
    
    // Evaluate nodal basis functions at quadrature points
    // N(i,j) value of basis function j at quadrature node i
    N.Shape(quadnodes.M(),2);
    for (i=0; i<quadnodes.M(); i++) {
      N(i,0) = 1.0 - quadnodes(i,0);
      N(i,1) = quadnodes(i,0);
    }

    // Evaluate integrals of 4 products N(i)*N(j).
    NN.Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        compproduct(product, N[i], N[j]);
        NN(i,j) = blas.DOT(quadweights.M(), quadweights.A(), product.A());
      }
    }

    // Evaluate integrals of 8 products N(i)*N(j)*N(k).
    NNN = new Epetra_SerialDenseMatrix[2];
    NNN[0].Shape(2,2); NNN[1].Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          compproduct(product, N[i], N[j], N[k]);
          NNN[k](i,j) = blas.DOT(quadweights.M(), quadweights.A(),
                                 product.A());
        }
      }
    }

    Epetra_IntSerialDenseVector robin_nodes(2);
    Epetra_SerialDenseVector robin_b_add(2);
    Epetra_SerialDenseMatrix robin_A_add(2,2);
    double l;
    for (i=0; i<robin.M(); i++) {
      robin_nodes(0) = robin(i,0); robin_nodes(1) = robin(i,1);
      
      robin_b_add(0) = pcoords(overlapmap.LID(robin_nodes(0)),0)
                      -pcoords(overlapmap.LID(robin_nodes(1)),0);
      robin_b_add(1) = pcoords(overlapmap.LID(robin_nodes(0)),1)
                      -pcoords(overlapmap.LID(robin_nodes(1)),1);
      l = blas.NRM2(robin_b_add.M(), robin_b_add.A());
      robin_b_add.Multiply('N', 'N', 1.0, NN, g, 0.0);
      robin_b_add.Scale(l);
      err = b->SumIntoGlobalValues(robin_nodes, robin_b_add); 
      if (err<0) return(err);

      NNN[0].Scale(sigma(0)); NNN[1].Scale(sigma(1));
      robin_A_add += NNN[0]; robin_A_add += NNN[1];
      robin_A_add.Scale(l);
      err = A->InsertGlobalValues(robin_nodes, robin_A_add, format);
      if (err<0) return(err);
    }

    delete [] NNN;
  }

  // Get Dirichlet boundary edges.
  Epetra_IntSerialDenseMatrix dirichlet(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==1) {
      dirichlet(j,0) = e(i,0);  dirichlet(j,1) = e(i,1);
      j++;
    }
  }
  dirichlet.Reshape(j,2);
  // DIRICHLET NOT DONE! DO THIS LATER!!!!

  /* ************************  Done building A and b.  ************************ */



  /* ************************  Second, we build M.  ************************ */

  Epetra_SerialDenseMatrix Mt;

  for (i=0; i< numNodesPerElem; i++) k(i)=0.0;
  for (i=0; i< numNodesPerElem; i++) { c(i,0)=0.0; c(i,1)=0.0; }
  for (i=0; i< numNodesPerElem; i++) r(i)=1.0;
  for (i=0; i< numNodesPerElem; i++) f(i)=0.0;
  g(0)=0.0; g(1)=0.0;
  sigma(0)=0.0; sigma(1)=0.0;
  for(i=0; i<numLocElems; i++) {
    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }
    lassembly(vertices, k, c, r, f, usr_par, Mt, bt);
    err = M->InsertGlobalValues(epetra_nodes, Mt, format);
    if (err<0) return(err);
  }

  /* ************************  Done building M.  ************************ */



  // Call global assemble and FillComplete at the same time.

  err = A->GlobalAssemble();
  if (err<0) return(err);

  err = b->GlobalAssemble();
  if (err<0) return(err);

  err = M->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}

/* \brief Computes local stiffness matrix and local RHS vector for simplex
           (triangles in two dimensions).
                            
  \param  vertices   [in]  - Matrix containing the global coordinates of the current simplex.
  \param  k          [in]  - Vector containing the value of the diffusion \f$k\f$ at each vertex
                             (\f$k\f$ is assumed to be piecewise linear), where \f$k\f$ is the coeff in the
                             term \f$ \nabla \cdot (k \nabla(u)) \f$ in the advection-diffusion equation.
  \param  c          [in]  - Matrix containing the value of the advection field \f$ \mathbf{c} \f$ at each
                             vertex (\f$ \mathbf{c} \f$ is assumed to be piecewise linear), where
                             \f$ \mathbf{c} \f$ is the 2-vector of coeffs in the term
                             \f$ \mathbf{c}(x) \cdot \nabla y(x) \f$ in the advection-diffusion equation.
  \param  r          [in]  - Vector containing the value of \f$ r \f$ at each vertex (\f$ r \f$ is assumed
                             to be piecewise linear), where \f$ r \f$ is the coefficient in the term
                             \f$ ru \f$ in the advection-diffusion equation.
  \param  f          [in]  - Vector containing the value of \f$ f \f$ at each vertex (\f$ f \f$ is assumed to be
                             piecewise linear), where \f$ f \f$ is the right hand side in the adv-diff eq
  \param  usr_par    [in]  - class containing: \n
                              - S1, S2, S3 (3x3) the integrals of various combinations of partials
                                of local basis functions
                              - N (1x3) integrals of local basis functions
                              - NNN[3] (3x3), integrals of products of local basis functions N(i)*N(j)*N(k)
                              - etc.
  \param  At         [out] - stiffness matrix contribution for the simplex
  \param  bt         [out] - right-hand side contribution for the simplex

  \return 0                 if successful.

  \par Detailed Description:

     Computes the local (per simplex) contributions to the FE volume stiffness matrix \e A and the
     right-hand side \e b for the advection-diffusion equation
    \f{align*}
       - \nabla \cdot \left( K(x) \nabla y(x) \right) + \mathbf{c}(x) \cdot \nabla y(x) + r(x) y(x)
       &= f(x), & x &\in \Omega,\;  \\
       y(x) &= d(x),    & x &\in {\partial \Omega}_d, \\
       K(x) \frac{\partial}{\partial \mathbf{n}}  y(x) &= g(x), & x &\in {\partial \Omega}_n, \\
       \sigma_0 K(x) \frac{\partial}{\partial \mathbf{n}} y(x)
       + \sigma_1 y(x) &= g(x), & x &\in {\partial \Omega}_r.
    \f}
*/
int GLpApp::lassembly(const Epetra_SerialDenseMatrix & vertices,
              const Epetra_SerialDenseVector & k,
              const Epetra_SerialDenseMatrix & c,
              const Epetra_SerialDenseVector & r,
              const Epetra_SerialDenseVector & f,
              const Usr_Par & usr_par,
              Epetra_SerialDenseMatrix & At,
              Epetra_SerialDenseVector & bt)
{
  // Note that the constructors below initialize all entries to 0.
  Epetra_SerialDenseMatrix B(2,2);
  Epetra_SerialDenseMatrix b(2,2);
  Epetra_SerialDenseMatrix BtB(2,2);  
  Epetra_SerialDenseMatrix C(2,2);  
  Epetra_SerialDenseMatrix M1(3,3);
  Epetra_SerialDenseMatrix M2(3,3);
  Epetra_SerialDenseMatrix M3(3,3);
  Epetra_SerialDenseMatrix tmp(3,3);
  double dB, adB;
  At.Shape(3,3);
  bt.Size(3);

  // Construct affine transformation matrix.
  for(int i=0; i<2; i++) {
    B(i,0) = vertices(1,i)-vertices(0,i);
    B(i,1) = vertices(2,i)-vertices(0,i);
  }
  dB  = determinant(B);
  adB = abs(dB);

  // Compute matrix C = inv(B'*B).
  BtB.Multiply('T', 'N', 1.0, B, B, 0.0);
  inverse(BtB, C);

  inverse(B, b); b.Scale(dB);
  
  // Compute the part corresponding to div(K*grad(u)).
  tmp = usr_par.S1; tmp.Scale(C(0,0));
  M1 += tmp;
  tmp = usr_par.S2; tmp.Scale(C(0,1));
  M1 += tmp;
  tmp = usr_par.S3; tmp.Scale(C(1,1));
  M1 += tmp;
  M1.Scale( (k(0)*usr_par.Nw(0) + k(1)*usr_par.Nw(1) +
             k(2)*usr_par.Nw(2)) * adB );

  // Compute the part corresponding to c'*grad(u).
  tmp = usr_par.NdNdx1Nw[0]; tmp.Scale(b(0,0)*c(0,0)+b(0,1)*c(0,1));
  M2 += tmp;
  tmp = usr_par.NdNdx1Nw[1]; tmp.Scale(b(0,0)*c(1,0)+b(0,1)*c(1,1));
  M2 += tmp;
  tmp = usr_par.NdNdx1Nw[2]; tmp.Scale(b(0,0)*c(2,0)+b(0,1)*c(2,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[0]; tmp.Scale(b(1,0)*c(0,0)+b(1,1)*c(0,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[1]; tmp.Scale(b(1,0)*c(1,0)+b(1,1)*c(1,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[2]; tmp.Scale(b(1,0)*c(2,0)+b(1,1)*c(2,1));
  M2 += tmp;
  M2.Scale(adB/dB);

  // Compute the part corresponding to r*u.
  tmp = usr_par.NNNw[0]; tmp.Scale(r(0));
  M3 += tmp;
  tmp = usr_par.NNNw[1]; tmp.Scale(r(1));
  M3 += tmp;
  tmp = usr_par.NNNw[2]; tmp.Scale(r(2));
  M3 += tmp;
  M3.Scale(adB);

  At += M1;
  At += M2;
  At += M3;

  bt.Multiply('T', 'N', adB, usr_par.NNw, f, 0.0);  
  
  return(0);
}

/*  \brief Computes the inverse of a dense matrix.

  \param  mat  [in]  - the matrix that is to be inverted
  \param  inv  [in]  - the inverse

  \return 0            if successful
*/
int GLpApp::inverse(const Epetra_SerialDenseMatrix & mat,
            Epetra_SerialDenseMatrix & inv)
{
  Epetra_LAPACK lapack;
  int dim = mat.M();
  int info;
  Epetra_IntSerialDenseVector ipiv(dim);
  Epetra_SerialDenseVector work(2*dim);

  inv.Shape(dim,dim);
  inv = mat;

  lapack.GETRF(dim, dim, inv.A(), dim, ipiv.A(), &info);
  lapack.GETRI(dim, inv.A(), dim, ipiv.A(), work.A(), &dim, &info);
  
  return(0);
}

/*  \brief Computes the determinant of a dense matrix.

  \param  mat  [in]  - the matrix

  \return the determinant
*/
double GLpApp::determinant(const Epetra_SerialDenseMatrix & mat)
{
  //Teuchos::LAPACK<int,double> lapack;
  Epetra_LAPACK lapack;
  double det;
  int swaps; 
  int dim = mat.M();
  int info;
  Epetra_IntSerialDenseVector ipiv(dim);
 
  Epetra_SerialDenseMatrix mymat(mat);
  lapack.GETRF(dim, dim, mymat.A(), dim, ipiv.A(), &info);

  det = 1.0;
  for (int i=0; i<dim; i++) {
    det *= mymat(i,i);
  }

  swaps = 0;
  for (int i=0; i<dim; i++) {
    if ((ipiv[i]-1) != i)
      swaps++;
  }

  det *= pow((double)(-1.0),swaps);

  return(det);
}

int GLpApp::meshreader(const Epetra_Comm & Comm,
               Epetra_IntSerialDenseVector & ipindx,
               Epetra_SerialDenseMatrix & ipcoords,
               Epetra_IntSerialDenseVector & pindx,
               Epetra_SerialDenseMatrix & pcoords,
               Epetra_IntSerialDenseMatrix & t,
               Epetra_IntSerialDenseMatrix & e,
               const char geomFileBase[],
               const bool trace,
               const bool dumpAll
               )
{
  int MyPID = Comm.MyPID();

  const int FileNameSize = 120;
  char FileName[FileNameSize];
  TEST_FOR_EXCEPT(static_cast<int>(std::strlen(geomFileBase) + 5) > FileNameSize);
  sprintf(FileName, "%s.%03d", geomFileBase, MyPID);

  {
    std::ifstream file_in(FileName);
    TEST_FOR_EXCEPTION(
      file_in.eof(), std::logic_error
      ,"Error, the file \""<<FileName<<"\" could not be opened for input!"
      );
  }

  FILE* fid;
  fid = fopen(FileName, "r");

  if(trace) printf("\nReading node info from %s ...\n", FileName);
  int numip = 0, numcp = 0; // # owned nodes, # shared nodes
  fscanf(fid, "%d %d", &numip, &numcp);
  const int nump = numip + numcp; // # total nodes
  if(trace) printf( "\nnumip = %d, numcp = %d, nump = %d\n", numip, numcp, nump );
  ipindx.Size(numip);
  ipcoords.Shape(numip, 2);
  pindx.Size(nump);
  pcoords.Shape(nump, 2);
  for (int i=0; i<numip; i++) {
    fscanf(fid,"%d %lf %lf %*d",&ipindx(i),&ipcoords(i,0),&ipcoords(i,1));
    if(trace&&dumpAll) printf("%d %lf %lf\n",ipindx(i),ipcoords(i,0),ipcoords(i,1));
    pindx(i) = ipindx(i);
    pcoords(i,0) = ipcoords(i,0); pcoords(i,1) = ipcoords(i,1);
  }
  for (int i=numip; i<nump; i++) {
    fscanf(fid,"%d %lf %lf %*d",&pindx(i),&pcoords(i,0),&pcoords(i,1));
  }

  fscanf(fid,"%*[^\n]");   // Skip to the End of the Line
  fscanf(fid,"%*1[\n]");   // Skip One Newline

  fscanf(fid,"%*[^\n]");   // Skip to the End of the Line
  fscanf(fid,"%*1[\n]");   // Skip One Newline

  for (int i=0; i<nump; i++) {
    fscanf(fid,"%*[^\n]"); // Skip to the End of the Line
    fscanf(fid,"%*1[\n]"); // Skip One Newline
  }

  if(trace) printf("\nReading element info from %s ...\n", FileName);
  int numelems = 0; // # elements on this processor
  fscanf(fid, "%d", &numelems);
  if(trace) printf( "\nnumelems = %d\n", numelems );
  t.Shape(numelems,3);
  for (int i=0; i<numelems; i++) {
    fscanf(fid, "%d %d %d", &t(i,0), &t(i,1), &t(i,2));
    if(trace&&dumpAll) printf("%d %d %d\n", t(i,0), t(i,1), t(i,2));
  }

  if(trace) printf("\nReading boundry edge info from %s ...\n", FileName);
  int numedges = 0; // # boundry edges on this processor
  fscanf(fid,"%d",&numedges);
  if(trace) printf( "\nnumedges = %d\n", numedges );
  e.Shape(numedges,3);
  for (int i=0; i<numedges; i++) {
    fscanf(fid, "%d %d %d", &e(i,0), &e(i,1), &e(i,2));
    if(trace&&dumpAll) printf("%d %d %d\n", e(i,0), e(i,1), e(i,2));
  }

  fclose(fid);
  if(trace) printf("Done reading mesh.\n");

  return(0);

}

/*  \brief Performs finite-element assembly of the Hessian of the nonlinear form times an arbitrary vector.

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the j-th vertex in triangle i, where i = 1, ..., ELE
  \param  y         [in]  - Reference-counting pointer to the Epetra_MultiVector at which the nonlinear
                            term is evaluated.
  \param  s         [in]  - Reference-counting pointer to the Epetra_MultiVector that is multiplied by the
                            Hessian of the nonlinear form.
  \param  lambda    [in]  - Reference-counting pointer to the Epetra_MultiVector of Lagrange Multipliers.
  \param  hessvec   [out] - Reference-counting pointer to the Epetra_FECrsMatrix containing Hessian of
                            the nonlinear form times vector product.
  \return 0                 if successful.

  \par Detailed Description:

  Assembles the nonlinear term \e hessvec, represented by
     \f[
       \{N''(y,\lambda)s\}_{j} = \langle g''(y_h) \lambda s,\phi_j \rangle
                               = \int_{\Omega} g''(y_h(x)) \lambda(x) s(x) \phi_j(x) dx,
     \f]
  where \f$ g''(y_h) \f$ is given in the local function \e g2pfctn, and \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the
  piecewise linear nodal basis for the state space.
*/
int GLpApp::nonlinhessvec(const Epetra_Comm & Comm,
                  const Epetra_IntSerialDenseVector & ipindx,
                  const Epetra_SerialDenseMatrix & ipcoords,
                  const Epetra_IntSerialDenseVector & pindx,
                  const Epetra_SerialDenseMatrix & pcoords,
                  const Epetra_IntSerialDenseMatrix & t,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & y,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & s,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & lambda,
                  Teuchos::RefCountPtr<Epetra_FEVector> & hessvec)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int numLocNodes     = pindx.M();
  int numMyLocNodes   = ipindx.M();
  int numLocElems     = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, numMyLocNodes, (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, numLocNodes, (int*)pindx.A(), indexBase, Comm);

  hessvec = rcp(new Epetra_FEVector(standardmap,1));

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;
  
  // get quadrature nodes and weights
  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;
  quadrature(2,3,Nodes,Weights);
  int numQuadPts = Nodes.M();

  // Evaluate nodal basis functions and their derivatives at quadrature points
  // N(i,j) = value of the j-th basis function at quadrature node i.
  Epetra_SerialDenseMatrix N;
  N.Shape(numQuadPts,3);
  for (int i=0; i<numQuadPts; i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }

  // Declare quantities needed for the call to the local assembly routine.
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());

  Epetra_SerialDenseVector ly;        // local entries of y
  Epetra_SerialDenseVector Nly;       // N*ly
  Epetra_SerialDenseVector ls;        // local entries of s
  Epetra_SerialDenseVector Nls;       // N*ls
  Epetra_SerialDenseVector llambda;   // local entries of lambda
  Epetra_SerialDenseVector Nllambda;  // N*llambda
  Epetra_SerialDenseVector lgfctn;    // gfctn(Nly)
  Epetra_SerialDenseVector lgfctnNi;  // lgfctn.*N(:,i)
  Epetra_SerialDenseVector lgfctnNls; // lgfctn.*N(:,i).*llambda.*ls
  Epetra_SerialDenseVector lhessvec;  // local contribution
  // Size and init to zero.
  ly.Size(numNodesPerElem);
  Nly.Size(numQuadPts);
  ls.Size(numNodesPerElem);
  Nls.Size(numQuadPts);
  llambda.Size(numNodesPerElem);
  Nllambda.Size(numQuadPts);
  lgfctn.Size(numQuadPts);
  lgfctnNi.Size(numQuadPts);
  lgfctnNls.Size(numQuadPts);
  lhessvec.Size(numNodesPerElem);
  
  Epetra_SerialDenseMatrix B(2,2);
  double adB;
  
  for(i=0; i<numLocElems; i++) {

    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }

    // Construct affine transformation matrix.
    for(int i=0; i<2; i++) {
      B(i,0) = vertices(1,i)-vertices(0,i);
      B(i,1) = vertices(2,i)-vertices(0,i);
    }
    adB  = abs(determinant(B));

    // Construct local (to each processor) element view of y, s, l.
    for (j=0; j<numNodesPerElem; j++) {
      ly(j) = (*((*y)(0)))[overlapmap.LID(nodes[j])];
      ls(j) = (*((*s)(0)))[overlapmap.LID(nodes[j])];
      llambda(j) = (*((*lambda)(0)))[overlapmap.LID(nodes[j])];
    }

    Nly.Multiply('N', 'N', 1.0, N, ly, 0.0);            //N*ly
    Nls.Multiply('N', 'N', 1.0, N, ls, 0.0);            //N*ls
    Nllambda.Multiply('N', 'N', 1.0, N, llambda, 0.0);  //N*llambda
    g2pfctn(Nly, lgfctn);
    
    for (int i=0; i<numNodesPerElem; i++) {
      compproduct(lgfctnNi, lgfctn.A(), N[i]);
      compproduct(lgfctnNls, lgfctnNi.A(), Nls.A(), Nllambda.A());  // g2pfctn(Nly).*Nls.*Nllambda.*N(:,i)
      lhessvec(i) = adB*lgfctnNls.Dot(Weights);
    }

    err = hessvec->SumIntoGlobalValues(epetra_nodes, lhessvec);
    if (err<0) return(err);
  }

  // Call global assemble.

  err = hessvec->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}

/*  \brief Componentwise evaluation of the second derivative of the nonlinear reaction term.
  \param  v   [in]  - Vector at which the second derivative is evaluated.
  \param  gv  [out] - Vector value.
*/
void GLpApp::g2pfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv) {
  for (int i=0; i<v.M(); i++) {
    gv(i) = 6.0*v(i);
  }  
}

/*  \brief Performs finite-element assembly of the Jacobian of the nonlinear form.

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the j-th vertex in triangle i, where i = 1, ..., ELE
  \param  y         [in]  - Reference-counting pointer to the Epetra_MultiVector at which the nonlinear
                            term is evaluated.
  \param  Gp        [out] - Reference-counting pointer to the Epetra_FECrsMatrix containing the Jacobian
                            of the nonlinear form.
  \return 0                 if successful.

  \par Detailed Description:

  Assembles the nonlinear term \e Gp, represented by
     \f[
       \{N'(y)\}_{jk} = \langle g'(y_h) \phi_k,\phi_j \rangle =  \int_{\Omega} g'(y_h(x)) \phi_k(x) \phi_j(x) dx,
     \f]
  where \f$ g'(y_h) \f$ is given in the local function \e gpfctn, and \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the
  piecewise linear nodal basis for the state space.
*/
int GLpApp::nonlinjac(const Epetra_Comm & Comm,
              const Epetra_IntSerialDenseVector & ipindx,
              const Epetra_SerialDenseMatrix & ipcoords,
              const Epetra_IntSerialDenseVector & pindx,
              const Epetra_SerialDenseMatrix & pcoords,
              const Epetra_IntSerialDenseMatrix & t,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> & y,
              Teuchos::RefCountPtr<Epetra_FECrsMatrix> & Gp)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int numLocNodes     = pindx.M();
  int numMyLocNodes   = ipindx.M();
  int numLocElems     = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, numMyLocNodes, (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, numLocNodes, (int*)pindx.A(), indexBase, Comm);

  int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
  Gp = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;
  
  // get quadrature nodes and weights
  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;
  quadrature(2,3,Nodes,Weights);
  int numQuadPts = Nodes.M();

  // Evaluate nodal basis functions and their derivatives at quadrature points
  // N(i,j) = value of the j-th basis function at quadrature node i.
  Epetra_SerialDenseMatrix N;
  N.Shape(numQuadPts,3);
  for (int i=0; i<numQuadPts; i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }
  
  // Declare quantities needed for the call to the local assembly routine.
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());

  Epetra_SerialDenseVector ly;          // local entries of y
  Epetra_SerialDenseVector Nly;         // N*ly
  Epetra_SerialDenseVector lgfctn;      // gfctn(Nly)
  Epetra_SerialDenseVector lgfctnNiNj;  // lgfctn.*N(:,i).*N(:,j)
  Epetra_SerialDenseMatrix lGp;         // local contribution
  // Size and init to zero.
  ly.Size(numNodesPerElem);
  Nly.Size(numQuadPts);
  lgfctn.Size(numQuadPts);
  lgfctnNiNj.Size(numQuadPts);
  lGp.Shape(numNodesPerElem, numNodesPerElem);
  
  Epetra_SerialDenseMatrix B(2,2);
  double adB;
  
  for(i=0; i<numLocElems; i++) {

    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }

    // Construct affine transformation matrix.
    for(int i=0; i<2; i++) {
      B(i,0) = vertices(1,i)-vertices(0,i);
      B(i,1) = vertices(2,i)-vertices(0,i);
    }
    adB  = abs(determinant(B));

    // Construct local (to each processor) element view of y. 
    for (j=0; j<numNodesPerElem; j++) {
      ly(j) = (*((*y)(0)))[overlapmap.LID(nodes[j])];
    }

    Nly.Multiply('N', 'N', 1.0, N, ly, 0.0);
    gpfctn(Nly, lgfctn);
    
    for (int i=0; i<numNodesPerElem; i++) {
      for (int j=0; j<numNodesPerElem; j++) {
        compproduct(lgfctnNiNj, lgfctn.A(), N[i], N[j]);
        lGp(i,j) = adB*lgfctnNiNj.Dot(Weights);
      }
    }
    
    err = Gp->InsertGlobalValues(epetra_nodes, lGp, format);
    if (err<0) return(err);
  }

  // Call global assemble.

  err = Gp->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}

/*  \brief Componentwise evaluation of the first derivative of the nonlinear reaction term.
  \param  v   [in]  - Vector at which the first derivative is evaluated.
  \param  gv  [out] - Vector value.
*/
void GLpApp::gpfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv) {
  for (int i=0; i<v.M(); i++) {
    gv(i) = 3.0*pow(v(i),2)-1.0;
  }  
}

/*  \brief Performs finite-element assembly of the nonlinear reaction term.

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the j-th vertex in triangle i, where i = 1, ..., ELE
  \param  y         [out] - Reference-counting pointer to the Epetra_MultiVector at which the nonlinear
                            form is evaluated.
  \param  g         [out] - Reference-counting pointer to the Epetra_FEVector containing the value
                            of the nonlinear form.
  \return 0                 if successful.

  \par Detailed Description:

  Assembles the nonlinear term \e g, represented by
     \f[
       \{N(y)\}_{j} = \langle g(y_h),\phi_j \rangle =  \int_{\Omega} g(y_h(x)) \phi_j(x) dx,
     \f]
  where \f$ g(y_h) \f$ is given in the local function \e gfctn, and \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the
  piecewise linear nodal basis for the state space.
*/
int GLpApp::nonlinvec(const Epetra_Comm & Comm,
              const Epetra_IntSerialDenseVector & ipindx,
              const Epetra_SerialDenseMatrix & ipcoords,
              const Epetra_IntSerialDenseVector & pindx,
              const Epetra_SerialDenseMatrix & pcoords,
              const Epetra_IntSerialDenseMatrix & t,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> & y,
              Teuchos::RefCountPtr<Epetra_FEVector> & g)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int numLocNodes     = pindx.M();
  int numMyLocNodes   = ipindx.M();
  int numLocElems     = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, numMyLocNodes, (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, numLocNodes, (int*)pindx.A(), indexBase, Comm);

  g = rcp(new Epetra_FEVector(standardmap,1));

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;
  
  // get quadrature nodes and weights
  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;
  quadrature(2,3,Nodes,Weights);
  int numQuadPts = Nodes.M();

  // Evaluate nodal basis functions and their derivatives at quadrature points
  // N(i,j) = value of the j-th basis function at quadrature node i.
  Epetra_SerialDenseMatrix N;
  N.Shape(numQuadPts,3);
  for (int i=0; i<numQuadPts; i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }

  // Declare quantities needed for the call to the local assembly routine.
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());

  Epetra_SerialDenseVector ly;        // local entries of y
  Epetra_SerialDenseVector Nly;       // N*ly
  Epetra_SerialDenseVector lgfctn;    // gfctn(Nly)
  Epetra_SerialDenseVector lgfctnNi;  // lgfctn.*N(:,i)
  Epetra_SerialDenseVector lg;        // local contribution
  // Size and init to zero.
  ly.Size(numNodesPerElem);
  Nly.Size(numQuadPts);
  lgfctn.Size(numQuadPts);
  lgfctnNi.Size(numQuadPts);
  lg.Size(numNodesPerElem);
  
  Epetra_SerialDenseMatrix B(2,2);
  double adB;
  
  for(i=0; i<numLocElems; i++) {

    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }

    // Construct affine transformation matrix.
    for(int i=0; i<2; i++) {
      B(i,0) = vertices(1,i)-vertices(0,i);
      B(i,1) = vertices(2,i)-vertices(0,i);
    }
    adB  = abs(determinant(B));

    // Construct local (to each processor) element view of y. 
    for (j=0; j<numNodesPerElem; j++) {
      ly(j) = (*((*y)(0)))[overlapmap.LID(nodes[j])];
    }

    Nly.Multiply('N', 'N', 1.0, N, ly, 0.0);
    gfctn(Nly, lgfctn);
    
    for (int i=0; i<numNodesPerElem; i++) {
      compproduct(lgfctnNi, lgfctn.A(), N[i]);
      lg(i) = adB*lgfctnNi.Dot(Weights);
    }
    
    err = g->SumIntoGlobalValues(epetra_nodes, lg);
    if (err<0) return(err);
  }

  // Call global assemble.

  err = g->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}


/*  \brief Componentwise evaluation of the nonlinear reaction term.
  \param  v   [in]  - Vector at which the nonlinear function is evaluated.
  \param  gv  [out] - Vector value.
*/
void GLpApp::gfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv) {
  for (int i=0; i<v.M(); i++) {
    gv(i) = pow(v(i),3)-v(i);
  }  
}

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the distributed CrsMatrix to 
 *                     print out
 * - ostream &         reference to output stream 
 */

bool GLpApp::CrsMatrix2MATLAB(const Epetra_CrsMatrix & A, ostream & outfile) 
{

  int MyPID = A.Comm().MyPID(); 
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) { 
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.TransformToLoca()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return false;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    outfile << "A = spalloc(";
    outfile << NumGlobalRows << ',' << NumGlobalRows;
    outfile << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    A.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "\n\n% On proc " << Proc << ": ";
      outfile << NumMyRows << " rows and ";
      outfile << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = A.GRID(MyRow);

        NumNzRow = A.NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        A.ExtractMyRowCopy(MyRow, NumNzRow, 
                           NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          outfile << "A(" << GlobalRow  + IndexBase 
                  << "," << A.GCID(Indices[j]) + IndexBase
                  << ") = " << Values[j] << ";\n";
        }

        delete Values;
        delete Indices;
      }
      
    }
    A.Comm().Barrier();
  }

  return true;

}


/* ======== ============= *
 * function Vector2MATLAB *
 * ======== ============= *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_Vector     reference to vector
 * - ostream &         reference to output stream 
 */

bool GLpApp::Vector2MATLAB( const Epetra_Vector & v, ostream & outfile)
{
  
  int MyPID = v.Comm().MyPID(); 
  int NumProc = v.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();
  
  // print out on cout if no filename is provided

  // write on file the dimension of the matrix

  if( MyPID == 0 ) outfile << "v = zeros(" << GlobalLength << ",1)\n";

  int NumMyElements = v.Map().NumMyElements();
  // get update list
  int * MyGlobalElements = v.Map().MyGlobalElements( );
  
  int Row;

  int IndexBase = v.Map().IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    v.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "% On proc " << Proc << ": ";
      outfile << MyLength << " rows of ";
      outfile << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        outfile << "v(" << MyGlobalElements[Row] + IndexBase
             << ") = " << v[Row] << ";\n";
      }
      
    }
      
    v.Comm().Barrier();
  }

  return true;

} /* Vector2MATLAB */


/* ======== =============== *
 * function FEVector2MATLAB *
 * ======== =============== *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_FEVector   reference to FE vector
 * - ostream &         reference to output stream 
 */

bool GLpApp::FEVector2MATLAB( const Epetra_FEVector & v, ostream & outfile)
{
  
  int MyPID = v.Comm().MyPID(); 
  int NumProc = v.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();
  
  // print out on cout if no filename is provided

  // write on file the dimension of the matrix

  if( MyPID == 0 ) outfile << "v = zeros(" << GlobalLength << ",1)\n";

  int NumMyElements = v.Map().NumMyElements();
  // get update list
  int * MyGlobalElements = v.Map().MyGlobalElements( );
  
  int Row;

  int IndexBase = v.Map().IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    v.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "% On proc " << Proc << ": ";
      outfile << MyLength << " rows of ";
      outfile << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        outfile << "v(" << MyGlobalElements[Row] + IndexBase
             << ") = " << v[0][Row] << ";\n";
      }
      
    }
      
    v.Comm().Barrier();
  }

  return true;

} /* FEVector2MATLAB */


/*  \brief  Returns the nodes and weights for the integration \n
             on the interval [0,1] (dim = 1) \n
             on the triangle with vertices (0,0), (1,0), (0,1) (if dim = 2) \n
             on the tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1) (if dim = 3).
             
  \param  dim     [in]   - spatial dimension (dim = 1, 2)
  \param  order   [in]   - required degree of polynomials that integrate exactly
  \param  nodes   [out]  - Matrix in which the i-th row of nodes gives coordinates of the i-th quadrature node
  \param  weights [out]  - quadrature weights

  \return 0                if successful
*/
int GLpApp::quadrature(const int dim, const int order,
               Epetra_SerialDenseMatrix & nodes,
               Epetra_SerialDenseVector & weights)
{
  
  if (dim == 1) {

    // Gauss quadrature nodes and weights on the interval [0,1]

    if (order == 1) {
      nodes.Shape(1,1);
      nodes(0,0) = 0.5;
      weights.Size(1);
      weights(0) = 1.0;
    }
    else if (order == 2) {
      nodes.Shape(2,1);
      nodes(0,0) = (1.0-1.0/sqrt(3.0))/2.0;
      nodes(1,0) = (1.0+1.0/sqrt(3.0))/2.0;
      weights.Size(2);
      weights(0) = 0.5;
      weights(1) = 0.5;
    }
    else if (order == 3) {
      nodes.Shape(3,1);
      nodes(0,0) = (1.0-sqrt(3.0/5.0))/2.0;
      nodes(1,0) = 0.5;
      nodes(2,0) = (1.0+sqrt(3.0/5.0))/2.0;
      weights.Size(3);
      weights(0) = 5.0/18.0;
      weights(1) = 4.0/9.0;
      weights(2) = 5.0/18.0;
    }
    else {
      cout << "Quadrature for dim = " << dim << " and order = ";
      cout << order << " not available.\n";
      exit(-1);
    }

  }
  else if (dim == 2) {
    
    // Quadrature nodes and weights on the unit simplex with
    // vertices (0,0), (1,0), and (0,1).

    if (order == 1) {
      nodes.Shape(1,2);
      nodes(0,0) = 1.0/3.0; nodes (0,1) = 1.0/3.0;
      weights.Size(1);
      weights(0) = 0.5;
    }
    else if (order == 2) {
      nodes.Shape(3,2);
      nodes(0,0) = 1.0/6.0; nodes (0,1) = 1.0/6.0;
      nodes(1,0) = 2.0/3.0; nodes (1,1) = 1.0/6.0;
      nodes(2,0) = 1.0/6.0; nodes (2,1) = 2.0/3.0;
      weights.Size(3);
      weights(0) = 1.0/6.0;
      weights(1) = 1.0/6.0;
      weights(2) = 1.0/6.0;
    }
    else if (order == 3) {
      nodes.Shape(4,2);
      nodes(0,0) = 1.0/3.0; nodes (0,1) = 1.0/3.0;
      nodes(1,0) = 3.0/5.0; nodes (1,1) = 1.0/5.0;
      nodes(2,0) = 1.0/5.0; nodes (2,1) = 3.0/5.0;
      nodes(3,0) = 1.0/5.0; nodes (3,1) = 1.0/5.0;
      weights.Size(4);
      weights(0) = -9.0/32.0;
      weights(1) = 25.0/96.0;
      weights(2) = 25.0/96.0;
      weights(3) = 25.0/96.0;
    }
    else {
      cout << "Quadrature for dim = " << dim << " and order = ";
      cout << order << " not available.\n";
      exit(-1);
    }

  }
  else {
    cout << "Quadrature for dim = " << dim << " not available.\n";
    exit(-1);
  }

  return(0);
}
