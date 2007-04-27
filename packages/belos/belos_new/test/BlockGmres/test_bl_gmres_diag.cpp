// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test generates diagonal matrices for block GMRES to solve.  
//
// NOTE: No preconditioner is used in this case.
//

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosEpetraOperator.h>
#include <BelosBlockGmresSolMgr.hpp>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_Time.hpp>

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

using std::vector;
using namespace Belos;

//************************************************************************************************

class Vector_Operator
{
public:
  
  Vector_Operator(int m, int n) : m(m), n(n) {};
  
  virtual ~Vector_Operator() {};
  
  virtual void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y) = 0;
  
  int size (int dim) const { return (dim == 1) ? m : n; };
  
protected:
  
  int m, n;        // an (m x n) operator 
  
private:
  
  // Not allowing copy construction.
  Vector_Operator( const Vector_Operator& ) {};
  Vector_Operator* operator=( const Vector_Operator& ) { return NULL; };
  
};

//************************************************************************************************

class Diagonal_Operator : public Vector_Operator
{
public:
  
  Diagonal_Operator(int n, double v) : Vector_Operator(n, n), v(v) { };
  
  ~Diagonal_Operator() { };
  
  void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y)
  {
    y.Scale( v, x );
  };
  
private:
  
  double v;
};

//************************************************************************************************

class Diagonal_Operator_2 : public Vector_Operator
{
public:
  
  Diagonal_Operator_2(int n, double v) : Vector_Operator(n, n), v(v) { };
  
  ~Diagonal_Operator_2() { };
  
  void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y)
  {
    int myCols = y.MyLength();
    for (int j=0; j < x.NumVectors(); ++j) {
      for (int i=0; i < myCols; ++i) (*y(j))[i] = (i+1)*v*(*x(j))[i];  // NOTE: square operator!
    }
  };
  
private:
  
  double v;
};

//************************************************************************************************

class Composed_Operator : public Vector_Operator
{
public:
  
  Composed_Operator(int n, 
		    const Teuchos::RefCountPtr<Vector_Operator>& pA, 
		    const Teuchos::RefCountPtr<Vector_Operator>& pB);
  
  virtual ~Composed_Operator() {};
  
  virtual void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y);
  
private:
  
  Teuchos::RefCountPtr<Vector_Operator> pA; 
  Teuchos::RefCountPtr<Vector_Operator> pB; 
};

Composed_Operator::Composed_Operator(int n, 
                                     const Teuchos::RefCountPtr<Vector_Operator>& pA, 
                                     const Teuchos::RefCountPtr<Vector_Operator>& pB) 
  : Vector_Operator(n, n), pA(pA), pB(pB) 
{
}

void Composed_Operator::operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y)
{
  Epetra_MultiVector ytemp(y.Map(), y.NumVectors(), false);
  (*pB)( x, ytemp );
  (*pA)( ytemp, y );
}

//************************************************************************************************

class Trilinos_Interface : public Epetra_Operator
{
public:
  
  Trilinos_Interface(const Teuchos::RefCountPtr<Vector_Operator>   pA,
		     const Teuchos::RefCountPtr<const Epetra_Comm> pComm,
		     const Teuchos::RefCountPtr<const Epetra_Map>  pMap)
    : pA(pA), pComm(pComm), pMap(pMap) 
  {
  };
  
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  {
    return(Apply(X,Y));  // No inverse
  };
  
  virtual ~Trilinos_Interface() {};
  
  const char * Label() const {return("Trilinos_Interface, an Epetra_Operator implementation");}; 
  
  bool UseTranspose() const {return(use_transpose);};      // always set to false
  
  int SetUseTranspose(bool UseTranspose_) { use_transpose = false; return(-1); };
  
  bool HasNormInf() const {return(false);};                // cannot return inf-norm
  
  double NormInf() const {return(0.0);};                   
  
  virtual const Epetra_Comm & Comm() const {return *pComm; }
  
  virtual const Epetra_Map & OperatorDomainMap() const {return *pMap; }
  
  virtual const Epetra_Map & OperatorRangeMap() const {return *pMap; }
  
private:
  
  Teuchos::RefCountPtr<Vector_Operator>   pA;
  Teuchos::RefCountPtr<const Epetra_Comm> pComm;
  Teuchos::RefCountPtr<const Epetra_Map>  pMap;
  
  bool use_transpose;
};

int Trilinos_Interface::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  (*pA)(X,Y);
  
  return(0);
}

//************************************************************************************************

class Iterative_Inverse_Operator : public Vector_Operator
{
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator    OP;
  
public:
  
  Iterative_Inverse_Operator(int n, int blocksize, 
			     const Teuchos::RefCountPtr<Vector_Operator>& pA, 
			     string opString="Iterative Solver", bool print=false);              
  
  virtual ~Iterative_Inverse_Operator() {}
  
  virtual void operator () (const Epetra_MultiVector &b, Epetra_MultiVector &x);
  
private:
  
  Teuchos::RefCountPtr<Vector_Operator> pA;       // operator which will be inverted 
  // supplies a matrix vector multiply
  const bool print;
  
  Teuchos::Time timer;
  Teuchos::RefCountPtr<Epetra_Comm> pComm;
  Teuchos::RefCountPtr<Epetra_Map>  pMap;
  
  Teuchos::RefCountPtr<OP> pPE;   
  Teuchos::RefCountPtr<Teuchos::ParameterList>         pList;
  Teuchos::RefCountPtr<LinearProblem<double,MV,OP> >   pProb;
  Teuchos::RefCountPtr<BlockGmresSolMgr<double,MV,OP> >      pBelos;
};

Iterative_Inverse_Operator::Iterative_Inverse_Operator(int n, int blocksize,
                                                       const Teuchos::RefCountPtr<Vector_Operator>& pA, 
                                                       string opString, bool print)
  : Vector_Operator(n, n),      // square operator
    pA(pA), 
    print(print),
    timer(opString)
{
  
  int n_global;
  
#ifdef EPETRA_MPI
  MPI_Allreduce(&n, &n_global, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  pComm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  pComm = Teuchos::rcp( new Epetra_SerialComm() );
  n_global = n;
#endif
  pMap =  Teuchos::rcp( new Epetra_Map(n_global, n, 0, *pComm) );
  
  pPE = Teuchos::rcp( new Trilinos_Interface(pA, pComm, pMap ) );
  
  pProb = Teuchos::rcp( new LinearProblem<double,MV,OP>() );
  pProb->setOperator( pPE );
  
  int restart  = 10; 
  int max_iter = 100;
  double tol = 1.0e-10;
  int verbosity = Belos::Errors + Belos::Warnings;
  if (print)
    verbosity += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
  
  pList = Teuchos::rcp( new Teuchos::ParameterList );
  pList->set( "Num Blocks", n/blocksize );
  pList->set( "Block Size", blocksize );  
  pList->set( "Maximum Iterations", max_iter ); 
  pList->set( "Maximum Restarts", restart );
  pList->set( "Convergence Tolerance", tol ); 
  pList->set( "Verbosity", verbosity );

  pBelos = Teuchos::rcp( new BlockGmresSolMgr<double,MV,OP>(pProb, *pList) );
}

void Iterative_Inverse_Operator::operator () (const Epetra_MultiVector &b, Epetra_MultiVector &x)
{
  // Initialize the solution to zero
  x.PutScalar( 0.0 );

  // Reset the solver, problem, and status test for next solve (HKT)
  pProb->setProblem( Teuchos::rcp(&x, false), Teuchos::rcp(&b, false) );
  
  timer.start();
  Belos::ReturnType ret = pBelos->solve();
  timer.stop();
  
  int pid = pComm->MyPID();
  
  if (pid == 0 && print)
    if (ret == Belos::Converged)
      {
	cout << endl << "pid[" << pid << "] Block GMRES converged" << endl;
	cout << "Solution time: " << timer.totalElapsedTime() << endl;
	
      }
    else 
      cout << endl << "pid[" << pid << "] Block GMRES did not converge" << endl;
}

//************************************************************************************************
//************************************************************************************************

int main(int argc, char *argv[])
{
  
  int pid = -1;
  
#ifdef EPETRA_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  pid = Comm.MyPID();

  int n(10);
  int numRHS=1;
  
  Epetra_Map Map = Epetra_Map(n, 0, Comm);
  Epetra_MultiVector X(Map, numRHS, false), Y(Map, numRHS, false);
  X.PutScalar( 1.0 );

  // Inner computes inv(D2)*y
  Teuchos::RefCountPtr<Diagonal_Operator_2> D2 = Teuchos::rcp(new Diagonal_Operator_2(n, 1.0));
  Iterative_Inverse_Operator A2(n, 1, D2, "Belos (inv(D2))", true); 
  
  // should return x=(1, 1/2, 1/3, ..., 1/10)
  A2(X,Y);
  
  if (pid==0) {
    std::cout << "Vector Y should have all entries [1, 1/2, 1/3, ..., 1/10]" << std::endl;
  }
  Y.Print(std::cout);
  
  // Inner computes inv(D)*x
  Teuchos::RefCountPtr<Diagonal_Operator> D = Teuchos::rcp(new Diagonal_Operator(n, 4.0));
  Teuchos::RefCountPtr<Iterative_Inverse_Operator> Inner = 
    Teuchos::rcp(new Iterative_Inverse_Operator(n, 1, D, "Belos (inv(D))", false)); 
  
  // Composed_Operator computed inv(D)*B*x
  Teuchos::RefCountPtr<Diagonal_Operator> B = Teuchos::rcp(new Diagonal_Operator(n, 4.0));
  Teuchos::RefCountPtr<Composed_Operator> C = Teuchos::rcp(new Composed_Operator(n, Inner, B)); 
  
  // Outer computes inv(C) = inv(inv(D)*B)*x = inv(B)*D*x = x
  Teuchos::RefCountPtr<Iterative_Inverse_Operator> Outer =
    Teuchos::rcp(new Iterative_Inverse_Operator(n, 1, C, "Belos (inv(C)=inv(inv(D)*B))", true)); 
  
  // should return x=1/4
  (*Inner)(X,Y);
  
  if (pid==0) {
    std::cout << std::endl << "Vector Y should have all entries [1/4, 1/4, 1/4, ..., 1/4]" << std::endl;
  }
  Y.Print(std::cout);
  
  // should return x=1
  (*Outer)(X,Y);
  
  if (pid==0) {
    std::cout << "Vector Y should have all entries [1, 1, 1, ..., 1]" << std::endl;
  } 
  Y.Print(std::cout); 
 
  // Compute the norm of Y - 1.0
  std::vector<double> norm_Y(Y.NumVectors());
  Y.Update(-1.0, X, 1.0);
  Y.Norm2(&norm_Y[0]);

  if (pid==0)
    std::cout << "Two-norm of vector (Y-1.0) : "<< norm_Y[0] << std::endl;
  
#ifdef EPETRA_MPI
  MPI_Finalize(); 
#endif
  
  if (norm_Y[0] > 1e-10 || Teuchos::ScalarTraits<double>::isnaninf( norm_Y[0] ) ) {
    if (pid==0)
      cout << "End Result: TEST FAILED" << endl;
    return -1;
  }
  
  //
  // Default return value
  //
  if (pid==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
}
