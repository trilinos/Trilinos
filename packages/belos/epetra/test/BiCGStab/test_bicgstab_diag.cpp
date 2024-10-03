// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include <BelosBiCGStabSolMgr.hpp>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_Time.hpp>
#include "Teuchos_StandardCatchMacros.hpp"

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

    Vector_Operator(int m_in, int n_in) : m(m_in), n(n_in) {};

    virtual ~Vector_Operator() {};

    virtual void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y) = 0;

    int size (int dim) const { return (dim == 1) ? m : n; };

  protected:

    int m, n;        // an (m x n) operator

  private:

    // Not allowing copy construction.
    Vector_Operator( const Vector_Operator& ): m(0), n(0) {};
    Vector_Operator* operator=( const Vector_Operator& ) { return NULL; };

};

//************************************************************************************************

class Diagonal_Operator : public Vector_Operator
{
  public:

    Diagonal_Operator(int n_in, double v_in) : Vector_Operator(n_in, n_in), v(v_in) { };

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

    Diagonal_Operator_2(int n_in, int min_gid_in, double v_in)
    : Vector_Operator(n_in, n_in), min_gid(min_gid_in), v(v_in) {}

    ~Diagonal_Operator_2() { };

    void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y)
    {
      int myCols = y.MyLength();
      for (int j=0; j < x.NumVectors(); ++j) {
        for (int i=0; i < myCols; ++i) (*y(j))[i] = (min_gid+i+1)*v*(*x(j))[i];  // NOTE: square operator!
      }
    };

  private:

    int min_gid;
    double v;
};

//************************************************************************************************

class Composed_Operator : public Vector_Operator
{
  public:

    Composed_Operator(int n,
        const Teuchos::RCP<Vector_Operator>& pA_in,
        const Teuchos::RCP<Vector_Operator>& pB_in);

    virtual ~Composed_Operator() {};

    virtual void operator () (const Epetra_MultiVector &x, Epetra_MultiVector &y);

  private:

    Teuchos::RCP<Vector_Operator> pA;
    Teuchos::RCP<Vector_Operator> pB;
};

Composed_Operator::Composed_Operator(int n_in,
    const Teuchos::RCP<Vector_Operator>& pA_in,
    const Teuchos::RCP<Vector_Operator>& pB_in)
: Vector_Operator(n_in, n_in), pA(pA_in), pB(pB_in)
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

    Trilinos_Interface(const Teuchos::RCP<Vector_Operator>   pA_in,
        const Teuchos::RCP<const Epetra_Comm> pComm_in,
        const Teuchos::RCP<const Epetra_Map>  pMap_in)
      : pA (pA_in),
      pComm (pComm_in),
      pMap (pMap_in),
      use_transpose (false)
  {}

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      return(Apply(X,Y));  // No inverse
    };

    virtual ~Trilinos_Interface() {};

    const char * Label() const {return("Trilinos_Interface, an Epetra_Operator implementation");};

    bool UseTranspose() const {return(use_transpose);};      // always set to false

    int SetUseTranspose(bool UseTranspose_in) { use_transpose = false; return(-1); };

    bool HasNormInf() const {return(false);};                // cannot return inf-norm

    double NormInf() const {return(0.0);};

    virtual const Epetra_Comm & Comm() const {return *pComm; }

    virtual const Epetra_Map & OperatorDomainMap() const {return *pMap; }

    virtual const Epetra_Map & OperatorRangeMap() const {return *pMap; }

  private:

    Teuchos::RCP<Vector_Operator>   pA;
    Teuchos::RCP<const Epetra_Comm> pComm;
    Teuchos::RCP<const Epetra_Map>  pMap;

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

  Iterative_Inverse_Operator(int n_in, int blocksize,
      const Teuchos::RCP<Vector_Operator>& pA_in,
      std::string opString="Iterative Solver", bool print_in=false);

  virtual ~Iterative_Inverse_Operator() {}

  virtual void operator () (const Epetra_MultiVector &b, Epetra_MultiVector &x);

  private:

  Teuchos::RCP<Vector_Operator> pA;       // operator which will be inverted
  // supplies a matrix std::vector multiply
  const bool print;

  Teuchos::Time timer;
  Teuchos::RCP<Epetra_Comm> pComm;
  Teuchos::RCP<Epetra_Map>  pMap;

  Teuchos::RCP<OP> pPE;
  Teuchos::RCP<Teuchos::ParameterList>         pList;
  Teuchos::RCP<LinearProblem<double,MV,OP> >   pProb;
  Teuchos::RCP<BiCGStabSolMgr<double,MV,OP> >      pBelos;
};

Iterative_Inverse_Operator::Iterative_Inverse_Operator(int n_in, int blocksize,
    const Teuchos::RCP<Vector_Operator>& pA_in,
    std::string opString, bool print_in)
: Vector_Operator(n_in, n_in),      // square operator
  pA(pA_in),
  print(print_in),
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

  int max_iter = 100;
  double tol = 1.0e-10;
  int verbosity = Belos::Errors + Belos::Warnings;
  if (print)
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;

  pList = Teuchos::rcp( new Teuchos::ParameterList );
  pList->set( "Maximum Iterations", max_iter );
  pList->set( "Convergence Tolerance", tol );
  pList->set( "Verbosity", verbosity );

  pBelos = Teuchos::rcp( new BiCGStabSolMgr<double,MV,OP>(pProb, pList) );
}

void Iterative_Inverse_Operator::operator () (const Epetra_MultiVector &b, Epetra_MultiVector &x)
{
  int pid = pComm->MyPID();

  // Initialize the solution to zero
  x.PutScalar( 0.0 );

  // Reset the solver, problem, and status test for next solve (HKT)
  pProb->setProblem( Teuchos::rcp(&x, false), Teuchos::rcp(&b, false) );

  timer.start();
  Belos::ReturnType ret = pBelos->solve();
  timer.stop();

  if (pid == 0 && print) {
    if (ret == Belos::Converged)
    {
      std::cout << std::endl << "pid[" << pid << "] BiCGStab converged" << std::endl;
      std::cout << "Solution time: " << timer.totalElapsedTime() << std::endl;

    }
    else
      std::cout << std::endl << "pid[" << pid << "] BiCGStab did not converge" << std::endl;
  }
}

//************************************************************************************************
//************************************************************************************************

int main(int argc, char *argv[])
{
  int pid = -1;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  bool success = true;
  try {

    pid = Comm.MyPID();

    int n(10);
    int numRHS=1;

    Epetra_Map Map = Epetra_Map(n, 0, Comm);
    Epetra_MultiVector X(Map, numRHS, false), Y(Map, numRHS, false);
    X.PutScalar( 1.0 );

    // Inner computes inv(D2)*y
    Teuchos::RCP<Diagonal_Operator_2> D2 = Teuchos::rcp(new Diagonal_Operator_2(n, Map.MinMyGID(), 1.0));
    Iterative_Inverse_Operator A2(n, 1, D2, "Belos (inv(D2))", true);

    // should return x=(1, 1/2, 1/3, ..., 1/10)
    A2(X,Y);

    if (pid==0) {
      std::cout << "Vector Y should have all entries [1, 1/2, 1/3, ..., 1/10]" << std::endl;
    }
    Y.Print(std::cout);

    // Inner computes inv(D)*x
    Teuchos::RCP<Diagonal_Operator> D = Teuchos::rcp(new Diagonal_Operator(n, 4.0));
    Teuchos::RCP<Iterative_Inverse_Operator> Inner =
      Teuchos::rcp(new Iterative_Inverse_Operator(n, 1, D, "Belos (inv(D))", false));

    // Composed_Operator computed inv(D)*B*x
    Teuchos::RCP<Diagonal_Operator> B = Teuchos::rcp(new Diagonal_Operator(n, 4.0));
    Teuchos::RCP<Composed_Operator> C = Teuchos::rcp(new Composed_Operator(n, Inner, B));

    // Outer computes inv(C) = inv(inv(D)*B)*x = inv(B)*D*x = x
    Teuchos::RCP<Iterative_Inverse_Operator> Outer =
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
      std::cout << "Two-norm of std::vector (Y-1.0) : "<< norm_Y[0] << std::endl;

    success = (norm_Y[0] < 1e-10 && !Teuchos::ScalarTraits<double>::isnaninf( norm_Y[0] ) );

    if (success) {
      if (pid==0)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (pid==0)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
