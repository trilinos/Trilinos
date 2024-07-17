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

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosTpetraOperator.hpp>
#include <BelosBiCGStabSolMgr.hpp>

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_Vector_fwd.hpp>
#include <Tpetra_CrsMatrix_fwd.hpp>
#include <Tpetra_MultiVector_fwd.hpp>

// Teuchos
#include <Teuchos_Time.hpp>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


using std::vector;
using Teuchos::RCP;
using Teuchos::rcp;

using namespace Belos;

//************************************************************************************************

template<class MV>
class Vector_Operator
{
  public:

    Vector_Operator(int m_in, int n_in) : m(m_in), n(n_in) {};

    virtual ~Vector_Operator() {};

    virtual void operator () (const MV &x, MV &y) = 0;

    int size (int dim) const { return (dim == 1) ? m : n; };

  protected:

    int m, n;        // an (m x n) operator

  private:

    // Not allowing copy construction.
    Vector_Operator( const Vector_Operator& ): m(0), n(0) {};
    Vector_Operator* operator=( const Vector_Operator& ) { return nullptr; };

};

//************************************************************************************************

template<class ST, class MV>
class Diagonal_Operator : public Vector_Operator<MV>
{
  public:

    Diagonal_Operator(int n_in, ST v_in) : Vector_Operator<MV>(n_in, n_in), v(v_in) { };

    ~Diagonal_Operator() { };

    void operator () (const MV &x, MV &y)
    {
      y.scale( v, x );
    };

  private:

    ST v;
};

//************************************************************************************************

template<class ST, class MV>
class Diagonal_Operator_2 : public Vector_Operator<MV>
{
  public:

    Diagonal_Operator_2<ST,MV>(int n_in, int min_gid_in, ST v_in)
    : Vector_Operator<MV>(n_in, n_in), min_gid(min_gid_in), v(v_in) {}

    ~Diagonal_Operator_2() { };

    void operator () (const MV &x, MV &y)
    {
      auto yLocalData = y.getLocalViewHost(Tpetra::Access::ReadWrite);
      auto xLocalData = x.getLocalViewHost(Tpetra::Access::ReadOnly);

      for (size_t j = 0; j < x.getNumVectors(); ++j) {
        for (size_t i = 0; i < x.getLocalLength(); ++i) {
            yLocalData(i, j) = (min_gid + i + 1) * v * xLocalData(i, j); // NOTE: square operator!
        }
      }
    }

  private:

    int min_gid;
    ST v;
};

//************************************************************************************************

template<class MV>
class Composed_Operator : public Vector_Operator<MV>
{
  public:

    Composed_Operator(int n,
        const RCP<Vector_Operator<MV>>& pA_in,
        const RCP<Vector_Operator<MV>>& pB_in);

    virtual ~Composed_Operator() {};

    virtual void operator () (const MV &x, MV &y);

  private:

    RCP<Vector_Operator<MV>> pA;
    RCP<Vector_Operator<MV>> pB;
};

template<class MV>
Composed_Operator<MV>::Composed_Operator(int n_in,
    const RCP<Vector_Operator<MV>>& pA_in,
    const RCP<Vector_Operator<MV>>& pB_in)
: Vector_Operator<MV>(n_in, n_in), pA(pA_in), pB(pB_in)
{
}

template<class MV>
void Composed_Operator<MV>::operator () (const MV &x, MV &y)
{
  MV ytemp(y.getMap(), y.getNumVectors(), false);
  (*pB)( x, ytemp );
  (*pA)( ytemp, y );
}

//************************************************************************************************

template<class OP, class ST, class MP, class MV>
class Trilinos_Interface : public OP
{
  public:

    Trilinos_Interface(const Teuchos::RCP<Vector_Operator<MV>> pA_in,
        const Teuchos::RCP<const Teuchos::Comm<int>> pComm_in,
        const Teuchos::RCP<const MP> pMap_in)
      : pA (pA_in),
      pComm (pComm_in),
      pMap (pMap_in),
      use_transpose (false)
  {}

    void apply (const MV &X,
                MV &Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                ST alpha = Teuchos::ScalarTraits<ST>::one(),
                ST beta = Teuchos::ScalarTraits<ST>::zero()) const override;

    virtual ~Trilinos_Interface() {};

    bool hasTransposeApply() const {return(use_transpose);};      // always set to false (in fact the default)

    Teuchos::RCP<const MP> getDomainMap() const override {return pMap; }
    Teuchos::RCP<const MP> getRangeMap() const override {return pMap; }

  private:

    Teuchos::RCP<Vector_Operator<MV>> pA;
    Teuchos::RCP<const Teuchos::Comm<int>> pComm;
    Teuchos::RCP<const MP> pMap;

    bool use_transpose;
};

template<class OP, class ST, class MP, class MV>
void Trilinos_Interface<OP, ST, MP, MV>::apply (const MV &X,
            MV &Y,
            Teuchos::ETransp mode,
            ST alpha,
            ST beta) const {
    (*pA)(X,Y);
}

//************************************************************************************************

template<class OP, class ST, class MP, class MV>
class Iterative_Inverse_Operator : public Vector_Operator<MV>
{
  public:

  Iterative_Inverse_Operator(int n_in, int blocksize,
      const RCP<Vector_Operator<MV>>& pA_in,
      std::string opString="Iterative Solver", bool print_in=false);

  virtual ~Iterative_Inverse_Operator() {}

  virtual void operator () (const MV &b, MV &x);

  private:

  RCP<Vector_Operator<MV>> pA;       // operator which will be inverted
  // supplies a matrix std::vector multiply
  const bool print;

  Teuchos::Time timer;
  RCP<const Teuchos::Comm<int>> pComm;
  RCP<MP> pMap;

  RCP<OP> pPE;
  RCP<Teuchos::ParameterList> pList;
  RCP<LinearProblem<ST,MV,OP> > pProb;
  RCP<BiCGStabSolMgr<ST,MV,OP> > pBelos;
};

template<class OP, class ST, class MP, class MV>
Iterative_Inverse_Operator<OP,ST,MP,MV>::Iterative_Inverse_Operator(int n_in, int blocksize,
    const RCP<Vector_Operator<MV>>& pA_in,
    std::string opString, bool print_in)
: Vector_Operator<MV>(n_in, n_in),      // square operator
  pA(pA_in),
  print(print_in),
  timer(opString)
{
  int n_global;
  const int count = 1;
  const auto reductType = Teuchos::REDUCE_SUM;

  pComm = Tpetra::getDefaultComm();
  Teuchos::reduceAll<int,int>( *pComm, reductType, count, &n_in, &n_global );

  pMap =  rcp( new MP(n_global, n_in, 0, pComm) );

  pPE = rcp( new Trilinos_Interface<OP,ST,MP,MV>(pA, pComm, pMap) );

  pProb = rcp( new LinearProblem<ST,MV,OP>() );
  pProb->setOperator( pPE );

  int max_iter = 100;
  ST tol = 1.0e-10;
  int verbosity = Belos::Errors + Belos::Warnings;
  if (print)
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;

  pList = rcp( new Teuchos::ParameterList );
  pList->set( "Maximum Iterations", max_iter );
  pList->set( "Convergence Tolerance", tol );
  pList->set( "Verbosity", verbosity );

  pBelos = rcp( new BiCGStabSolMgr<ST,MV,OP>(pProb, pList) );
}

template<class OP, class ST, class MP, class MV>
void Iterative_Inverse_Operator<OP,ST,MP,MV>::operator () (const MV &b, MV &x)
{
  int pid = pComm->getRank();

  // Initialize the solution to zero
  x.putScalar( 0.0 );

  // Reset the solver, problem, and status test for next solve (HKT)
  pProb->setProblem( rcp(&x, false), rcp(&b, false) );

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

template<class ScalarType>
int run(int argc, char *argv[])
{
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  using MP = typename Tpetra::Map<LO,GO,NT>;
  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;

  int pid = -1;

  Teuchos::GlobalMPISession session(&argc, &argv, nullptr);

  const auto comm = Tpetra::getDefaultComm();

  bool verbose = false;
  bool success = true;

  ST tol = 1.0e-10;

  try {

    pid = comm->getRank();

    int n(10);
    int numRHS=1;

    RCP<MP> Map = rcp(new MP(n, 0, comm));
    MV X(Map, numRHS, false), Y(Map, numRHS, false);
    X.putScalar( 1.0 );

    // Inner computes inv(D2)*y
    RCP<Diagonal_Operator_2<ST,MV>> D2 = rcp(new Diagonal_Operator_2<ST,MV>(n, Map->getMinGlobalIndex(), 1.0));
    Iterative_Inverse_Operator<OP,ST,MP,MV> A2(n, 1, D2, "Belos (inv(D2))", true);

    // should return x=(1, 1/2, 1/3, ..., 1/10)
    A2(X,Y);

    if (pid==0) {
      std::cout << "Vector Y should have all entries [1, 1/2, 1/3, ..., 1/10]" << std::endl;
    }
    Y.print(std::cout);

    // Inner computes inv(D)*x
    RCP<Diagonal_Operator<ST,MV>> D = rcp(new Diagonal_Operator<ST,MV>(n, 4.0));
    RCP<Iterative_Inverse_Operator<OP,ST,MP,MV>> Inner =
      rcp(new Iterative_Inverse_Operator<OP,ST,MP,MV>(n, 1, D, "Belos (inv(D))", false));

    // Composed_Operator computed inv(D)*B*x
    RCP<Diagonal_Operator<ST,MV>> B = rcp(new Diagonal_Operator<ST,MV>(n, 4.0));
    RCP<Composed_Operator<MV>> C = rcp(new Composed_Operator<MV>(n, Inner, B));

    // Outer computes inv(C) = inv(inv(D)*B)*x = inv(B)*D*x = x
    RCP<Iterative_Inverse_Operator<OP,ST,MP,MV>> Outer =
      rcp(new Iterative_Inverse_Operator<OP,ST,MP,MV>(n, 1, C, "Belos (inv(C)=inv(inv(D)*B))", true));

    // should return x=1/4
    (*Inner)(X,Y);

    if (pid==0) {
      std::cout << std::endl << "Vector Y should have all entries [1/4, 1/4, 1/4, ..., 1/4]" << std::endl;
    }
    Y.print(std::cout);

    // should return x=1
    (*Outer)(X,Y);

    if (pid==0) {
      std::cout << "Vector Y should have all entries [1, 1, 1, ..., 1]" << std::endl;
    }
    Y.print(std::cout);

    // Compute the norm of Y - 1.0
    std::vector<ST> norm_Y(Y.getNumVectors());
    Y.update(-1.0, X, 1.0);
    Y.norm2(norm_Y);

    if (pid==0)
      std::cout << "Two-norm of std::vector (Y-1.0) : "<< norm_Y[0] << std::endl;

    success = (norm_Y[0] < tol && !Teuchos::ScalarTraits<ST>::isnaninf( norm_Y[0] ) );

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

int main(int argc, char *argv[]) {
  // run with different ST (which may require different tolerances)
  return run<double>(argc, argv);
  // return run<float>(argc, argv); // FAILS
}

