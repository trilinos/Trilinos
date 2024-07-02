// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test generates diagonal matrices for TFQMR to solve.
//
// NOTE: No preconditioner is used in this case.
//

// Teuchos
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTpetraOperator.hpp"
#include "BelosTFQMRSolMgr.hpp"

using namespace Belos;

//************************************************************************************************

template<class MV>
class VectorOperator
{
  public:

    VectorOperator(int m_in, int n_in) : m(m_in), n(n_in) {};

    virtual ~VectorOperator() {};

    virtual void operator () (const MV &x, MV &y) = 0;

    int size (int dim) const { return (dim == 1) ? m : n; };

  protected:

    int m, n;        // an (m x n) operator

  private:

    // Not allowing copy construction.
    VectorOperator( const VectorOperator& ): m(0), n(0) {};
    VectorOperator* operator=( const VectorOperator& ) { return NULL; };

};

//************************************************************************************************

template<class ST, class MV>
class DiagonalOperator : public VectorOperator<MV>
{
  public:

    DiagonalOperator(int n_in, ST v_in) : VectorOperator<MV>(n_in, n_in), v(v_in) { };

    ~DiagonalOperator() { };

    void operator () (const MV &x, MV &y)
    {
      y.scale(v, x);
    };

  private:

    ST v;
};

//************************************************************************************************

template<class ST, class MV>
class DiagonalOperator2 : public VectorOperator<MV>
{
  public:

    DiagonalOperator2<ST,MV>(int n_in, int min_gid_in, ST v_in)
    : VectorOperator<MV>(n_in, n_in), min_gid(min_gid_in), v(v_in) {}

    ~DiagonalOperator2() { };

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
class ComposedOperator : public VectorOperator<MV>
{
  public:

    ComposedOperator(int n,
        const Teuchos::RCP<VectorOperator<MV>>& pA_in,
        const Teuchos::RCP<VectorOperator<MV>>& pB_in);

    virtual ~ComposedOperator() {};

    virtual void operator () (const MV &x, MV &y);

  private:

    Teuchos::RCP<VectorOperator<MV>> pA;
    Teuchos::RCP<VectorOperator<MV>> pB;
};

template<class MV>
ComposedOperator<MV>::ComposedOperator(int n_in,
    const Teuchos::RCP<VectorOperator<MV>>& pA_in,
    const Teuchos::RCP<VectorOperator<MV>>& pB_in)
: VectorOperator<MV>(n_in, n_in), pA(pA_in), pB(pB_in)
{
}

template<class MV>
void ComposedOperator<MV>::operator () (const MV &x, MV &y)
{
  MV ytemp(y.getMap(), y.getNumVectors(), false);
  (*pB)( x, ytemp );
  (*pA)( ytemp, y );
}

//************************************************************************************************

template<class OP, class ST, class MP, class MV>
class TrilinosInterface : public OP
{
  public:

    TrilinosInterface(const Teuchos::RCP<VectorOperator<MV>> pA_in,
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

    virtual ~TrilinosInterface() {};

    bool hasTransposeApply() const {return(use_transpose);};

    Teuchos::RCP<const MP> getDomainMap() const override {return pMap; }

    Teuchos::RCP<const MP> getRangeMap() const override {return pMap; }

  private:

    Teuchos::RCP<VectorOperator<MV>> pA;
    Teuchos::RCP<const Teuchos::Comm<int>> pComm;
    Teuchos::RCP<const MP> pMap;

    bool use_transpose;
};

template<class OP, class ST, class MP, class MV>
void TrilinosInterface<OP, ST, MP, MV>::apply (const MV &X, MV &Y,
                                                Teuchos::ETransp mode, ST alpha, ST beta) const
{
  (*pA)(X,Y);
}

//************************************************************************************************

template<class OP, class ST, class MP, class MV>
class IterativeInverseOperator : public VectorOperator<MV>
{
  public:

  IterativeInverseOperator(int n_in, int blocksize,
    const Teuchos::RCP<VectorOperator<MV>>& pA_in,
    std::string opString="Iterative Solver", bool print_in=false);

  virtual ~IterativeInverseOperator() {}

  virtual void operator () (const MV &b, MV &x);

  private:

  Teuchos::RCP<VectorOperator<MV>> pA;       // operator which will be inverted
  // supplies a matrix std::vector multiply
  const bool print;

  Teuchos::Time timer;
  Teuchos::RCP<const Teuchos::Comm<int>> pComm;
  Teuchos::RCP<MP> pMap;

  Teuchos::RCP<OP> pPE;
  Teuchos::RCP<Teuchos::ParameterList>   pList;
  Teuchos::RCP<LinearProblem<ST,MV,OP> > pProb;
  Teuchos::RCP<TFQMRSolMgr<ST,MV,OP> >   pBelos;
};

template<class OP, class ST, class MP, class MV>
IterativeInverseOperator<OP, ST, MP, MV>::IterativeInverseOperator(int n_in, int blocksize,
  const Teuchos::RCP<VectorOperator<MV>>& pA_in,
  std::string opString, bool print_in)
: VectorOperator<MV>(n_in, n_in),      // square operator
  pA(pA_in),
  print(print_in),
  timer(opString)
{
  int nGlobal;

  pComm = Tpetra::getDefaultComm();
  Teuchos::reduceAll<int, int>(*pComm, Teuchos::REDUCE_SUM, 1, &n_in, &nGlobal);

  pMap = Teuchos::rcp( new MP(nGlobal, n_in, 0, pComm) );
  pPE = Teuchos::rcp( new TrilinosInterface<OP, ST, MP, MV>(pA, pComm, pMap ) );

  pProb = Teuchos::rcp( new LinearProblem<ST,MV,OP>() );
  pProb->setOperator( pPE );

  int max_iter = 100;
  ST tol = sqrt(std::numeric_limits<ST>::epsilon());
  int verbosity = Belos::Errors + Belos::Warnings;
  if (print)
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;

  pList = Teuchos::rcp( new Teuchos::ParameterList );
  pList->set( "Maximum Iterations", max_iter );
  pList->set( "Convergence Tolerance", tol );
  pList->set( "Verbosity", verbosity );

  pBelos = Teuchos::rcp( new TFQMRSolMgr<ST,MV,OP>(pProb, pList) );
}

template<class OP, class ST, class MP, class MV>
void IterativeInverseOperator<OP, ST, MP, MV>::operator () (const MV &b, MV &x)
{
  int pid = pComm->getRank();

  // Initialize the solution to zero
  x.putScalar( 0.0 );

  // Reset the solver, problem, and status test for next solve (HKT)
  pProb->setProblem( Teuchos::rcp(&x, false), Teuchos::rcp(&b, false) );

  timer.start();
  Belos::ReturnType ret = pBelos->solve();
  timer.stop();

  if (pid == 0 && print) {
    if (ret == Belos::Converged)
    {
      std::cout << std::endl << "pid[" << pid << "] TFQMR converged" << std::endl;
      std::cout << "Solution time: " << timer.totalElapsedTime() << std::endl;

    }
    else
      std::cout << std::endl << "pid[" << pid << "] TFQMR did not converge" << std::endl;
  }
}

//************************************************************************************************
//************************************************************************************************

template<class ScalarType>
int run(int argc, char *argv[])
{
  // Get default Tpetra template types
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  // Init Tpetra types
  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MP = typename Tpetra::Map<LO,GO,NT>;

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int pid = rank(*comm);

  bool verbose = false;
  bool success = true;

  ST tol = sqrt(std::numeric_limits<ST>::epsilon());

  try {
    int n(10);
    int numRHS=1;

    RCP<MP> map = RCP(new MP(n, 0, comm));
    MV X(map, numRHS), Y(map, numRHS);

    X.putScalar(1.0);

    // Inner computes inv(D2)*y
    RCP<DiagonalOperator2<ST,MV>> D2 = rcp(new DiagonalOperator2<ST,MV>(n, map->getMinGlobalIndex(), 1.0));
    IterativeInverseOperator<OP, ST, MP, MV> A2(n, 1, D2, "Belos (inv(D2))", true);

    // should return x=(1, 1/2, 1/3, ..., 1/10)
    A2(X,Y);

    if (pid==0) {
      std::cout << "Vector Y should have all entries [1, 1/2, 1/3, ..., 1/10]" << std::endl;
    }
    Y.print(std::cout);

    // Inner computes inv(D)*x
    RCP<DiagonalOperator<ST,MV>> D = rcp(new DiagonalOperator<ST,MV>(n, 4.0));
    RCP<IterativeInverseOperator<OP, ST, MP, MV>> Inner =
      rcp(new IterativeInverseOperator<OP, ST, MP, MV>(n, 1, D, "Belos (inv(D))", false));

    // Composed_Operator computed inv(D)*B*x
    RCP<DiagonalOperator<ST,MV>> B = rcp(new DiagonalOperator<ST,MV>(n, 4.0));
    RCP<ComposedOperator<MV>> C = rcp(new ComposedOperator<MV>(n, Inner, B));

    // Outer computes inv(C) = inv(inv(D)*B)*x = inv(B)*D*x = x
    RCP<IterativeInverseOperator<OP, ST, MP, MV>> Outer =
      rcp(new IterativeInverseOperator<OP, ST, MP, MV>(n, 1, C, "Belos (inv(C)=inv(inv(D)*B))", true));

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
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // run<float>(argc, argv);
}
