#include <Teuchos_TimeMonitor.hpp>

/** \file CrsMatrix_CrsTimingTests.hpp
    \brief A file unit-testing and demonstrating comparison of different CRS object types.
 */

  namespace Test {
    int numRows = 1000;
    int numIters = 1000;
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&Test::numRows,"Vector length for tests.");
    clp.setOption("test-iters",&Test::numIters,"Number of mat-vecs to time.");
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsTiming, PowerTriDiag, GRPH, MAT )
  {
    typedef typename MAT::ScalarType        Scalar;
    typedef typename MAT::OrdinalType       Ordinal;
    typedef typename MAT::NodeType          Node;
    typedef typename MAT::LocalMatOpsType   SparseOps;
    typedef MultiVector<Scalar,Node>                                 MV;
    typedef Teuchos::ScalarTraits<Scalar>                            ST;
    

    out << "Testing " << Teuchos::TypeNameTraits<GRPH>::name() << std::endl 
        << "    and " << Teuchos::TypeNameTraits<MAT>::name()  << std::endl;

    RCP<Node> node = getNode<Node>();

    GRPH G(Test::numRows,node);
    MAT  A(G);
    // 
    // allocate buffers and fill the graph and matrix
    // 
    {
      const size_t totalNNZ = 3*Test::numRows - 2;
      ArrayRCP<size_t> offsets(Test::numRows+1);
      ArrayRCP<Ordinal>   inds(totalNNZ);
      ArrayRCP<Scalar>    vals(totalNNZ);
      size_t NNZsofar = 0;
      offsets[0] = NNZsofar;
      inds[NNZsofar] = 0; inds[NNZsofar+1] =  1;
      vals[NNZsofar] = 2; vals[NNZsofar+1] = -1;
      NNZsofar += 2;
      for (int i=1; i != Test::numRows-1; ++i) {
        offsets[i] = NNZsofar;
        inds[NNZsofar] = i-1; inds[NNZsofar+1] = i; inds[NNZsofar+2] = i+1;
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 3; vals[NNZsofar+2] =  -1;
        NNZsofar += 3;
      }
      offsets[Test::numRows-1] = NNZsofar;
      inds[NNZsofar] = Test::numRows-2; inds[NNZsofar+1] = Test::numRows-1;
      vals[NNZsofar] =  -1;           vals[NNZsofar+1] = 2;
      NNZsofar += 2;
      offsets[Test::numRows]   = NNZsofar;
      TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
      G.set1DStructure(inds, offsets, offsets.persistingView(1,Test::numRows));
      offsets = Teuchos::null;
      inds    = Teuchos::null;
      A.set1DValues(vals);
      vals    = Teuchos::null;
      A.finalize(true);
    }
    // 
    // fill the matvec
    // 
    typename SparseOps::template rebind<Scalar>::other matvec(node);
    matvec.initializeStructure(G);
    matvec.initializeValues(A);
    // 
    // time the matvec
    // 
    MV X(node), Y(node);
    X.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
    Y.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
    DefaultArithmetic<MV>::Init( X, ST::one() );
    Teuchos::RCP<Teuchos::Time> matvectime = Teuchos::TimeMonitor::getNewTimer("LocalTimer");
    {
      Teuchos::TimeMonitor lcltimer(*matvectime);
      for (int i=0; i<Test::numIters; ++i) {
        if ( (i%2) == 0 ) {
          matvec.multiply(Teuchos::NO_TRANS,ST::one(),X,Y);
        }
        else {
          matvec.multiply(Teuchos::NO_TRANS,ST::one(),Y,X);
        }
      }
    }
    out << "Time is " << matvectime->totalElapsedTime() << std::endl;
  }

