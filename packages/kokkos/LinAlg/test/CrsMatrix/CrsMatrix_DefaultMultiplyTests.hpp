/** \file CrsMatrix_DefaultMultiplyTests.hpp
    \brief A file unit-testing and demonstrating the Kokkos default sparse kernel provider.
 */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( DefaultSparseOps, NodeTest, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps DSM;
    typedef CrsGraph<Ordinal,Node,DSM>                             GRPH;
    typedef CrsMatrix<Scalar,Ordinal,Node,DSM>                      MAT;
    typedef MultiVector<Scalar,Node>                                 MV;
    typedef Teuchos::ScalarTraits<Scalar>                            ST;
    typename DSM::template rebind<Scalar>::other dsm(node);
    TEST_EQUALITY(dsm.getNode(), node);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( DefaultSparseOps, ResubmitMatrix, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps DSM;
    typedef CrsGraph<Ordinal,Node,DSM>                             GRPH;
    typedef CrsMatrix<Scalar,Ordinal,Node,DSM>                      MAT;
    typedef MultiVector<Scalar,Node>                                 MV;
    typedef Teuchos::ScalarTraits<Scalar>                            ST;
    const size_t numRows = 5;
    GRPH G(numRows,node);
    MAT  A(G);
    A.finalize(true);
    out << "\n**\n** Can't submit the values before the structure\n**\n";
    {
      typename DSM::template rebind<Scalar>::other dsm(node);
      TEST_THROW( dsm.initializeValues(A), std::runtime_error );
    }
    out << "\n**\n** Can't submit the graph twice\n**\n";
    {
      typename DSM::template rebind<Scalar>::other dsm(node);
      dsm.initializeStructure(G);
      TEST_THROW( dsm.initializeStructure(G), std::runtime_error );
    }
    out << "\n**\n** Can submit the graph again after calling clear\n**\n";
    {
      typename DSM::template rebind<Scalar>::other dsm(node);
      dsm.initializeStructure(G);
      dsm.clear();
      TEST_NOTHROW( dsm.initializeStructure(G) );
    }
    out << "\n**\n** Can submit the values twice\n**\n";
    {
      typename DSM::template rebind<Scalar>::other dsm(node);
      dsm.initializeStructure(G);
      dsm.initializeValues(A);
      TEST_NOTHROW( dsm.initializeValues(A) );
    }
  }
