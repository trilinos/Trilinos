////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonlocalAfterResume, LO, GO, Scalar, Node )
{
  RCP<Node> node = getNode<Node>();
  // test that an exception is thrown when we exceed statically allocated memory
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t numImages = size(*comm);
  const size_t myImageID = rank(*comm);
  // create a row Map, 5 rows per processor
  const size_t numLocal = 5;
  RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm,node);
  RCP<const Map<LO,GO,Node> > cmap;
  // create a column Map, with super- and sub-diagonal blocks
  {
    Array<GO> cols;
    for (GO c=rmap->getMinGlobalIndex(); c <= rmap->getMaxGlobalIndex(); ++c) cols.push_back(c);
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex() + numLocal) {
      for (GO c = rmap->getMinGlobalIndex()-numLocal; c < rmap->getMinGlobalIndex(); ++c) {
        cols.push_back(c);
      }
    }
    if (rmap->getMaxGlobalIndex()+numLocal <= rmap->getMaxAllGlobalIndex()) {
      for (GO c = rmap->getMaxGlobalIndex()+1; c <= rmap->getMaxGlobalIndex()+numLocal; ++c) {
        cols.push_back(c);
      }
    }
    cmap = createNonContigMapWithNode<LO,GO,Node>(cols(), comm, node);
  }
  {
    //----------------------------------------------------------------------
    // put in diagonal, locally
    //----------------------------------------------------------------------
    Tpetra::CrsMatrix<Scalar,LO,GO,Node> matrix(rmap,cmap,3,DynamicProfile);
    for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
      matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
    }
    // fill, but do not pack, because we will add new entries below
    TEST_NOTHROW       ( matrix.fillComplete( DoNotOptimizeStorage ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), numLocal*numImages );
    TEST_EQUALITY      ( matrix.getNodeNumEntries(),   numLocal           );

    //----------------------------------------------------------------------
    // add super-diagonal, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex()+numLocal < rmap->getMaxAllGlobalIndex()) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r+numLocal,tuple(r),tuple(ST::one()));
      }
    }
    // fill, but do not pack, because we will add new entries below
    TEST_NOTHROW       ( matrix.fillComplete( DoNotOptimizeStorage ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 2*numLocal*numImages-numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0) expected += numLocal; // super-diagonal
      TEST_EQUALITY( matrix.getNodeNumEntries(), expected );
    }

    //----------------------------------------------------------------------
    // add sub-diagonal block, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex()+numLocal) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r-numLocal,tuple(r),tuple(ST::one()));
      }
    }
    // fill; it is okay to pack now
    TEST_NOTHROW       ( matrix.fillComplete( DoOptimizeStorage ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 3*numLocal*numImages-2*numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0)           expected += numLocal; // super-diagonal
      if (myImageID < numImages-1) expected += numLocal; // sub-diagonal
      TEST_EQUALITY( matrix.getNodeNumEntries(), expected );
    }
  }
  // All procs fail if any node fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );
}
