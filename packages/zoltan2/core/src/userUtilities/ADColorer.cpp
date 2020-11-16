Zoltan2::TpetraCrsColorer<matrix_t> colorer(J);
{

  if (Jacobian)
    if IanNotReady && !testZoltan2anyway
      ZoltanCrsColorer(PartialDistance2)
    else if IanReady || testZoltan2anyway
      Zoltan2CrsColorer(PartialDistance2)

  else if (Hessian)
    // Hessian is already symmetric; 
    // code currently still uses partial distance 2.
    // Should use Distance2 instead?
    if IanNotReady && !testZoltan2anyway
      ZoltanCrsColorer(Distance2)
    else if IanReady || testZoltan2anyway
      Zoltan2CrsColorer(Distance2)


}

class ColorerTest {
public:
  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using matrix_t = Tpetra::CrsMatrix<zscalar_t>;
  using multivector_t = Tpetra::MultiVector<zscalar_t>;

  ///////////////////////////////////////////////////////////
  ColorerTest(Teuchos::RCP<const Teucho::Comm<int> > &comm) {
    // Process command line arguments
    bool distributeInput = true;
    std::string filename = "";
    size_t xdim = 1, ydim = 1, zdim = 1;

    Teuchos::CommandLineProcessor cmdp(false, false);
    cmdp.setOption("file", &filename, 
                   "Name of the Matrix Market file to use");
    cmdp.setOption("xdim", &xdim, 
                   "Number of nodes in x-direction for generated matrix");
    cmdp.setOption("ydim", &ydim, 
                   "Number of nodes in y-direction for generated matrix");
    cmdp.setOption("zdim", &zdim, 
                   "Number of nodes in z-direction for generated matrix");
    cmdp.setOption("distribute", "no-distribute", &distributeInput, 
                   "Should Zoltan2 symmetrize the matrix?"
    cmdp.parse(narg, arg);

    // Get and store a matrix
    if (filename != "") {
      // Read from a file
      UserInputForTests uinput(".", filename, comm, true, distributeInput);
      J = uinput.getUITpetraCrsMatrix();
    }
    else {
      // Generate a matrix
      UserInputForTests uinput(xdim, ydim, zdim, string("Laplace3D"), comm,
                               true, distributeInput);
      J = uinput.getUITpetraCrsMatrix();
    }

    vMapFitted = J->getDomainMap();
    wMapFitted = J->getRangeMap();

    size_t nIndices = std::max(J->getGlobalNumCols(), J->getGlobalNumRows()));
    Teuchos::Array<gno_t> indices(nIndices);

    vMapNotFitted = getNotFittedMap(J->getGlobalNumCols(), indices);
    wMapNotFitted = getNotFittedMap(J->getGlobalNumRows(), indices);
  }

  ////////////////////////////////////////////////////////////////
  bool run(const char* testname, Teuchos::ParameterList *params) {

    bool ok = true;

    // Create a colorer
    Zoltan2::TpetraCrsADColorer<matrix_t> colorer(J);

    colorer.computeColoring(params);

    // Check coloring
    if (! colorer.checkColoring()) {
      if (me == 0)
        std::cout << testname << " FAILED: invalid coloring returned"
                  << std::endl;
      return false;
    }

    // Compute seed matrix V -- the application wants this matrix
    // Dense matrix of 0/1 indicating the compression via coloring

    const int numColors = colorer.getNumColors();

    // test with fitted maps
    ok = buildAndCheckSeedMatrix(testname, colorer, params, true);
    printCheckStatus(testname, ok, params, true);

    // test with non-fitted maps
    ok = buildAndCheckSeedMatrix(testname, colorer, params, false);
    printCheckStatus(testname, ok, params);
  }
    
  
private:

  ////////////////////////////////////////////////////////////////
  // Return a map that is cyclic (like dealing rows to processors)
  Teuchos::RCP<const map_t> getNotFittedMap(
    size_t nIndices, 
    Teuchos::Array<gno_t> &indices)
  {
    size_t cnt = 0;
    for (size_t i = 0; i < nIndices; i++) 
      if (me == int(i % np)) indices[cnt++] = i;

    Tpetra::global_size_t dummy =
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

    return rcp(new map_t(dummy, indices(0,cnt), 0, comm);
  }

  ///////////////////////////////////////////////////////////////
  bool buildAndCheckSeedMatrix(
    const char *testname,
    Zoltan2::TpetraCrsADColorer<matrix_t> &colorer,
    bool useFitted
  )
  {
    bool ok = true;

    // Pick maps for V, W depending on useFitted flag
    Teuchos::RCP<const map_t> vmap, wmap;
    if (useFitted) {
      vmap = vMapFitted;
      wmap = wMapFitted;
    }
    else {
      vmap = vMapNotFitted;
      wmap = wMapNotFitted;
    } 

    // Compute the seed matrix; applications want this matrix

    int numColors = colorer.getNumColors();
    multivector_t V(vmap, numColors);

    colorer.computeSeedMatrix(V);

    // To test the result...
    // Compute compression vector
    multivector_t W(wmap, numColors);  // W is the compressed matrix
  
    J->apply(W, V);

    // Reconstruct matrix from compression vector
    // TODO matrix_t Jp();
    // TODO colorer.reconstructMatrix(W, Jp);
  
    // Check J = Jp somehow
    // KDD is there a way to do this comparison in Tpetra?

    if (!ok) {
      if (me == 0) {
        std::cout << testname << " FAILED:  "
                  << " fitted=" << useFitted << std::endl;
        params.print();
      }
    }

    return ok;
  }

  ////////////////////////////////////////////////////////////////
  // Input matrix -- built in Constructor
  Teuchos::RCP<const matrix_t> J;

  // Test with V and W having maps that are LocallyFitted
  Teuchos::RCP<const map_t> vMapFitted;
  Teuchos::RCP<const map_t> wMapFitted;

  // Test with V and W having NOT maps that are LocallyFitted
  Teuchos::RCP<const map_t> vMapNotFitted;
  Teuchos::RCP<const map_t> wMapNotFitted;
}

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(narg, arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  bool ok = true;
  int ierr = 0;

  ColorTest<zscalar_t> testColorer(J);

  // Set parameters and compute coloring
  {
    Teuchos::ParameterList coloring_params;
    bool symmetric = true;
    bool symmetrize = true;
    std::string matrixType = "Jacobian";

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetric", symmetric);
    coloring_params.set("symmetrize", symmetrize);

    testColorer.run("Test One", coloring_params);
    if (!ok) ierr++;
  }



}

Test cases -- UserInputForTests can generate Galeri or read files:
-  tri-diagonal matrix -- can check the number of colors
-  galeri matrix
-  read from file:  symmetric matrix and non-symmetric matrix

Test with fitted and non-fitted maps
Call regular and fitted versions of functions

Test both with and without Symmetrize -- 
test both to exercise both sets of callbacks in Zoltan
 --symmetric, --no-symmetric
 --symmetrize, --no-symmetrize

Test both with and without distributeInput
 --distribute, --no-distribute

