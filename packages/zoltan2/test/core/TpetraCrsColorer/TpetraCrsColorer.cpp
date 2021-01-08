
#include "Tpetra_Core.hpp"
#include "Kokkos_Random.hpp"
#include "Zoltan2_TestHelpers.hpp"
#include "Zoltan2_TpetraCrsColorer.hpp"

// Class to test the Colorer utility
class ColorerTest {
public:
  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using matrix_t = Tpetra::CrsMatrix<zscalar_t>;
  using multivector_t = Tpetra::MultiVector<zscalar_t>;
  using execution_space_t = typename matrix_t::device_type::execution_space;

  ///////////////////////////////////////////////////////////
  // Construct the test:  
  //   Read or generate a matrix (JBlock) with default range and domain maps
  //   Construct identical matrix (JCyclic) with cyclic range and domain maps

  ColorerTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              int narg, char**arg)
  : symmetric(false)
  {
    int me = comm->getRank();
    int np = comm->getSize();

    // Process command line arguments
    bool distributeInput = true;
    std::string filename = "";
    size_t xdim = 10, ydim = 11, zdim = 12;

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
                   "Should Zoltan2 distribute the matrix as it is read?");
    cmdp.setOption("symmetric", "non-symmetric", &symmetric,
                   "Is the matrix symmetric?");
    cmdp.parse(narg, arg);

    // Get and store a matrix
    if (filename != "") {
      // Read from a file
      UserInputForTests uinput(".", filename, comm, true, distributeInput);
      JBlock = uinput.getUITpetraCrsMatrix();
    }
    else {
      // Generate a matrix
      UserInputForTests uinput(xdim, ydim, zdim, string("Laplace3D"), comm,
                               true, distributeInput);
      JBlock = uinput.getUITpetraCrsMatrix();
    }

    // Build same matrix with cyclic domain and range maps
    // To make range and domain maps differ for square matrices,
    // keep some processors empty in the cyclic maps

    size_t nIndices = std::max(JBlock->getGlobalNumCols(),
                               JBlock->getGlobalNumRows());
    Teuchos::Array<gno_t> indices(nIndices);

    Teuchos::RCP<const map_t> vMapCyclic = 
                 getCyclicMap(JBlock->getGlobalNumCols(), indices, np-1, comm);
    Teuchos::RCP<const map_t> wMapCyclic =
                 getCyclicMap(JBlock->getGlobalNumRows(), indices, np-2, comm);

    // Fill JBlock with random numbers for a better test.
    JBlock->resumeFill();
    auto local_matrix = JBlock->getLocalMatrix();
    auto local_graph = JBlock->getCrsGraph()->getLocalGraph();

    using IST = typename Kokkos::Details::ArithTraits<zscalar_t>::val_type;
    using pool_type = 
          Kokkos::Random_XorShift64_Pool<execution_space_t>;
    pool_type rand_pool(static_cast<uint64_t>(me));

    Kokkos::fill_random(local_matrix.values, rand_pool, 
                        static_cast<IST>(1.), static_cast<IST>(9999.));
    JBlock->fillComplete();

    // Make JCyclic:  same matrix with different Domain and Range maps
    auto lclMatrix = JBlock->getLocalMatrix();
    JCyclic = rcp(new matrix_t(JBlock->getLocalMatrix(),
                               JBlock->getRowMap(), JBlock->getColMap(),
                               vMapCyclic, wMapCyclic));
  }

  ////////////////////////////////////////////////////////////////
  bool run(const char* testname, Teuchos::ParameterList &params) {

    bool ok = true;

    params.set("symmetric", symmetric);

    // test with default maps
    ok = buildAndCheckSeedMatrix(testname, params, true);

    // test with cyclic maps
    ok &= buildAndCheckSeedMatrix(testname, params, false);

    return ok;
  }
    
  ///////////////////////////////////////////////////////////////
  bool buildAndCheckSeedMatrix(
    const char *testname,
    Teuchos::ParameterList &params,
    const bool useBlock
  )
  {
    int ierr = 0;

    // Pick matrix depending on useBlock flag
    Teuchos::RCP<matrix_t> J = (useBlock ? JBlock : JCyclic);
    int me = J->getRowMap()->getComm()->getRank();

    std::cout << "Running " << testname << " with "
              << (useBlock ? "Block maps" : "Cyclic maps")
              << std::endl;

    // Create a colorer
    Zoltan2::TpetraCrsColorer<matrix_t> colorer(J);

    colorer.computeColoring(params);

    // Check coloring
    if (!colorer.checkColoring()) {
      std::cout << testname << " with "
                << (useBlock ? "Block maps" : "Cyclic maps")
                << " FAILED: invalid coloring returned"
                << std::endl;
      return false;
    }

    // Compute seed matrix V -- the application wants this matrix
    // Dense matrix of 0/1 indicating the compression via coloring

    const int numColors = colorer.getNumColors();

    // Compute the seed matrix; applications want this seed matrix

    multivector_t V(J->getDomainMap(), numColors);
    colorer.computeSeedMatrix(V);

    // To test the result...
    // Compute the compressed matrix W
    multivector_t W(J->getRangeMap(), numColors);
  
    J->apply(V, W);

    // Reconstruct matrix from compression vector
    Teuchos::RCP<matrix_t> Jp = rcp(new matrix_t(*J, Teuchos::Copy));
    Jp->resumeFill();
    Jp->setAllToScalar(static_cast<zscalar_t>(-1.));
    Jp->fillComplete();

    colorer.reconstructMatrix(W, *Jp);

    // Check that values of J = values of Jp
    auto J_local_matrix = J->getLocalMatrix();
    auto Jp_local_matrix = Jp->getLocalMatrix();
    const size_t num_local_nz = J->getNodeNumEntries();

    Kokkos::parallel_reduce(
      "TpetraCrsColorer::testReconstructedMatrix()",
      Kokkos::RangePolicy<execution_space_t>(0, num_local_nz),
      KOKKOS_LAMBDA(const size_t nz, int &errorcnt) {
        if (J_local_matrix.values(nz) != Jp_local_matrix.values(nz)) {
          printf("Error in nonzero comparison %zu:  %g != %g", 
                  nz, J_local_matrix.values(nz), Jp_local_matrix.values(nz));
          errorcnt++;
        }
      }, 
      ierr);
   

    if (ierr > 0) {
      if (me == 0) {
        std::cout << testname << " FAILED with "
                  << (useBlock ? "Block maps" : "Cyclic maps")
                  << std::endl;
        params.print();
      }
    }

    return (ierr == 0);
  }
  
private:

  ////////////////////////////////////////////////////////////////
  // Return a map that is cyclic (like dealing rows to processors)
  Teuchos::RCP<const map_t> getCyclicMap(
    size_t nIndices, 
    Teuchos::Array<gno_t> &indices,
    int mapNumProc, 
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  {
    size_t cnt = 0;
    int me = comm->getRank();
    int np = comm->getSize();
    if (mapNumProc > np) mapNumProc = np; // corner case: bad input
    if (mapNumProc <= 0) mapNumProc = 1;  // corner case: np is too small

    for (size_t i = 0; i < nIndices; i++) 
      if (me == int(i % np)) indices[cnt++] = i;

    Tpetra::global_size_t dummy =
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

    return rcp(new map_t(dummy, indices(0,cnt), 0, comm));
  }

  ////////////////////////////////////////////////////////////////
  // Input matrix -- built in Constructor
  bool symmetric;               // User can specify whether matrix is symmetric
  Teuchos::RCP<matrix_t> JBlock;   // has Trilinos default domain and range maps
  Teuchos::RCP<matrix_t> JCyclic;  // has cyclic domain and range maps
};


///////////////////////////////////////////////////////////////////////////////
int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  bool ok = true;
  int ierr = 0;

  ColorerTest testColorer(comm, narg, arg);

  // Set parameters and compute coloring
  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";
    bool symmetrize = true;  // Zoltan's TRANSPOSE symmetrization, if needed

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetrize", symmetrize);

    testColorer.run("Test One", coloring_params);
    if (!ok) ierr++;
  }

  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";
    bool symmetrize = false;  // colorer's BIPARTITE symmetrization, if needed

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetrize", symmetrize);

    testColorer.run("Test Two", coloring_params);
    if (!ok) ierr++;
  }

  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";

    coloring_params.set("MatrixType", matrixType);

    testColorer.run("Test Three", coloring_params);
    if (!ok) ierr++;
  }

  if (ierr == 0)
    std::cout << "TEST PASSED" << std::endl;

//Through cmake...
//Test cases -- UserInputForTests can generate Galeri or read files:
//-  tri-diagonal matrix -- can check the number of colors
//-  galeri matrix
//-  read from file:  symmetric matrix and non-symmetric matrix

//Through code ...
//Test with fitted and non-fitted maps 
//Call regular and fitted versions of functions

//Through code ...
//Test both with and without Symmetrize -- 
//test both to exercise both sets of callbacks in Zoltan
// --matrixType = Jacobian/Hessian
// --symmetric, --no-symmetric
// --symmetrize, --no-symmetrize

//Through cmake
//Test both with and without distributeInput
// --distribute, --no-distribute

}
