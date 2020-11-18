#include "Tpetra_Core.hpp"
#include "Zoltan2_TestHelpers.hpp"
#include "Zoltan2_TpetraCrsADColorer.hpp"

// Class to test the Colorer utility
class ColorerTest {
public:
  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using matrix_t = Tpetra::CrsMatrix<zscalar_t>;
  using multivector_t = Tpetra::MultiVector<zscalar_t>;

  ///////////////////////////////////////////////////////////
  // Construct the test:  
  //   Read or generate a matrix with fitted range and domain maps
  //   Construct identical matrix with non-fitted range and domain maps

  ColorerTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              int narg, char**arg)
  {
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
                   "Should Zoltan2 symmetrize the matrix?");
    cmdp.parse(narg, arg);

    // Get and store a matrix
    if (filename != "") {
      // Read from a file
      UserInputForTests uinput(".", filename, comm, true, distributeInput);
      JFitted = uinput.getUITpetraCrsMatrix();
    }
    else {
      // Generate a matrix
      UserInputForTests uinput(xdim, ydim, zdim, string("Laplace3D"), comm,
                               true, distributeInput);
      JFitted = uinput.getUITpetraCrsMatrix();
    }

    // Build same matrix with non-fitted domain and range maps

    size_t nIndices = std::max(JFitted->getGlobalNumCols(),
                               JFitted->getGlobalNumRows());
    Teuchos::Array<gno_t> indices(nIndices);

    Teuchos::RCP<const map_t> vMapNotFitted = 
                 getNotFittedMap(JFitted->getGlobalNumCols(), indices, comm);
    Teuchos::RCP<const map_t> wMapNotFitted =
                 getNotFittedMap(JFitted->getGlobalNumRows(), indices, comm);

    JNotFitted = rcp(new matrix_t(*JFitted));
    JNotFitted->resumeFill();
    JNotFitted->fillComplete(vMapNotFitted, wMapNotFitted);
  }

  ////////////////////////////////////////////////////////////////
  bool run(const char* testname, Teuchos::ParameterList &params) {

    bool ok = true;

    // test with fitted maps
    ok = buildAndCheckSeedMatrix(testname, params, true);

    // test with non-fitted maps
    ok &= buildAndCheckSeedMatrix(testname, params, false);

    return ok;
  }
    
  
private:

  ////////////////////////////////////////////////////////////////
  // Return a map that is cyclic (like dealing rows to processors)
  Teuchos::RCP<const map_t> getNotFittedMap(
    size_t nIndices, 
    Teuchos::Array<gno_t> &indices,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  {
    size_t cnt = 0;
    int me = comm->getRank();
    int np = comm->getSize();

    for (size_t i = 0; i < nIndices; i++) 
      if (me == int(i % np)) indices[cnt++] = i;

    Tpetra::global_size_t dummy =
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

    return rcp(new map_t(dummy, indices(0,cnt), 0, comm));
  }

  ///////////////////////////////////////////////////////////////
  bool buildAndCheckSeedMatrix(
    const char *testname,
    Teuchos::ParameterList &params,
    const bool useFitted
  )
  {
    bool ok = true;

    // Pick matrix depending on useFitted flag
    Teuchos::RCP<matrix_t> J = (useFitted ? JFitted : JNotFitted);
    int me = J->getRowMap()->getComm()->getRank();

    // Create a colorer
    Zoltan2::TpetraCrsADColorer<matrix_t> colorer(J);

    colorer.computeColoring(params);

    // Check coloring
    if (!colorer.checkColoring()) {
      if (me == 0)
        std::cout << testname << " FAILED: invalid coloring returned"
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
    // Compute compression vector

    multivector_t W(J->getRangeMap(), numColors);  // W is the compressed matrix
  
    J->apply(W, V);

    // Reconstruct matrix from compression vector
    // TODO matrix_t Jp();
    // TODO colorer.reconstructMatrix(W, Jp);
    // TODO KDD:  I do not understand the operation here
  
    // Check J = Jp somehow
    // KDD is there a way to do this comparison in Tpetra?

    if (!ok) {
      if (me == 0) {
        std::cout << testname << " FAILED "
                  << (useFitted ? "with fitted maps" : "with non-fitted maps")
                  << std::endl;
        params.print();
      }
    }

    return ok;
  }

  ////////////////////////////////////////////////////////////////
  // Input matrix -- built in Constructor
  Teuchos::RCP<matrix_t> JFitted;     // has locally fitted domain and 
                                      // range maps
  Teuchos::RCP<matrix_t> JNotFitted;  // does NOT have locally fitted
                                      // domain and range maps
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
    bool symmetric = true;
    bool symmetrize = true;
    std::string matrixType = "Jacobian";

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetric", symmetric);
    coloring_params.set("symmetrize", symmetrize);

    testColorer.run("Test One", coloring_params);
    if (!ok) ierr++;
  }

  if (ok)
    std::cout << "TEST PASSED" << std::endl;

//Through cmake...
//Test cases -- UserInputForTests can generate Galeri or read files:
//-  tri-diagonal matrix -- can check the number of colors
//-  galeri matrix
//-  read from file:  symmetric matrix and non-symmetric matrix

//Through code ...
//Test with fitted and non-fitted maps  DONE
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
