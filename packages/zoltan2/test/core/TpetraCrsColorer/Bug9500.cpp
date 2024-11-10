// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Kokkos_Random.hpp"
#include "Zoltan2_TestHelpers.hpp"
#include "Zoltan2_TpetraCrsColorer.hpp"

// Class to test the Colorer utility
class ColorerTest {
public:
  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using graph_t = Tpetra::CrsGraph<>;
  using matrix_t = Tpetra::CrsMatrix<zscalar_t>;
  using multivector_t = Tpetra::MultiVector<zscalar_t>;
  using execution_space_t = typename matrix_t::device_type::execution_space;

  ///////////////////////////////////////////////////////////
  // Construct the test:  

  ColorerTest(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, int multiple)
  {
    int me = comm->getRank();
    int np = comm->getSize();

    // Create non-symmetrix matrix with non-contiguous row map -- only even GIDs
    size_t myNrows = 4;
    Teuchos::Array<gno_t> myRows(myNrows);
    for (size_t i = 0; i < myNrows; i++) {
      myRows[i] = multiple * (me * myNrows + i);  
    }

    Tpetra::global_size_t dummy =
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, myRows, 0, comm));

    size_t nnz = 2;
    JBlock = rcp(new matrix_t(map, nnz));
    Teuchos::Array<gno_t> myCols(nnz);
    Teuchos::Array<zscalar_t> myVals(nnz);

    for (size_t i = 0; i < myNrows; i++) {
      auto gid = map->getGlobalElement(i);
      size_t cnt = 0;
      myCols[cnt++] = gid;
      if (gid+multiple <= map->getMaxAllGlobalIndex()) 
        myCols[cnt++] = gid+multiple;
      JBlock->insertGlobalValues(gid, myCols(0,cnt), myVals(0, cnt));
    }
    JBlock->fillComplete();

    // Fill JBlock with random numbers for a better test.
    using IST = typename Kokkos::ArithTraits<zscalar_t>::val_type;
    using pool_type = 
          Kokkos::Random_XorShift64_Pool<execution_space_t>;
    pool_type rand_pool(static_cast<uint64_t>(me));

    Kokkos::fill_random(JBlock->getLocalMatrixDevice().values, rand_pool, 
                        static_cast<IST>(1.), static_cast<IST>(9999.));

    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));
    JBlock->describe(foo, Teuchos::VERB_EXTREME);

    // Build same matrix with cyclic domain and range maps
    // To make range and domain maps differ for square matrices,
    // keep some processors empty in the cyclic maps

    size_t nIndices = std::max(JBlock->getGlobalNumCols(),
                               JBlock->getGlobalNumRows());
    Teuchos::Array<gno_t> indices(nIndices);

    Teuchos::RCP<const map_t> vMapCyclic = 
                 getCyclicMap(JBlock->getGlobalNumCols(), indices, np-1, 
                              multiple, comm);
    Teuchos::RCP<const map_t> wMapCyclic =
                 getCyclicMap(JBlock->getGlobalNumRows(), indices, np-2, 
                              multiple, comm);

    // Make JCyclic:  same matrix with different Domain and Range maps
    RCP<const graph_t> block_graph = JBlock->getCrsGraph();
    RCP<graph_t> cyclic_graph = rcp(new graph_t(*block_graph));
    cyclic_graph->resumeFill();
    cyclic_graph->fillComplete(vMapCyclic, wMapCyclic);
    JCyclic = rcp(new matrix_t(cyclic_graph));
    JCyclic->resumeFill();
    TEUCHOS_ASSERT(block_graph->getLocalNumRows() == 
                   cyclic_graph->getLocalNumRows());
    {
      auto val_s = JBlock->getLocalMatrixHost().values;
      auto val_d = JCyclic->getLocalMatrixHost().values;
      TEUCHOS_ASSERT(val_s.extent(0) == val_d.extent(0));
      Kokkos::deep_copy(val_d, val_s);
    }
    JCyclic->fillComplete();
    JCyclic->describe(foo, Teuchos::VERB_EXTREME);
  }

  ////////////////////////////////////////////////////////////////
  bool run(const char* testname, Teuchos::ParameterList &params) {

    bool ok = true;

    params.set("symmetric", false);

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
    Jp->setAllToScalar(static_cast<zscalar_t>(-1.));

    colorer.reconstructMatrix(W, *Jp);

    // Check that values of J = values of Jp
    auto J_local_matrix = J->getLocalMatrixDevice();
    auto Jp_local_matrix = Jp->getLocalMatrixDevice();
    const size_t num_local_nz = J->getLocalNumEntries();

    Kokkos::parallel_reduce(
      "TpetraCrsColorer::testReconstructedMatrix()",
      Kokkos::RangePolicy<execution_space_t>(0, num_local_nz),
      KOKKOS_LAMBDA(const size_t nz, int &errorcnt) {
        if (J_local_matrix.values(nz) != Jp_local_matrix.values(nz)) {
          Kokkos::printf("Error in nonzero comparison %zu:  %g != %g",
                  nz, J_local_matrix.values(nz), Jp_local_matrix.values(nz));
          errorcnt++;
        }
      }, 
      ierr);
   

    if (ierr > 0) {
      std::cout << testname << " FAILED on rank " << me << " with "
                << (useBlock ? "Block maps" : "Cyclic maps")
                << std::endl;
      params.print();
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
    int multiple,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  {
    size_t cnt = 0;
    int me = comm->getRank();
    int np = comm->getSize();
    if (mapNumProc > np) mapNumProc = np; // corner case: bad input
    if (mapNumProc <= 0) mapNumProc = 1;  // corner case: np is too small

    for (size_t i = 0; i < nIndices; i++) 
      if (me == int(i % np)) indices[cnt++] = multiple*i; 

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
int doTheTest(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, int multiple)
{
  bool ok = true;
  int ierr = 0;

  ColorerTest testColorer(comm, multiple);

  // Set parameters and compute coloring
  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";
    bool symmetrize = true;  // Zoltan's TRANSPOSE symmetrization, if needed

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetrize", symmetrize);

    ok = testColorer.run("Test One", coloring_params);
    if (!ok) ierr++;
  }

  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";
    bool symmetrize = false;  // colorer's BIPARTITE symmetrization, if needed

    coloring_params.set("MatrixType", matrixType);
    coloring_params.set("symmetrize", symmetrize);

    ok = testColorer.run("Test Two", coloring_params);
    if (!ok) ierr++;
  }

  {
    Teuchos::ParameterList coloring_params;
    std::string matrixType = "Jacobian";

    coloring_params.set("MatrixType", matrixType);

    ok = testColorer.run("Test Three", coloring_params);
    if (!ok) ierr++;
  }
  return ierr;
}

///////////////////////////////////////////////////////////////////////////////
int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int ierr = 0;

  ierr += doTheTest(comm, 1);  // Contiguous row map
  ierr += doTheTest(comm, 2);  // Non-contiguous row map -- 
                               // only even-numbered indices
  ierr += doTheTest(comm, 5);  // Indices spaced wider than rows/proc

  int gerr;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (comm->getRank() == 0) {
    if (gerr == 0)
      std::cout << "TEST PASSED" << std::endl;
    else
      std::cout << "TEST FAILED" << std::endl;
  }
}
