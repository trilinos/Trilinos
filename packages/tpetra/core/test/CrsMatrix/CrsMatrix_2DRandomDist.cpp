// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This program tests matrix creation and matrix apply using matrices with
// arbitrarily distributed nonzeros (not necessarily row-based distribution).
//
// Create global matrix nonzeros randomly; store all global nonzeros on 
// each proc in a std::map.
// Create distributed vectors with randomized entries using Trilinos' default
// maps 
// For each test (linear row-wise distribution, linear column-wise distribution,
//                random distribution of nonzeros (2D) to processors)
//    distribute matrix nonzeros (each proc selects subset of global nonzeros)
//    create distributed CrsMatrix
//    perform SpMV (nMatvec SpMVs)
//    return result of SpMV
// Compare norms of results of SpMV from all distributions; they should be the
// same.
//
// NOTE:  timings are also reported but should be interpreted carefully.
// This test does not attempt to optimize the distribution of the vectors to
// minimize communication costs.  Moreover, 2D random distribution of nonzeros
// can lead to high communication volume; a better 2D approach would be a 
// block-based approach that better aligns vector entries with matrix entries.

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Class to generate, distribute and apply nonzeros 
template <typename scalar_t, typename gno_t>
class generatedNonzeros
{
public:
  using map_t = Tpetra::Map<>;
  using vector_t = Tpetra::Vector<scalar_t>;

  // Randomly generate all nonzeros for the entire matrix on every processor
  // Values are in the range 0.1 to 10.1;
  generatedNonzeros(
    size_t nRows_, size_t nCols_, size_t nnz,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_
  ) :
    nRows(nRows_), nCols(nCols_), comm(comm_)
  {
    // Synchronize RNG; want all processors to generate the same nonzeros
    srand(1);

    // use a map to remove duplicates and sort entries by row
    for (size_t n = 0; n < nnz; n++) {
      gno_t i = std::rand() % nRows;
      gno_t j = std::rand() % nCols;
      scalar_t val = 0.1 + scalar_t(std::rand() % 10);
      nzmap[std::make_pair(i,j)] = val;
    }

    nNz = nzmap.size();
  }

  // Select nonzeros from nzmap that are to be assigned to this processor
  // Create CRS data structures to ease matrix construction
  void getMyNonzeros(
    const int distribution,  // flag indicating how to distribute the nonzeros
                             //   == 1 --> row-wise (Trilinos default)
                             //   == 2 --> column-wise
                             //   == 3 --> randomly assign nonzeros to procs
    Teuchos::Array<gno_t> &rowIdx,    // output:  unique row indices of
                                      // nonzeros on this proc (sorted
                                      // in ascending order)
    Teuchos::Array<size_t> &nPerRow,  // output:  nPerRow[i] == 
                                      // number of nonzeros
                                      // in rowIdx[i] on this proc
    Teuchos::Array<size_t> &offsets,  // output:  CRS offset array;
                                      // length = length(rowIdx) + 1
    Teuchos::Array<gno_t> &colIdx,    // output:  CRS column indices
    Teuchos::Array<scalar_t> &val     // output:  nonzerovalues
  ) const
  {
    int np = comm->getSize();
    int me = comm->getRank();

    // Precompute values needed for distribution 1 (linear row-wise)
    gno_t nMyRows = nRows / np + (int(nRows % np) > me);
    gno_t myFirstRow = (me * (nRows / np) + std::min<gno_t>(nRows % np, me));
    gno_t myLastRow = myFirstRow + nMyRows - 1;

    // Precompute values needed for distribution 2 (linear column-wise)
    gno_t nMyCols = nCols / np + (int(nCols % np) > me);
    gno_t myFirstCol = (me * (nCols / np) + std::min<gno_t>(nCols % np, me));
    gno_t myLastCol = myFirstCol + nMyCols - 1;

    // Produce CRS formatted arrays of nonzeros assigned to this processor
    // given a requested Distribution
    gno_t prev_i = std::numeric_limits<gno_t>::max();

    // Loop over global nonzeros; insert those assigned to this processor in
    // CRS arrays.
    // Exploit fact that nzmap entries are sorted by row i to build CRS arrays
    for (auto nz = nzmap.begin(); nz != nzmap.end(); nz++) {

      gno_t i = nz->first.first;
      gno_t j = nz->first.second;
      scalar_t v = nz->second;

      // Check whether nonzero (i,j) should be stored on this processor
      bool mine = false;
      switch (distribution) {
      case 1:  // linear row-wise
        if (i >= myFirstRow && i <= myLastRow) mine = true;
        break;
      case 2:  // linear col-wise
        if (j >= myFirstCol && j <= myLastCol) mine = true;
        break;
      case 3:  // random
        int randomproc = std::rand() % np;
        if (me == randomproc) mine = true;
        break;
      }

      if (mine) {
        // nzmap entries are sorted by i; add a new i when different from prev i
        if (i != prev_i) {
          rowIdx.push_back(i);
          nPerRow.push_back(0);
          prev_i = i;
        }
        colIdx.push_back(j);
        val.push_back(v);
        nPerRow.back()++;
      }
    }

    // Compute prefix sum in offsets array
    offsets.resize(rowIdx.size() + 1);
    offsets[0] = 0;
    size_t nRowIdx = size_t(rowIdx.size());
    for (size_t row = 0; row < nRowIdx; row++)
      offsets[row+1] = offsets[row] + nPerRow[row];
  }

  // Distribute nonzeros to processors, create CrsMatrix, then apply it to
  // input vector x, giving y
  // Time the SpMV application
  void distributeAndApply(
    const int distribution,  // flag indicating how to distribute the nonzeros
                             //   == 1 --> row-wise (Trilinos default)
                             //   == 2 --> column-wise
                             //   == 3 --> randomly assign nonzeros to procs
    const int nMatvecs,      // Number of SpMV to do (for timing test)
    const vector_t &xvec,    // input:  domain vector
    vector_t &yvec           // output: range vector
  ) const
  {
    // Select this processor's nonzeros based on distribution
    Teuchos::Array<size_t> offsets;
    Teuchos::Array<size_t> nPerRow;
    Teuchos::Array<gno_t> rowIdx;
    Teuchos::Array<gno_t> colIdx;
    Teuchos::Array<scalar_t> val;

    getMyNonzeros(distribution, rowIdx, nPerRow, offsets, colIdx, val);
  
    // Build the CrsMatrix with the assigned nonzeros
    using matrix_t = Tpetra::CrsMatrix<scalar_t>;
  
    size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> rowMap =
             Teuchos::rcp(new map_t(dummy, rowIdx(), 0, comm));

    Teuchos::RCP<matrix_t> Amat = Teuchos::rcp(new matrix_t(rowMap, nPerRow()));
  
    size_t nRowIdx = size_t(rowIdx.size());
    for (size_t r = 0; r < nRowIdx; r++) {
      size_t tmp = offsets[r+1] - offsets[r];
      Amat->insertGlobalValues(rowIdx[r], 
                               colIdx(offsets[r],tmp), val(offsets[r],tmp));
    }
  
    std::string tname;
    {
      switch (distribution) {
        case 1:  tname = "fillComplete: 1 row-wise"; break;
        case 2:  tname = "fillComplete: 2 column-wise"; break;
        case 3:  tname = "fillComplete: 3 random 2D"; break;
      }
      auto timer = Teuchos::TimeMonitor::getNewTimer(tname);

      Teuchos::TimeMonitor tt(*timer);
      Amat->fillComplete(xvec.getMap(), yvec.getMap());
    }

    std::cout << comm->getRank() 
              << ": nRows " << Amat->getLocalNumRows()
              << "; nCols " << Amat->getLocalNumCols()
              << "; nnz " << Amat->getLocalNumEntries()
              << "; import " 
              << (Amat->getGraph()->getImporter() == Teuchos::null ? 0 :
                  Amat->getGraph()->getImporter()->getNumExportIDs())
              << "; export " 
              << (Amat->getGraph()->getExporter() == Teuchos::null ? 0 : 
                  Amat->getGraph()->getExporter()->getNumExportIDs())
              << std::endl;
  
    // Perform SpMV; do several iterations to get some timing info
    {
      switch (distribution) {
        case 1:  tname = "SpMV: 1 row-wise"; break;
        case 2:  tname = "SpMV: 2 column-wise"; break;
        case 3:  tname = "SpMV: 3 random 2D"; break;
      }
  
      scalar_t alpha = 2.;
      scalar_t beta = 3.;
      auto timer = Teuchos::TimeMonitor::getNewTimer(tname);
      for (int n = 0; n < nMatvecs; n++) {
        Teuchos::TimeMonitor tt(*timer);
        Amat->apply(xvec, yvec, Teuchos::NO_TRANS, alpha, beta);
      }
    }
  }

  // Distribute nonzeros to processors, create CrsMatrix, then apply its 
  // transpose to input vector x, giving y
  // Time the SpMV application
  void distributeAndApplyTranspose(
    const int distribution,  // flag indicating how to distribute the nonzeros
                             //   == 1 --> row-wise (Trilinos default)
                             //   == 2 --> column-wise
                             //   == 3 --> randomly assign nonzeros to procs
    const int nMatvecs,      // Number of SpMV to do (for timing test)
    vector_t &xvec,          // output:  domain vector
    const vector_t &yvec     // input: range vector
  ) const
  {
    // Select this processor's nonzeros based on distribution
    Teuchos::Array<size_t> offsets;
    Teuchos::Array<size_t> nPerRow;
    Teuchos::Array<gno_t> rowIdx;
    Teuchos::Array<gno_t> colIdx;
    Teuchos::Array<scalar_t> val;

    getMyNonzeros(distribution, rowIdx, nPerRow, offsets, colIdx, val);
  
    // Build the CrsMatrix with the assigned nonzeros
    using matrix_t = Tpetra::CrsMatrix<scalar_t>;
  
    size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> rowMap =
             Teuchos::rcp(new map_t(dummy, rowIdx(), 0, comm));

    Teuchos::RCP<matrix_t> Amat = Teuchos::rcp(new matrix_t(rowMap, nPerRow()));
  
    size_t nRowIdx = size_t(rowIdx.size());
    for (size_t r = 0; r < nRowIdx; r++) {
      size_t tmp = offsets[r+1] - offsets[r];
      Amat->insertGlobalValues(rowIdx[r], 
                               colIdx(offsets[r],tmp), val(offsets[r],tmp));
    }
  
    std::string tname;
    {
      switch (distribution) {
        case 1:  tname = "fillComplete: 1 row-wise"; break;
        case 2:  tname = "fillComplete: 2 column-wise"; break;
        case 3:  tname = "fillComplete: 3 random 2D"; break;
      }
      auto timer = Teuchos::TimeMonitor::getNewTimer(tname);

      Teuchos::TimeMonitor tt(*timer);
      Amat->fillComplete(xvec.getMap(), yvec.getMap());
    }

    std::cout << comm->getRank() 
              << ": nRows " << Amat->getLocalNumRows()
              << "; nCols " << Amat->getLocalNumCols()
              << "; nnz " << Amat->getLocalNumEntries()
              << "; import " 
              << (Amat->getGraph()->getImporter() == Teuchos::null ? 0 :
                  Amat->getGraph()->getImporter()->getNumExportIDs())
              << "; export " 
              << (Amat->getGraph()->getExporter() == Teuchos::null ? 0 : 
                  Amat->getGraph()->getExporter()->getNumExportIDs())
              << std::endl;
  
    // Perform transpose SpMV; do several iterations to get some timing info
    {
      switch (distribution) {
        case 1:  tname = "SpMV Transpose: 1 row-wise"; break;
        case 2:  tname = "SpMV Transpose: 2 column-wise"; break;
        case 3:  tname = "SpMV Transpose: 3 random 2D"; break;
      }
  
      scalar_t alpha = 2.;
      scalar_t beta = 3.;
      auto timer = Teuchos::TimeMonitor::getNewTimer(tname);
      for (int n = 0; n < nMatvecs; n++) {
        Teuchos::TimeMonitor tt(*timer);
        Amat->apply(yvec, xvec, Teuchos::TRANS, alpha, beta);
      }
    }
  }

private:

  using coord = std::pair<gno_t, gno_t>;
  struct compareCoord {  // sort nonzeros by row, then column
    bool operator() (const coord &lhs, const coord &rhs) const
    { if (lhs.first < rhs.first) return true;
      if ((lhs.first == rhs.first) && (lhs.second < rhs.second)) return true;
      return false;
    }
  };
  std::map<coord, scalar_t, compareCoord> nzmap;  // sorted global nonzeros
  
  size_t nRows, nCols, nNz;
  const Teuchos::RCP<const Teuchos::Comm<int> > comm;

};
 
////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[]) 
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using scalar_t = Tpetra::Details::DefaultTypes::scalar_type; 
  using gno_t = Tpetra::Map<>::global_ordinal_type;

  int me = comm->getRank();
  int np = comm->getSize();

  const int nMatvecs = 1000;
  size_t nRows = np * 100;
  size_t nCols = np * 200;
  size_t nNz = np * 1000;

  // Create random nonzeros -- all global nonzeros generated on every processor
  generatedNonzeros<scalar_t, gno_t> gNz(nRows, nCols, nNz, comm);

  // Create vectors; use Trilinos default range and domain maps
  // These vectors do not optimize communication for the random 2D distribution
  using map_t = Tpetra::Map<>;
  using vector_t = Tpetra::Vector<scalar_t>;

  Teuchos::RCP<const map_t> rangeMap = Teuchos::rcp(new map_t(nRows, 0, comm));
  vector_t yvec(rangeMap);

  Teuchos::RCP<const map_t> domainMap = Teuchos::rcp(new map_t(nCols, 0, comm));
  vector_t xvec(domainMap);
  xvec.randomize();

  // Row-wise 1D distribution
  yvec.putScalar(1000.);
  gNz.distributeAndApply(1, nMatvecs, xvec, yvec);
  scalar_t row1DNorm1 = yvec.norm1();
  scalar_t row1DNorm2 = yvec.norm2();
  scalar_t row1DNormInf = yvec.normInf();
  if (me == 0) 
    std::cout << "Row-wise 1D distribution:  norm1 " << row1DNorm1
              << "; norm2 " << row1DNorm2
              << "; norminf " << row1DNormInf
              << std::endl;

  // Column-wise 1D distribution
  yvec.putScalar(1000.);
  gNz.distributeAndApply(2, nMatvecs, xvec, yvec);
  scalar_t col1DNorm1 = yvec.norm1();
  scalar_t col1DNorm2 = yvec.norm2();
  scalar_t col1DNormInf = yvec.normInf();
  if (me == 0) 
    std::cout << "Col-wise 1D distribution:  norm1 " << col1DNorm1
              << "; norm2 " << col1DNorm2
              << "; norminf " << col1DNormInf
              << std::endl;

  // Random 2D distribution
  yvec.putScalar(1000.);
  gNz.distributeAndApply(3, nMatvecs, xvec, yvec);
  scalar_t random2DNorm1 = yvec.norm1();
  scalar_t random2DNorm2 = yvec.norm2();
  scalar_t random2DNormInf = yvec.normInf();
  if (me == 0) 
    std::cout << "Random 2D distribution:   norm1 " << random2DNorm1
              << "; norm2 " << random2DNorm2
              << "; norminf " << random2DNormInf
              << std::endl;

  // Check results
  int ierr = 0;

  const scalar_t epsilon = 0.0000001;
  if (std::abs(col1DNorm2 - row1DNorm2) > epsilon) {
    ierr++;
    if (me == 0) 
      std::cout << "FAIL:  column-wise 1D norm " << col1DNorm2
                << " - " << row1DNorm2 << " row-wise 1D norm" 
                << " = " << std::abs(col1DNorm2 - row1DNorm2) << std::endl;
  }

  if (std::abs(random2DNorm2 - row1DNorm2) > epsilon) {
    ierr++;
    if (me == 0) 
      std::cout << "FAIL:  random 2D norm " << random2DNorm2
                << " = " << row1DNorm2 << " row-wise 1D norm" 
                << " = " << std::abs(random2DNorm2 - row1DNorm2) << std::endl;
  }

  // Now do same with transpose
  yvec.randomize();

  // Row-wise 1D distribution
  xvec.putScalar(1000.);
  gNz.distributeAndApplyTranspose(1, nMatvecs, xvec, yvec);
  row1DNorm1 = yvec.norm1();
  row1DNorm2 = xvec.norm2();
  row1DNormInf = xvec.normInf();
  if (me == 0) 
    std::cout << "Row-wise 1D distribution Transpose:  norm1 " << row1DNorm1
              << "; norm2 " << row1DNorm2
              << "; norminf " << row1DNormInf
              << std::endl;

  // Column-wise 1D distribution
  xvec.putScalar(1000.);
  gNz.distributeAndApplyTranspose(2, nMatvecs, xvec, yvec);
  col1DNorm1 = xvec.norm1();
  col1DNorm2 = xvec.norm2();
  col1DNormInf = xvec.normInf();
  if (me == 0) 
    std::cout << "Col-wise 1D distribution Transpose:  norm1 " << col1DNorm1
              << "; norm2 " << col1DNorm2
              << "; norminf " << col1DNormInf
              << std::endl;

  // Random 2D distribution
  xvec.putScalar(1000.);
  gNz.distributeAndApplyTranspose(3, nMatvecs, xvec, yvec);
  random2DNorm1 = xvec.norm1();
  random2DNorm2 = xvec.norm2();
  random2DNormInf = xvec.normInf();
  if (me == 0) 
    std::cout << "Random 2D distribution Transpose:   norm1 " << random2DNorm1
              << "; norm2 " << random2DNorm2
              << "; norminf " << random2DNormInf
              << std::endl;

  // Check results
  if (std::abs(col1DNorm2 - row1DNorm2) > epsilon) {
    ierr++;
    if (me == 0) 
      std::cout << "FAIL:  Transpose column-wise 1D norm " << col1DNorm2
                << " - " << row1DNorm2 << " row-wise 1D norm" 
                << " = " << std::abs(col1DNorm2 - row1DNorm2) << std::endl;
  }

  if (std::abs(random2DNorm2 - row1DNorm2) > epsilon) {
    ierr++;
    if (me == 0) 
      std::cout << "FAIL:  Transpose random 2D norm " << random2DNorm2
                << " = " << row1DNorm2 << " row-wise 1D norm" 
                << " = " << std::abs(random2DNorm2 - row1DNorm2) << std::endl;
  }

  if (ierr == 0 && me == 0) 
    std::cout << "PASS" << std::endl;

  Teuchos::TimeMonitor::summarize();
  return ierr;
}

