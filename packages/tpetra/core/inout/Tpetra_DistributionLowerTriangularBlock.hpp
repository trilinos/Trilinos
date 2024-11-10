// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TPETRA_DISTRIBUTORLOWERTRIANGULARBLOCK_HPP
#define __TPETRA_DISTRIBUTORLOWERTRIANGULARBLOCK_HPP

// Needed by DistributionLowerTriangularBlock
#include "Tpetra_Distributor.hpp"

// Needed by LowerTriangularBlock operator
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra 
{

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class DistributionLowerTriangularBlock : public Distribution<gno_t,scalar_t> {
// Seher Acer's lower-triangular block decomposition for triangle counting
// See also:  LowerTriangularBlockOperator below that allows this distribution
// to be used in Tpetra SpMV.
// 
// Requirements:
//    Matrix must be square (undirected graph)
//    Number of processors np = q(q+1)/2 for some q.
// 
// Only the lower triangular portion of the matrix is stored.
// Processors are arranged logically as follows:
//    0
//    1  2
//    3  4  5
//    ...
//
// The lower triangular part of the matrix is stored in a 2D distribution.
// For example, the dense 7x7 lower triangular matrix below would be assigned
// to processors according to numbers shown as follows:
//    0 |   |
//    00|   |
//    ---------
//    11|2  |
//    11|22 |
//    11|222|
//    ---------
//    33|444|5
//    33|444|55
//    ...
// (Note that we expect the matrix to be sparse.  For dense matrices,
// CrsMatrix is the wrong tool.)
//    
// Matrix rows are assigned to processor rows greedily to roughly balance 
//   (# nonzeros in processor row / # processors in processor row)
// across processor rows.
// The same cuts are used to divide rows and columns among processors
// (that is, all processors have a square block).
// 
// The lower triangular algorithm:
// 1. distribute all matrix entries via 1D linear distribution
//    (this initial distribution is needed to avoid storing the entire 
//    matrix on one processor, while providing info about the nonzeros per row
//    needed in step 2.
// 2. (optional) sort rows in decreasing order wrt the number of nonzeros 
//    per row
// 3. find "chunk cuts":  divisions in row assignments such that 
//    (# nonzeros in processor row / # processors in processor row) is 
//    roughly equal for all processor rows
// 4. send nonzeros to their new processor assignment
// 
// Known issues:  (TODO)
// - The sorting in Step 2 and computation of chunk cuts in step 3 are 
//   currently done in serial and requires O(number of rows) storage each 
//   processor.  More effort could parallelize this computation, but parallel
//   load balancing algorithms are more appropriate in Zoltan2 than Tpetra.
// - The sorting in Step 2 renumbers the rows (assigns new Global Ordinals to 
//   the rows) to make them contiguous, as needed in Acer's triangle counting 
//   algorithm.  
//   (Acer's algorithm relies on local indexing from the chunk boundaries to
//   find neighbors needed for communication.)
//   The class currently provides a permutation matrix P describing the 
//   reordering.  Thus, the matrix stored in the lower triangular block
//   distribution is actually P A P -- the row and column permutation of 
//   matrix A in the Matrix Market file.
//   The fact that a permuted matrix is stored complicates use of the matrix
//   in algorithms other than Acer's triangle counting.  For SpMV with the
//   vector numbered according to the MatrixMarket numbering, for example,
//   P^T must be applied to the vector before SpMV, and P^T must be applied to
//   the result of SpMV.  See LowerTriangularBlockOperator to see how this
//   permutation matrix is used.
//
// Before addressing these issues, we will decide (TODO)
// -  Is this Distribution general enough to be in Tpetra? 
// -  Should we, instead, have a separate package for distributions (that could
//    use Zoltan2 and Tpetra without circular dependence)? 
// -  Or should we allow users (such as the triangle counting algorithm) to 
//    provide their own distributions (e.g., LowerTriangularBlock) that 
//    inherit from Tpetra's Distribution class?
// For now, we will push this Distribution into Tpetra, but we will revisit 
// this decision.

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::comm;
  using Distribution<gno_t,scalar_t>::nrows;
  using typename Distribution<gno_t,scalar_t>::NZindex_t;
  using typename Distribution<gno_t,scalar_t>::LocalNZmap_t;

  using map_t = Tpetra::Map<>;
  using matrix_t = Tpetra::CrsMatrix<scalar_t>;

  DistributionLowerTriangularBlock(size_t nrows_, 
                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm_,
                  const Teuchos::ParameterList &params) :
                  Distribution<gno_t,scalar_t>(nrows_, comm_, params),
                  initialDist(nrows_, comm_, params),
                  sortByDegree(false), permMatrix(Teuchos::null),
                  redistributed(false), chunksComputed(false), nChunks(0)
  {
    int npnp = 2 * np;
    nChunks = int(std::sqrt(float(npnp)));
    while (nChunks * (nChunks + 1) < npnp) nChunks++;

    TEUCHOS_TEST_FOR_EXCEPTION(nChunks * (nChunks+1) != npnp, std::logic_error,
                               "Number of processors np = " << np <<
                               " must satisfy np = q(q+1)/2 for some q" <<
                               " for LowerTriangularBlock distribution");
    nChunksPerRow = double(nChunks) / double(nrows);

    const Teuchos::ParameterEntry *pe = params.getEntryPtr("sortByDegree");
    if (pe != NULL) sortByDegree = pe->getValue<bool>(&sortByDegree);

    pe = params.getEntryPtr("readPerProcess");
    if (pe != NULL) redistributed = pe->getValue<bool>(&redistributed);

    if (me == 0) std::cout << "\n LowerTriangularBlock Distribution: "
                           << "\n     np      = " << np 
                           << "\n     nChunks = " << nChunks
                           << std::endl;
  }

  enum DistributionType DistType() { return LowerTriangularBlock; }

  bool areChunksComputed() {return chunksComputed; }

  Teuchos::Array<gno_t> getChunkCuts() { 
    if(chunksComputed)
      return chunkCuts; 
    else {
      throw std::runtime_error("Error:  Requested chunk cuts have not been computed yet.");
    }
  }

  // Return whether this rank owns vector entry i.
  // TODO:  for now, use same vector dist as 1DLinear;
  // TODO:  think about best distribution of Vectors
  inline bool VecMine(gno_t i) { return initialDist.VecMine(i); }

  // Return whether this rank owns nonzero (i,j)
  // Vector map and row map are the same in 1D distribution.
  // But keep only the lower Triangular entries
  bool Mine(gno_t i, gno_t j) {
    if (redistributed) {
      if (j > i) return false;  // Don't keep any upper triangular entries
      else return (procFromChunks(i,j) == me);
    }
    else
      return initialDist.Mine(i,j);
  }

  inline bool Mine(gno_t i, gno_t j, int p) {return Mine(i,j);}

  // How to redistribute according to chunk-based row distribution
  void Redistribute(LocalNZmap_t &localNZ)
  {
    // Going to do chunking and sorting serially for now; 
    // need to gather per-row information from each processor
    // TODO:  think about a parallel implementation

    gno_t myFirstRow = initialDist.getFirstRow(me);
    gno_t nMyRows = initialDist.getNumRow(me);
    Teuchos::Array<gno_t> nnzPerRow(nMyRows, 0);

    Teuchos::Array<int> rcvcnt(np);
    Teuchos::Array<int> disp(np);
    for (int sum = 0, p = 0; p < np; p++) {
      int prows = initialDist.getNumRow(p);
      rcvcnt[p] = prows;
      disp[p] = sum;
      sum += prows;
    }

    // If desire sortByDegree, first need to sort with respect to ALL entries
    // in matrix (not lower-triangular entries);
    // decreasing sort by number of entries per row in global matrix.
    // Generate permuteIndex for the sorted rows

    Teuchos::Array<gno_t> permuteIndex;  // This is the inverse permutation
    Teuchos::Array<gno_t> sortedOrder;   // This is the original permutation
 
    Teuchos::Array<gno_t> globalRowBuf;
    // TODO Dunno why there isn't a Teuchos::gatherAllv; 
    // TODO for now, compute and broadcast
    if (me == 0) { 
      globalRowBuf.resize(nrows, 0);  // TODO:  Ick!  Need parallel
    }

    if (sortByDegree) {
      // Compute nnzPerRow; distribution is currently 1D and includes all nz
      for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
        gno_t I = it->first.first;
        nnzPerRow[I-myFirstRow]++;
      }

      Teuchos::gatherv<int,gno_t>(nnzPerRow.getRawPtr(), nMyRows,
                                  globalRowBuf.getRawPtr(), 
                                  rcvcnt.getRawPtr(), disp.getRawPtr(), 
                                  0, *comm);

      permuteIndex.resize(nrows); // TODO:  Ick!  Need parallel
      sortedOrder.resize(nrows); // TODO:  Ick!  Need parallel

      if (me == 0) {  // TODO:  do on all procs once have allgatherv

        for (size_t i = 0 ; i != nrows; i++) sortedOrder[i] = i;

        std::sort(sortedOrder.begin(), sortedOrder.end(),
                  [&](const size_t& a, const size_t& b) {
                      return (globalRowBuf[a] > globalRowBuf[b]);
                  }
        );

        // Compute inverse permutation; it is more useful for our needs
        for (size_t i = 0; i < nrows; i++) {
          permuteIndex[sortedOrder[i]] = i;
        }
      }

      Teuchos::broadcast<int,gno_t>(*comm, 0, permuteIndex(0,nrows)); 
                                   // Ick!  Use a directory  TODO

      // Sorting is changing the global IDs associated
      // with rows/columns.  To make this distribution applicable beyond
      // triangle counting (e.g., in a Tpetra operator), we need a way 
      // to map from the original global IDs and back again.  
      // Create a permutation matrix for use in the operator; use
      // default Tpetra layout.
      Teuchos::Array<gno_t> myRows;
      for (size_t i = 0; i < nrows; i++) {
        if (VecMine(i)) myRows.push_back(i);
      }

      Tpetra::global_size_t dummy = 
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      Teuchos::RCP<const map_t> permMap = 
                                rcp(new map_t(dummy, myRows(), 0, comm));

      permMatrix = rcp(new matrix_t(permMap, 1)); // one nz / row in permMatrix

      Teuchos::Array<gno_t> cols(1);
      Teuchos::Array<scalar_t> vals(1); vals[0] = 1.;

      for (size_t i = 0; i < permMap->getLocalNumElements(); i++) {
        gno_t gid = permMap->getGlobalElement(i);
        cols[0] = permuteIndex[gid];
        permMatrix->insertGlobalValues(gid, cols(), vals());
      }

      permMatrix->fillComplete(permMap, permMap);
    }

    // Now, to determine the chunks, we care only about the number of 
    // nonzeros in the lower triangular matrix.
    // Compute nnzPerRow; distribution is currently 1D 
    nnzPerRow.assign(nMyRows, 0);
    size_t nnz = 0;
    for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
      gno_t I = (sortByDegree ? permuteIndex[it->first.first] 
                              : it->first.first);
      gno_t J = (sortByDegree ? permuteIndex[it->first.second]
                              : it->first.second);
      if (J <= I) {// Lower-triangular part 
        nnzPerRow[it->first.first - myFirstRow]++;
	nnz++;
      }
    }

    // TODO Dunno why there isn't a Teuchos::gatherAllv; 
    // TODO for now, compute and broadcast

    Teuchos::gatherv<int,gno_t>(nnzPerRow.getRawPtr(), nMyRows,
                                globalRowBuf.getRawPtr(), 
                                rcvcnt.getRawPtr(), disp.getRawPtr(), 
                                0, *comm);

    Teuchos::Array<int>().swap(rcvcnt);  // no longer needed
    Teuchos::Array<int>().swap(disp);    // no longer needed

    size_t gNnz;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &nnz, &gNnz);

    chunkCuts.resize(nChunks+1, 0);


    if (me == 0) { // TODO:  when have allgatherv, can do on all procs
                   // TODO:  or better, implement parallel version

      // Determine chunk cuts
      size_t target = gNnz / np;  // target nnz per processor
      size_t targetRunningTotal = 0;
      size_t currentRunningTotal = 0;
      gno_t I = gno_t(0);
      for (int chunkCnt = 0; chunkCnt < nChunks; chunkCnt++) {
        targetRunningTotal = (target * (chunkCnt+1));
	currentRunningTotal = 0;
        while (I < static_cast<gno_t>(nrows)) {
          size_t nextNnz = (sortByDegree ? globalRowBuf[sortedOrder[I]]
                                         : globalRowBuf[I]);
          if (currentRunningTotal + nextNnz <= targetRunningTotal) {
            currentRunningTotal += nextNnz;
            I++;
          }
          else
            break;
        } 
        chunkCuts[chunkCnt+1] = I;
      }
      chunkCuts[nChunks] = static_cast<gno_t>(nrows);
    }

    // Free memory associated with globalRowBuf
    Teuchos::Array<gno_t>().swap(globalRowBuf);

    Teuchos::broadcast<int,gno_t>(*comm, 0, chunkCuts(0,nChunks+1));
    chunksComputed = true;

    // Determine new owner of each nonzero; buffer for sending
    Kokkos::View<gno_t*, Kokkos::HostSpace> iOut("iOut", localNZ.size());
    Kokkos::View<gno_t*, Kokkos::HostSpace> jOut("jOut", localNZ.size());
    Kokkos::View<scalar_t*, Kokkos::HostSpace> vOut("vOut", localNZ.size());
    Teuchos::Array<int> pOut(localNZ.size());

    size_t sendCnt = 0;
    for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
      iOut[sendCnt] = (sortByDegree ? permuteIndex[it->first.first] 
                                    : it->first.first);
      jOut[sendCnt] = (sortByDegree ? permuteIndex[it->first.second]
                                    : it->first.second);
      if (jOut[sendCnt] <= iOut[sendCnt]) { // keep only lower diagonal entries
        vOut[sendCnt] = it->second;
        pOut[sendCnt] = procFromChunks(iOut[sendCnt], jOut[sendCnt]);

        sendCnt++;
      }
    }

    // Free memory associated with localNZ and permuteIndex
    LocalNZmap_t().swap(localNZ);
    if (sortByDegree) Teuchos::Array<gno_t>().swap(permuteIndex);

    // Use a Distributor to send nonzeros to new processors.
    Tpetra::Distributor plan(comm);
    size_t nrecvs = plan.createFromSends(pOut(0,sendCnt));
    Kokkos::View<gno_t*, Kokkos::HostSpace> iIn("iIn", nrecvs);
    Kokkos::View<gno_t*, Kokkos::HostSpace> jIn("jIn", nrecvs);
    Kokkos::View<scalar_t*, Kokkos::HostSpace> vIn("vIn", nrecvs);

    // TODO:  With more clever packing, could do only one round of communication
    auto sendIndices = std::make_pair(static_cast<size_t>(0), sendCnt);
    plan.doPostsAndWaits(Kokkos::subview(iOut, sendIndices), 1, iIn);
    plan.doPostsAndWaits(Kokkos::subview(jOut, sendIndices), 1, jIn);
    plan.doPostsAndWaits(Kokkos::subview(vOut, sendIndices), 1, vIn);

    // Put received nonzeros in map
    for (size_t n = 0; n < nrecvs; n++) {
      NZindex_t nz(iIn[n], jIn[n]);
      localNZ[nz] = vIn[n];
    }

    redistributed = true;
  }

  Teuchos::RCP<matrix_t> getPermutationMatrix() const { return permMatrix; } 

private:
  // Initially distribute nonzeros with a 1D linear distribution
  Distribution1DLinear<gno_t,scalar_t> initialDist;  

  // Flag indicating whether matrix should be reordered and renumbered 
  // in decreasing sort order of number of nonzeros per row in full matrix
  bool sortByDegree;

  // Column permutation matrix built only when sortByDegree = true;
  Teuchos::RCP<matrix_t> permMatrix;

  // Flag whether redistribution has occurred yet
  // This is true 
  //   i) after Tpetra performs the redistribution or 
  //  ii) when Tpetra reads already-distributed nonzeros by readPerProcess function
  bool redistributed;  

  // If we read the already-distributed nonzeros from per-process files, 
  //  this will remain false until a triangle counting code actually computes 
  //  the chunks when the need arises.
  bool chunksComputed;

  int nChunks;  // in np = q(q+1)/2 layout, nChunks = q
  double nChunksPerRow;
  Teuchos::Array<gno_t> chunkCuts;

  int findIdxInChunks(gno_t I) {
    int m = I * nChunksPerRow;
    while (I < chunkCuts[m]) m--;
    while (I >= chunkCuts[m+1]) m++;
    return m;
  }

  int procFromChunks(gno_t I, gno_t J) {
    int m = findIdxInChunks(I);
    int n = findIdxInChunks(J);
    int p = m*(m+1)/2 + n; 
    return p;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Tpetra::Operator that works with the DistributionLowerTriangularBlock

template <typename scalar_t, 
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class LowerTriangularBlockOperator : 
      public Tpetra::Operator<scalar_t, Tpetra::Map<>::local_ordinal_type,
                                        Tpetra::Map<>::global_ordinal_type,
                                        Node> 
{
public:
  using lno_t = Tpetra::Map<>::local_ordinal_type;
  using gno_t = Tpetra::Map<>::global_ordinal_type;
  using map_t = Tpetra::Map<>;
  using import_t = Tpetra::Import<>;
  using export_t = Tpetra::Export<>;
  using vector_t = Tpetra::Vector<scalar_t>;
  using mvector_t = Tpetra::MultiVector<scalar_t>;
  using matrix_t = Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, Node>;
  using dist_t = Tpetra::DistributionLowerTriangularBlock<gno_t, scalar_t>;

  LowerTriangularBlockOperator(
    const Teuchos::RCP<const matrix_t> &lowerTriangularMatrix_,
    const dist_t &dist)
  : lowerTriangularMatrix(lowerTriangularMatrix_),
    permMatrix(dist.getPermutationMatrix())
  {
    // LowerTriangularBlockOperator requires the range map and domain map 
    // to be the same.  Check it here.
    TEUCHOS_TEST_FOR_EXCEPTION(
      !lowerTriangularMatrix->getRangeMap()->isSameAs(
                             *lowerTriangularMatrix->getDomainMap()), 
       std::logic_error,
       "The Domain and Range maps of the LowerTriangularBlock matrix "
       "must be the same");
    
    // Extract diagonals

    vector_t diagByRowMap(lowerTriangularMatrix->getRowMap());
    lowerTriangularMatrix->getLocalDiagCopy(diagByRowMap);
    diag = Teuchos::rcp(new vector_t(lowerTriangularMatrix->getRangeMap()));
    Tpetra::Export<> exporter(lowerTriangularMatrix->getRowMap(), 
                              lowerTriangularMatrix->getRangeMap());
    diag->doExport(diagByRowMap, exporter, Tpetra::ADD);
  }

  void apply(const mvector_t &x, mvector_t &y, Teuchos::ETransp mode,
             scalar_t alpha, scalar_t beta) const
  {
    scalar_t ZERO =  Teuchos::ScalarTraits<scalar_t>::zero();
    scalar_t ONE =  Teuchos::ScalarTraits<scalar_t>::one();
    if (alpha == ZERO) {
      if (beta == ZERO) y.putScalar(ZERO);
      else y.scale(beta);
      return;
    }

    if (permMatrix == Teuchos::null) {

      // Multiply lower triangular
      lowerTriangularMatrix->apply(x, y, Teuchos::NO_TRANS, alpha, beta);

      // Multiply upper triangular
      lowerTriangularMatrix->apply(x, y, Teuchos::TRANS, alpha, ONE);

      // Subtract out duplicate diagonal terms
      y.elementWiseMultiply(-alpha, *diag, x, ONE);
    }
    else {

      // With sorting, the LowerTriangularBlock distribution stores (P^T A P)
      // in the CrsMatrix, for permutation matrix P. 
      // Thus, apply must compute
      // y = P (beta (P^T y) + alpha (P^T A P) (P^T x))

      vector_t xtmp(x.getMap(), x.getNumVectors());
      vector_t ytmp(y.getMap(), y.getNumVectors());

      permMatrix->apply(x, xtmp, Teuchos::TRANS);
      if (beta != ZERO) permMatrix->apply(y, ytmp, Teuchos::TRANS);

      // Multiply lower triangular
      lowerTriangularMatrix->apply(xtmp, ytmp, Teuchos::NO_TRANS, alpha, beta);

      // Multiply upper triangular
      lowerTriangularMatrix->apply(xtmp, ytmp, Teuchos::TRANS, alpha, ONE);

      // Subtract out duplicate diagonal terms
      ytmp.elementWiseMultiply(-alpha, *diag, xtmp, ONE);

      permMatrix->apply(ytmp, y, Teuchos::NO_TRANS);
    }
  }

  Teuchos::RCP<const map_t> getDomainMap() const {
    return lowerTriangularMatrix->getDomainMap();
  }

  Teuchos::RCP<const map_t> getRangeMap() const {
    return lowerTriangularMatrix->getRangeMap();
  }

  bool hasTransposeApply() const {return true;}  // Symmetric matrix

private:
  const Teuchos::RCP<const matrix_t > lowerTriangularMatrix;
  const Teuchos::RCP<const matrix_t > permMatrix;
  Teuchos::RCP<vector_t> diag;
};


}
#endif
