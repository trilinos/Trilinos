// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra 
{

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class DistributionLowerTriangularBlock : public Distribution<gno_t,scalar_t> {
// Seher Acer's lower-triangular block decomposition for triangle counting
// First distribute lower-triangular entries only via 1D linear distribution
// Then redistribute according to chunk-based algorithm

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::comm;
  using Distribution<gno_t,scalar_t>::nrows;
  using typename Distribution<gno_t,scalar_t>::NZindex_t;
  using typename Distribution<gno_t,scalar_t>::LocalNZmap_t;

  DistributionLowerTriangularBlock(gno_t nrows_, 
                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm_,
                  const Teuchos::ParameterList &params) :
                  Distribution<gno_t,scalar_t>(nrows_, comm_, params),
                  initialDist(nrows_, comm_, params),
                  sortByDegree(false), redistributed(false), nChunks(0)
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

    if (me == 0) std::cout << "\n LowerTriangularBlock Distribution: "
                           << "\n     np      = " << np 
                           << "\n     nChunks = " << nChunks
                           << std::endl;
  }

  enum DistributionType DistType() { return LowerTriangularBlock; }

  Teuchos::Array<gno_t> getChunkCuts() { return chunkCuts; }

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

    Teuchos::Array<gno_t> permuteIndex;

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

      if (me == 0) {  // TODO:  do on all procs once have allgatherv

        for (size_t i = 0 ; i != nrows; i++) permuteIndex[i] = i;

        std::sort(permuteIndex.begin(), permuteIndex.end(),
                  [&](const size_t& a, const size_t& b) {
                      return (globalRowBuf[a] > globalRowBuf[b]);
                  }
        );

        // Compute inverse permutation; it is more useful for our needs

        gno_t tmp;
        for (size_t i = 0, change = 0; i < permuteIndex.size(); i++) {
          globalRowBuf[permuteIndex[i]] = i;
        }
        globalRowBuf.swap(permuteIndex);
      }

      Teuchos::broadcast<int,gno_t>(*comm, 0, permuteIndex(0,nrows)); 
                                   // Ick!  Use a directory  TODO
    }

    // Now, to determine the chunks, we care only about the number of 
    // nonzeros in the lower triangular matrix.
    // Compute nnzPerRow; distribution is currently 1D 
    nnzPerRow.assign(nMyRows, 0);
    for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
      gno_t I = (sortByDegree ? permuteIndex[it->first.first] 
                              : it->first.first);
      gno_t J = (sortByDegree ? permuteIndex[it->first.second]
                              : it->first.second);
      if (J <= I) // Lower-triangular part 
        nnzPerRow[it->first.first - myFirstRow]++;
    }

    // TODO Dunno why there isn't a Teuchos::gatherAllv; 
    // TODO for now, compute and broadcast

    Teuchos::gatherv<int,gno_t>(nnzPerRow.getRawPtr(), nMyRows,
                                globalRowBuf.getRawPtr(), 
                                rcvcnt.getRawPtr(), disp.getRawPtr(), 
                                0, *comm);

    Teuchos::Array<int>().swap(rcvcnt);  // no longer needed
    Teuchos::Array<int>().swap(disp);    // no longer needed

    size_t nnz = localNZ.size(), gNnz;
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
        targetRunningTotal += (target * (chunkCnt+1));
        while (I < nrows) {
          size_t nextNnz = (sortByDegree ? globalRowBuf[permuteIndex[I]]
                                         : globalRowBuf[I]);
          if (currentRunningTotal + nextNnz < targetRunningTotal) {
            currentRunningTotal += nextNnz;
            I++;
          }
          else
            break;
        } 
        chunkCuts[chunkCnt+1] = I;
      }
      chunkCuts[nChunks] = nrows;
    }

    // Free memory associated with globalRowBuf
    Teuchos::Array<gno_t>().swap(globalRowBuf);

    Teuchos::broadcast<int,gno_t>(*comm, 0, chunkCuts(0,nChunks+1));

    //std::cout << comm->getRank() << " KDDKDD chunkCuts: ";
    //for (int kdd=0; kdd <= nChunks; kdd++) std::cout << chunkCuts[kdd] << " ";
    //std::cout << std::endl;

    // Determine new owner of each nonzero; buffer for sending
    Teuchos::Array<gno_t> iOut(localNZ.size());
    Teuchos::Array<gno_t> jOut(localNZ.size());
    Teuchos::Array<scalar_t> vOut(localNZ.size());
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

        //std::cout << comm->getRank() 
        //          << "    KDDKDD IJ (" 
        //          << it->first.first << "," << it->first.second
        //          << ") permuted to ("
        //          << iOut[sendCnt] << "," << jOut[sendCnt] 
        //          << ") being sent to " << pOut[sendCnt]
        //          << std::endl;
        sendCnt++;
      }
    }

    // Free memory associated with localNZ and permuteIndex
    LocalNZmap_t().swap(localNZ);
    if (sortByDegree) Teuchos::Array<gno_t>().swap(permuteIndex);

    // Use a Distributor to send nonzeros to new processors.
    Tpetra::Distributor plan(comm);
    size_t nrecvs = plan.createFromSends(pOut(0,sendCnt));
    Teuchos::Array<gno_t> iIn(nrecvs);
    Teuchos::Array<gno_t> jIn(nrecvs);
    Teuchos::Array<scalar_t> vIn(nrecvs);

    // TODO:  With more clever packing, could do only one round of communication
    plan.doPostsAndWaits<gno_t>(iOut(0,sendCnt), 1, iIn());
    plan.doPostsAndWaits<gno_t>(jOut(0,sendCnt), 1, jIn());
    plan.doPostsAndWaits<scalar_t>(vOut(0,sendCnt), 1, vIn());

    // Put received nonzeros in map
    for (size_t n = 0; n < nrecvs; n++) {
      NZindex_t nz(iIn[n], jIn[n]);
      localNZ[nz] = vIn[n];
    }

    redistributed = true;
  }

private:
  // Initially distribute nonzeros with a 1D linear distribution
  Distribution1DLinear<gno_t,scalar_t> initialDist;  

  // Flag indicating whether matrix should be reordered and renumbered 
  // in decreasing sort order of number of nonzeros per row in full matrix
  bool sortByDegree;

  // Flag whether redistribution has occurred yet
  bool redistributed;  

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

  LowerTriangularBlockOperator(
    const Teuchos::RCP<const matrix_t> &lowerTriangularMatrix_) 
  : lowerTriangularMatrix(lowerTriangularMatrix_)
  {
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

    // Multiply lower triangular
    lowerTriangularMatrix->apply(x, y, Teuchos::NO_TRANS, alpha, beta);

    // Multiply upper triangular
    lowerTriangularMatrix->apply(x, y, Teuchos::TRANS, alpha, ONE);

    // Subtract out duplicate diagonal terms
    y.elementWiseMultiply(-ONE, *diag, x, ONE);
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
  Teuchos::RCP<vector_t> diag;
};


}
#endif
