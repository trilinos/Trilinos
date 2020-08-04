#ifndef __TPETRA_DISTRIBUTORLOWERTRIANGULARBLOCK_HPP
#define __TPETRA_DISTRIBUTORLOWERTRIANGULARBLOCK_HPP

#include "Tpetra_Distributor.hpp"

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
                  redistributed(false), nChunks(0)
  {
    int npnp = 2 * np;
    nChunks = int(std::sqrt(float(npnp)));
    while (nChunks * (nChunks + 1) < npnp) nChunks++;

    TEUCHOS_TEST_FOR_EXCEPTION(nChunks * (nChunks+1) != npnp, std::logic_error,
                               "Number of processors np = " << np <<
                               " must satisfy np = q(q+1)/2 for some q" <<
                               " for LowerTriangularBlock distribution");
    nChunksPerRow = double(nChunks) / double(nrows);

    if (me == 0) std::cout << "\n LowerTriangularBlock Distribution: "
                           << "\n     np      = " << np 
                           << "\n     nChunks = " << nChunks
                           << std::endl;
  }

  enum DistributionType DistType() { return LowerTriangularBlock; }

  // Return whether this rank owns vector entry i.
  // TODO:  for now, use same vector dist as 1DLinear;
  // TODO:  think about best distribution of Vectors
  inline bool VecMine(gno_t i) { return initialDist.VecMine(i); }

  // Return whether this rank owns nonzero (i,j)
  // Vector map and row map are the same in 1D distribution.
  // But keep only the lower Triangular entries
  bool Mine(gno_t i, gno_t j) {
    if (j > i) return false;  // Don't keep any upper triangular entries
    if (redistributed)
      return (procFromChunks(i,j) == me);
    else
      return initialDist.Mine(i,j);
  }

  inline bool Mine(gno_t i, gno_t j, int p) {return Mine(i,j);}

  // How to redistribute according to chunk-based row distribution
  void Redistribute(LocalNZmap_t &localNZ)
  {
std::cout << comm->getRank() << " KDDKDD begin Redistribute " << std::endl;
    // Compute nnzPerRow; distribution is currently 1D and lower triangular
    // Exploit fact that map has entries sorted by I, then J
    // Simultaneously, store everything in buffers for communication
    gno_t myFirstRow = initialDist.getFirstRow(me);
    gno_t nMyRows = initialDist.getNumRow(me);
    Teuchos::Array<size_t> nnzPerRow(nMyRows, 0);

    gno_t prevI = std::numeric_limits<gno_t>::max();
    for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
      gno_t I = it->first.first;
      nnzPerRow[I-myFirstRow]++;
    }

    // TODO For now, determine the chunks serially; can parallelize later
    {

    Teuchos::Array<int> rcvcnt(np);
    Teuchos::Array<int> disp(np);
    for (int sum = 0, p = 0; p < np; p++) {
      int prows = initialDist.getNumRow(p);
      rcvcnt[p] = prows;
      disp[p] = sum;
      sum += prows;
    }

    Teuchos::Array<size_t> globalNnzPerRow;
    if (me == 0) {
      globalNnzPerRow.resize(nrows, 0);  // TODO:  Ick!  Need parallel
    }

    // TODO Dunno why there isn't a Teuchos::gatherAllv; 
    // TODO for now, compute and broadcast
    Teuchos::gatherv<int,size_t>(nnzPerRow.getRawPtr(), nMyRows,
                                 globalNnzPerRow.getRawPtr(), 
                                 rcvcnt.getRawPtr(), disp.getRawPtr(), 
                                 0, *comm);

    size_t nnz = localNZ.size(), gNnz;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &nnz, &gNnz);

    
    // Determine chunk cuts
    chunkCuts.resize(nChunks+1, 0);
    if (me == 0) { // TODO:  when have allgatherv, can do on all procs
      size_t target = gNnz / np;  // target nnz per processor
      size_t targetRunningTotal = 0;
      size_t currentRunningTotal = 0;
      gno_t I = gno_t(0);
      for (int chunkCnt = 0; chunkCnt < nChunks; chunkCnt++) {
        targetRunningTotal += (target * (chunkCnt+1));
        while (I < nrows) {
          if (currentRunningTotal + globalNnzPerRow[I] < targetRunningTotal) {
            currentRunningTotal += globalNnzPerRow[I++];
          }
          else
            break;
        } 
        chunkCuts[chunkCnt+1] = I;
      }
      chunkCuts[nChunks] = nrows;
    }

    }

    Teuchos::broadcast<int,gno_t>(*comm, 0, chunkCuts(0,nChunks+1));

    std::cout << comm->getRank() << " KDDKDD chunkCuts: ";
    for (int kdd = 0; kdd <= nChunks; kdd++) std::cout << chunkCuts[kdd] << " ";
    std::cout << std::endl;

    // Determine new owner of each nonzero; buffer for sending
    Teuchos::Array<gno_t> iOut(localNZ.size());
    Teuchos::Array<gno_t> jOut(localNZ.size());
    Teuchos::Array<scalar_t> vOut(localNZ.size());
    Teuchos::Array<int> pOut(localNZ.size());

std::cout << comm->getRank() << " KDDKDD buffers done " << localNZ.size() << std::endl;
    size_t cnt = 0;
    for (auto it = localNZ.begin(); it != localNZ.end(); it++, cnt++) {
      iOut[cnt] = it->first.first;
      jOut[cnt] = it->first.second;
      vOut[cnt] = it->second;
std::cout << comm->getRank() << "    KDDKDD IJ " << iOut[cnt] << " " << jOut[cnt] << std::endl;
      pOut[cnt] = procFromChunks(iOut[cnt], jOut[cnt]);
    }

std::cout << comm->getRank() << " KDDKDD buffers filled " << localNZ.size() << std::endl;
    // Free memory associated with localNZ
    LocalNZmap_t().swap(localNZ);

    // Use a Distributor to send nonzeros to new processors.
    Tpetra::Distributor plan(comm);
    size_t nrecvs = plan.createFromSends(pOut());
    Teuchos::Array<gno_t> iIn(nrecvs);
    Teuchos::Array<gno_t> jIn(nrecvs);
    Teuchos::Array<scalar_t> vIn(nrecvs);

    // TODO:  With more clever packing, could do only one round of communication
    plan.doPostsAndWaits<gno_t>(iOut(), 1, iIn());
    plan.doPostsAndWaits<gno_t>(jOut(), 1, jIn());
    plan.doPostsAndWaits<scalar_t>(vOut(), 1, vIn());

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
std::cout << "    KDDKDD procFromChunks (" << I << "," << J << "): " << p << std::endl;
    return p;
  }
};
#endif
