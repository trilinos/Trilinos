// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file mj_backwardcompat.cpp
    \brief Generate a test to backward compatibility of MJ wrt adapters
*/

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Tpetra_Map.hpp>
#include <vector>
#include <cstdlib>

/////////////////////////////////////////////
// Simple adapter with contiguous coordinates
template <typename User>
class OldSchoolVectorAdapterContig : public Zoltan2::VectorAdapter<User>
{
public:
  typedef typename Zoltan2::InputTraits<User>::gno_t gno_t;
  typedef typename Zoltan2::InputTraits<User>::scalar_t scalar_t;

  OldSchoolVectorAdapterContig(
    const size_t nids_,
    const gno_t *gids_,
    const int dim_,
    const scalar_t *coords_,
    const scalar_t *weights_ = NULL)
  : nids(nids_), gids(gids_), dim(dim_), coords(coords_), weights(weights_)
  { }

  size_t getLocalNumIDs() const { return nids; }

  void getIDsView(const gno_t *&ids) const { ids = gids; }

  int getNumWeightsPerID() const { return (weights != NULL); }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const 
  { wgt = weights; stride = 1; }

  int getNumEntriesPerID() const { return dim; }

  void getEntriesView(const scalar_t *&coo, int &stride, int idx = 0) const
  {
    coo = &(coords[idx*nids]);
    stride = 1;
  }

private:
  const size_t nids;
  const gno_t *gids;
  const int dim;
  const scalar_t *coords;
  const scalar_t *weights;
};

//////////////////////////////////////////
// Simple adapter with strided coordinates
template <typename User>
class OldSchoolVectorAdapterStrided : public Zoltan2::VectorAdapter<User>
{
public:
  typedef typename Zoltan2::InputTraits<User>::gno_t gno_t;
  typedef typename Zoltan2::InputTraits<User>::scalar_t scalar_t;

  OldSchoolVectorAdapterStrided(
    const size_t nids_,
    const gno_t *gids_,
    const int dim_,
    const scalar_t *coords_,
    const scalar_t *weights_ = NULL)
  : nids(nids_), gids(gids_), dim(dim_), coords(coords_), weights(weights_)
  { }

  size_t getLocalNumIDs() const { return nids; }

  void getIDsView(const gno_t *&ids) const { ids = gids; }

  int getNumWeightsPerID() const { return (weights != NULL); }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const 
  { wgt = weights; stride = 1; }

  int getNumEntriesPerID() const { return dim; }

  void getEntriesView(const scalar_t *&coo, int &stride, int idx = 0) const
  {
    coo = &(coords[idx]);
    stride = dim;
  }

private:
  const size_t nids;
  const gno_t *gids;
  const int dim;
  const scalar_t *coords;
  const scalar_t *weights;
};

//////////////////////////////////////////
// Same adapter as above but now we supply the Kokkos calls, not the original
// calls. This is to verify that backwards and forward compatibility are all
// working properly.
template <typename User>
class KokkosVectorAdapter : public Zoltan2::VectorAdapter<User>
{
public:
  typedef typename Zoltan2::InputTraits<User>::gno_t gno_t;
  typedef typename Zoltan2::InputTraits<User>::scalar_t scalar_t;
  typedef Tpetra::Map<>::node_type node_t;

  KokkosVectorAdapter(
    const size_t nids_,
    const gno_t *gids_,
    const int dim_,
    const scalar_t *coords_,
    const scalar_t *weights_ = NULL)
  : nids(nids_), dim(dim_)
  {
    // create kokkos_gids in default memory space
    {
      typedef Kokkos::View<gno_t *, typename node_t::device_type> view_ids_t;
      view_ids_t gids = view_ids_t(
        Kokkos::ViewAllocateWithoutInitializing("gids"), nids);
      typename view_ids_t::HostMirror host_gids =
        Kokkos::create_mirror_view(gids);
      for(size_t n = 0; n < nids; ++n) {
        host_gids(n) = gids_[n];
      }
      Kokkos::deep_copy(gids, host_gids);
      kokkos_gids = gids; // to const View
    }

    // create kokkos_weights in default memory space
    if(weights_ != NULL)
    {
      typedef Kokkos::View<scalar_t **, typename node_t::device_type> view_weights_t;
      kokkos_weights = view_weights_t(
        Kokkos::ViewAllocateWithoutInitializing("weights"), nids, 0);
      typename view_weights_t::HostMirror host_kokkos_weights =
        Kokkos::create_mirror_view(kokkos_weights);
      for(size_t n = 0; n < nids; ++n) {
        host_kokkos_weights(n,0) = weights_[n];
      }
      Kokkos::deep_copy(kokkos_weights, host_kokkos_weights);
    }

    // create kokkos_coords in default memory space
    {
      typedef Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
        typename node_t::device_type> kokkos_coords_t;
      kokkos_coords = kokkos_coords_t(
        Kokkos::ViewAllocateWithoutInitializing("coords"), nids, dim);
      typename kokkos_coords_t::HostMirror host_kokkos_coords =
        Kokkos::create_mirror_view(kokkos_coords);
      for(size_t n = 0; n < nids; ++n) {
        for(int idx = 0; idx < dim; ++idx) {
          host_kokkos_coords(n,idx) = coords_[n+idx*nids];
        }
      }
      Kokkos::deep_copy(kokkos_coords, host_kokkos_coords);
    }
  }

  size_t getLocalNumIDs() const { return nids; }

  virtual void getIDsKokkosView(Kokkos::View<const gno_t *,
    typename node_t::device_type> &ids) const {
    ids = kokkos_gids;
  }

  int getNumWeightsPerID() const { return (kokkos_weights.size() != 0); }

  virtual void getWeightsKokkosView(Kokkos::View<scalar_t **,
    typename node_t::device_type> & wgt) const {
    wgt = kokkos_weights;
  }

  int getNumEntriesPerID() const { return dim; }

  virtual void getEntriesKokkosView(Kokkos::View<scalar_t **,
    Kokkos::LayoutLeft, typename node_t::device_type> & coo) const {
    coo = kokkos_coords;
  }

private:
  const size_t nids;
  Kokkos::View<const gno_t *, typename node_t::device_type> kokkos_gids;
  const int dim;
  Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type>
    kokkos_coords;
  Kokkos::View<scalar_t **, typename node_t::device_type> kokkos_weights;
};

//////////////////////////////////////////////
int run_test_strided_versus_contig(const std::string & algorithm) {

  int nFail = 0;

  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef double scalar_t;

  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;
  typedef OldSchoolVectorAdapterStrided<myTypes> stridedAdapter_t;
  typedef OldSchoolVectorAdapterContig<myTypes> contigAdapter_t;
  typedef KokkosVectorAdapter<myTypes> kokkosAdapter_t;
  typedef Zoltan2::EvaluatePartition<stridedAdapter_t> quality_t;

  ///////////////////////////////////////////////////////////////////////
  // Create input data.

  size_t localCount = 40;
  int dim = 3; // no higher since we are testing rcb as a control which supports dim <= 3

  // Create coordinates strided
  scalar_t *cStrided = new scalar_t [dim * localCount];
  size_t cnt = 0;
  for (size_t i = 0; i < localCount; i++)
    for (int d = 0; d < dim; d++)
      cStrided[cnt++] = d*1000 + rank*100 + i;

  // Create same coords, stored contiguously
  scalar_t *cContig = new scalar_t [dim * localCount];
  cnt = 0;
  for (int d = 0; d < dim; d++)
    for (size_t i = 0; i < localCount; i++)
      cContig[cnt++] = d*1000 + rank*100 + i;

  // Create global ids for the coordinates.
  globalId_t *globalIds = new globalId_t [localCount];
  globalId_t offset = rank * localCount;
  for (size_t i=0; i < localCount; i++) globalIds[i] = offset++;
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an MJ problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("error_check_level", "debug_mode_assertions");

  params.set("algorithm", algorithm); // test runs multijagged and rcb
  params.set("num_global_parts", nprocs+1);

  ///////////////////////////////////////////////////////////////////////
  // Test one:  No weights

  // Partition using strided coords
  stridedAdapter_t *ia1 = 
                    new stridedAdapter_t(localCount,globalIds,dim,cStrided);

  Zoltan2::PartitioningProblem<stridedAdapter_t> *problem1 =
           new Zoltan2::PartitioningProblem<stridedAdapter_t>(ia1, &params);
   
  problem1->solve();

  quality_t *metricObject1 = new quality_t(ia1, &params, comm,
					   &problem1->getSolution());
  if (rank == 0){

    metricObject1->printMetrics(std::cout);

    double imb = metricObject1->getObjectCountImbalance();
    if (imb <= 1.03)  // Should get perfect balance
      std::cout << "no weights -- balance satisfied: " << imb << std::endl;
    else {
      std::cout << "no weights -- balance failure: " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }

  // Partition using contiguous coords
  contigAdapter_t *ia2 = new contigAdapter_t(localCount,globalIds,dim,cContig);

  Zoltan2::PartitioningProblem<contigAdapter_t> *problem2 =
           new Zoltan2::PartitioningProblem<contigAdapter_t>(ia2, &params);
   
  problem2->solve();

  // Partition using contiguous coords to generate kokkos adapter
  kokkosAdapter_t *ia3 = new kokkosAdapter_t(localCount,globalIds,dim,cContig);

  Zoltan2::PartitioningProblem<kokkosAdapter_t> *problem3 =
           new Zoltan2::PartitioningProblem<kokkosAdapter_t>(ia3, &params);

  problem3->solve();

  // compare strided vs contiguous vs kokkos
  size_t ndiff = 0;
  for (size_t i = 0; i < localCount; i++) {
    if((problem1->getSolution().getPartListView()[i] !=
        problem2->getSolution().getPartListView()[i]) ||
       (problem2->getSolution().getPartListView()[i] !=
        problem3->getSolution().getPartListView()[i]))
    {
      std::cout << rank << " Error: differing parts for index " << i 
                << problem1->getSolution().getPartListView()[i] << " "
                << problem2->getSolution().getPartListView()[i] << " "
                << problem3->getSolution().getPartListView()[i] << std::endl;

      ndiff++;
    }
  }
  if (ndiff > 0) nFail++;
  else if (rank == 0) std::cout << "no weights -- comparisons OK " << std::endl;

  delete metricObject1;
  delete problem1;
  delete problem2;
  delete problem3;
  delete ia1;
  delete ia2;
  delete ia3;
   
  ///////////////////////////////////////////////////////////////////////
  // Test two:  weighted
  // Create a Zoltan2 input adapter that includes weights.
  
  scalar_t *weights = new scalar_t [localCount];
  for (size_t i=0; i < localCount; i++) weights[i] = 1 + scalar_t(rank);

  // Test with strided coords
  ia1 = new stridedAdapter_t(localCount, globalIds, dim, cStrided, weights);

  problem1 = new Zoltan2::PartitioningProblem<stridedAdapter_t>(ia1, &params);

  problem1->solve();

  metricObject1 = new quality_t(ia1, &params, comm, &problem1->getSolution());

  if (rank == 0){

    metricObject1->printMetrics(std::cout);

    double imb = metricObject1->getWeightImbalance(0);
    if (imb <= 1.03)
      std::cout << "weighted -- balance satisfied " << imb << std::endl;
    else {
      std::cout << "weighted -- balance failed " << imb << std::endl;
      nFail++;
    }
    std::cout << std::endl;
  }

  // Partition using contiguous coords
  ia2 = new contigAdapter_t(localCount, globalIds, dim, cContig, weights);

  problem2 = new Zoltan2::PartitioningProblem<contigAdapter_t>(ia2, &params);
   
  problem2->solve();

  // compare strided vs contiguous
  ndiff = 0;
  for (size_t i = 0; i < localCount; i++) {
    if (problem1->getSolution().getPartListView()[i] != 
        problem2->getSolution().getPartListView()[i]) {
      std::cout << rank << " Error: differing parts for index " << i 
                << problem1->getSolution().getPartListView()[i] << " "
                << problem2->getSolution().getPartListView()[i] << std::endl;
            
      ndiff++;
    }
  }
  if (ndiff > 0) nFail++;
  else if (rank == 0) std::cout << "weighted -- comparisons OK " << std::endl;

  delete metricObject1;
  delete problem1;
  delete problem2;
  delete ia1;
  delete ia2;

  // Test with strided coords
  if (weights) delete [] weights;
  if (cStrided) delete [] cStrided;
  if (cContig) delete [] cContig;
  if (globalIds) delete [] globalIds;

  // check result

  int gnFail;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &nFail, &gnFail);
  return gnFail;
}


int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int err = 0;

  // err += run_test_strided_versus_contig("multijagged");
  err += run_test_strided_versus_contig("rcb");

  if (comm->getRank() == 0) {
    if (err == 0) std::cout << "PASS" << std::endl;
    else  std::cout << "FAIL:  " << err << " tests failed" << std::endl;
  }

  return 0;
}

