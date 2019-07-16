#ifndef _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_
#define _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include <KokkosKernels_Handle.hpp>
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"

#include "zz_rand.h"
/////////////////////////////////////////////////
//! \file Zoltan2_Distance1_2GhostLayer.hpp
//! \brief A Communication Avoidant Distance-1 Coloring Algorithm


namespace Zoltan2{

template <typename Adapter>
class AlgDistance1TwoGhostLayer : public Algorithm<Adapter> {

  public:
    
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;
    typedef typename Adapter::base_adapter_t base_adapter_t;
    typedef Tpetra::Map<lno_t,gno_t> map_t;
    typedef int femv_scalar_t;
    typedef Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t> femv_t; 
    typedef Tpetra::Map<>::device_type device_type;
    typedef Tpetra::Map<>::execution_space execution_space;
    typedef Tpetra::Map<>::memory_space memory_space;
  private:
    
    void buildModel(modelFlag_t &flags);
    void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,device_type> adjs_view,
                       Kokkos::View<offset_t*, device_type> offset_view,
                       Teuchos::ArrayRCP<int> colors,
                       Teuchos::RCP<femv_t> femv);
    
    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    int numColors;

  public:
    AlgDistance1TwoGhostLayer(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_){

      numColors = 4;
      modelFlag_t flags;
      flags.reset();
      buildModel(flags);
    }
    
    void getGhostOwners(const std::vector<gno_t>& ghosts,
                        RCP<const map_t> mapOwned,
                        std::vector<int>& owners){
      
      int ghostsendcounts[comm->getSize()];
      int ghostrecvcounts[comm->getSize()];
      for(int i = 0; i < comm->getSize(); i++){
        ghostsendcounts[i] = ghosts.size();
      }

      //send the length of the ghosts on each processor to each processor
      MPI_Alltoall(ghostsendcounts,1,MPI_INT,ghostrecvcounts,1,MPI_INT, MPI_COMM_WORLD);

      int ghostsenddispls[comm->getSize()+1];
      int ghostrecvdispls[comm->getSize()+1];
      ghostsenddispls[0] = 0;
      ghostrecvdispls[0] = 0;
      
      int recvcount = 0;
      for(int i = 1; i < comm->getSize()+1; i ++){
        ghostsenddispls[i] = ghostsenddispls[i-1] + ghostsendcounts[i-1];
        ghostrecvdispls[i] = ghostrecvdispls[i-1] + ghostrecvcounts[i-1];
        recvcount += ghostrecvcounts[i-1];
      }

      int* ghostsendbuf = new int[comm->getSize() * ghosts.size()];
      int* ghostrecvbuf = new int[recvcount];
      
      for(int i = 0; i< comm->getSize(); i++){
        for(int j = 0; j < ghosts.size(); j++){
          ghostsendbuf[i*ghosts.size()+j] = ghosts[j];
        }
      }
      
      MPI_Alltoallv(ghostsendbuf,ghostsendcounts,ghostsenddispls,MPI_INT,
                    ghostrecvbuf,ghostrecvcounts,ghostrecvdispls,MPI_INT, MPI_COMM_WORLD);
    
      lno_t notFound = Teuchos::OrdinalTraits<lno_t>::invalid();
      //At this point, ghostrecvbuf should contain all GIDs that need owners.
      for(int i = 0; i <recvcount; i++){
        if(mapOwned->getLocalElement(ghostrecvbuf[i]) != notFound){
          ghostrecvbuf[i] = comm->getRank();
        } else {
          ghostrecvbuf[i] = -1;
        }
      }
      
      MPI_Alltoallv(ghostrecvbuf,ghostrecvcounts,ghostrecvdispls,MPI_INT,
                    ghostsendbuf,ghostsendcounts,ghostsenddispls,MPI_INT, MPI_COMM_WORLD);

      //look through the ghostsenddispls array to find the owners of your ghosted nodes.
      for(int i = 0; i < comm->getSize(); i++){
        for(int j = 0; j < ghosts.size(); j++){
          if(ghostsendbuf[ghostsenddispls[i]+j] != -1 ){
            owners[j] = ghostsendbuf[ghostsenddispls[i]+j];
          }
        }
      }
      delete [] ghostsendbuf;
      delete [] ghostrecvbuf;
    }

    void constructSecondGhostLayer(const std::vector<gno_t>& ownedPlusGhosts,
                                   const std::vector<int>& owners,
                                   ArrayView<const gno_t> adjs,
                                   ArrayView<const offset_t> offsets,
                                   RCP<const map_t> mapOwned,
                                   ArrayView<const gno_t>& adjs_2GL,
                                   ArrayView<const offset_t>& offsets_2GL) {
      //loop through owners, count how many vertices we'll send to each processor
      //alltoall to send counts out
      //alltoallv send vertex GIDs (needed for converting back and forth later)
      //construct sdispls and rdispls
      //send degrees for each vertex
      //construct offests with the received vertex degrees
      //build adjacency lists to send
      //recvbuf should be the adjacency information
    }
    
    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution){
      //Timers I want in the future:
      // Setup, timing how long it takes to construct the local graph
      // Color, how long KokkosKernels takes
      // Communicate, how long the communication takes
      // Conflict-Resolve, how long the final resolution takes
      
      //convert from global graph to local graph
      
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      //don't really need the weights
      
      ArrayView<const gno_t> adjs;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);      
      
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                           <Tpetra::global_size_t>::invalid();
      
      


      std::unordered_map<gno_t,lno_t> globalToLocal;
      std::vector<gno_t> ownedPlusGhosts;
      for(int i = 0; i < vtxIDs.size(); i++){
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
      }
      int nGhosts = 0;
      std::vector<lno_t> local_adjs;
      for(int i = 0; i < adjs.size(); i++){
        if(globalToLocal.count(adjs[i])==0){
          ownedPlusGhosts.push_back(adjs[i]);
          globalToLocal[adjs[i]] = vtxIDs.size()+nGhosts;
          nGhosts++;
        }
        local_adjs.push_back(globalToLocal[adjs[i]]);
      }
      
      RCP<const map_t> mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));
      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy,
                                           Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                           0, comm));
      std::vector<gno_t> ghosts;
      std::vector<gno_t> owners;
      for(int i = nVtx; i < nVtx+nGhosts; i++){
        ghosts.push_back(ownedPlusGhosts[i]);
        owners.push_back(-1);
      }

      //Using this for now
      getGhostOwners(ghosts,mapOwned,owners);
      //printf("--Rank %d: numGlobalElements: %d\n",comm->getRank(), mapOwned->getGlobalNumElements());
      //printf("--Rank %d: description: %s\n", comm->getRank(), mapOwned->description().c_str());
      //ArrayView<int> owningProcs = Teuchos::arrayViewFromVector(ghosts);
      //ArrayView<const gno_t> gids = Teuchos::arrayViewFromVector(ghosts);
      //for( int i = 0; i < owningProcs.size(); i++){
      //  printf("--Rank %d, vertex %d is owned by proc %d\n",comm->getRank(),gids[i],owningProcs[i]);
      //}
      //Tpetra::LookupStatus ls = mapOwned->getRemoteIndexList(gids, owningProcs());
      
      printf("--Rank %d owns %d vertices\n",comm->getRank(),nVtx);
      for( int i = 0; i < owners.size(); i++){
        printf("--Rank %d, vertex %d is owned by proc %d\n",comm->getRank(),ghosts[i],owners[i]);
      }
      
      //use the mapOwned to find the owners of ghosts
      ArrayView<const gno_t> actual_adjs;
      ArrayView<const offset_t> actual_offsets;
      //constructSecondGhostLayer(ownedPlusGhosts, adjs, offsets, mapOwned, actual_adjs, actual_offsets);
      
      /*typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned,
                                                    importer, 1, true));
      
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
         colors[i] = 0;
      }
      
      //Create random numbers seeded on global IDs, as in AlgHybridGMB.
      //This may or may not have an effect on this algorithm, but we
      //might as well see.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(int i = 0; i < rand.size(); i++){
        Zoltan_Srand((unsigned int) ownedPlusGhosts[i], NULL);
        rand[i] = (int) (((double) Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
      }*/

      //call actual coloring algorithm

      //copy colors to the output array

    }

    void twoGhostLayer(const size_t n_local, const size_t n_total,
                       Teuchos::ArrayView<const lno_t> adjs,
                       Teuchos::ArrayView<const offset_t> offsets,
                       Teuchos::RCP<femv_t> femv,
                       std::vector<gno_t> gids, std::vector<int> rand){}
}; //end class

template <typename Adapter>
void AlgDistance1TwoGhostLayer<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
  
  this->env->debug(DETAILED_STATUS, "   building graph model");
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
  this->env->debug(DETAILED_STATUS, "   graph model built");
}

template<typename Adapter>
void AlgDistance1TwoGhostLayer<Adapter>::colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,device_type> adjs_view,
                       Kokkos::View<offset_t*, device_type> offset_view,
                       Teuchos::ArrayRCP<int> colors,
                       Teuchos::RCP<femv_t> femv) {
  typedef Kokkoskernels::Experimental::KokkosKernelsHandle
      <size_t, lno_t, lno_t, execution_space, memory_space, memory_space> KernelHandle;
  KernelHandle kh;
  
  kh.set_team_work_size(-1);
  kh.set_shmem_size(16128);
  kh.set_suggested_team_size(-1);
  kh.set_suggested_vector_size(-1);
  kh.set_dynamic_scheduling(0);
  kh.set_verbose(0);
  kh.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
  Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
  Kokkos::View<int*, device_type> sv = subview(femvColors, Kokkos::ALL, 0);
  kh.get_graph_coloring_handle()->set_vertex_colors(sv);
  
  KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx,nVtx, offset_view, adjs_view);
  
  Kokkos::fence();

  numColors = kh.get_graph_coloring_handle()->get_num_colors();
  
}//end colorInterior

}//end namespace Zoltan2

#endif
