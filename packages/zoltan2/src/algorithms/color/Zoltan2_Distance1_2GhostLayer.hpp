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
                       Teuchos::RCP<femv_t> femv,
                       bool recolor);
    
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

    void constructSecondGhostLayer(std::vector<gno_t>& ownedPlusGhosts, //this argument changes
                                   const std::vector<int>& owners,
                                   ArrayView<const gno_t> adjs,
                                   ArrayView<const offset_t> offsets,
                                   RCP<const map_t> mapOwned,
                                   std::vector< gno_t>& adjs_2GL,
                                   std::vector< offset_t>& offsets_2GL) {
      //loop through owners, count how many vertices we'll send to each processor
      int vertsendcounts[comm->getSize()];
      int vertrecvcounts[comm->getSize()];
      for(int i = 0; i < comm->getSize(); i++){
        vertsendcounts[i] = 0;
      } 
      for(size_t i = 0; i < owners.size(); i++){
        if(owners[i] != comm->getRank()) vertsendcounts[owners[i]]++;
      }
      //printf("--Rank %d: sending vert counts\n",comm->getRank());
      //alltoall to send counts out
      MPI_Alltoall(vertsendcounts,1,MPI_INT,vertrecvcounts,1,MPI_INT,MPI_COMM_WORLD);
      //construct sdispls and rdispls
      //printf("--Rank %d: computing GIDs to send\n",comm->getRank());
      int degreesdispls[comm->getSize()+1];
      int degreerdispls[comm->getSize()+1];
      degreesdispls[0]=0;
      degreerdispls[0]=0;
      int recvvertcount = 0;
      int sendvertcount = 0;
      for(int i = 1; i < comm->getSize()+1; i++){
        degreesdispls[i] = degreesdispls[i-1] + vertsendcounts[i-1];
        degreerdispls[i] = degreerdispls[i-1] + vertrecvcounts[i-1];
        recvvertcount += vertrecvcounts[i-1];
        sendvertcount += vertsendcounts[i-1];
      } 
      //reorder vertex GIDs in the order they'll be recv'd (needed for converting back and forth later)
      int idx[comm->getSize()];
      //printf("--Rank %d starting send offsets:\n\t",comm->getRank());
      for(int i = 0; i < comm->getSize(); i++) {
        idx[i]=degreesdispls[i];
        //printf("%d ",idx[i]);
      }
      //printf("\n");
      int *ghostGIDs = new int[sendvertcount];
      
      int *recvGIDs = new int[recvvertcount];
      for(size_t i = offsets.size()-1; i < owners.size(); i++){
        if(owners[i] != comm->getRank()){ //vertices this proc owns are not ghosts
          ghostGIDs[idx[owners[i]]++] = ownedPlusGhosts[i];
        }
      }
      /*printf("--Rank %d ending send offsets:\n\t",comm->getRank());
      for(int i = 0; i < comm->getSize(); i++) {
        printf("%d ",idx[i]);
      }
      printf("\n");*/
      
      for(int i = 0; i < sendvertcount; i++){
        ownedPlusGhosts[i + offsets.size()-1] = ghostGIDs[i];
      }
      
      //printf("--Rank %d: sending GIDs to owning procs\n",comm->getRank());
      //Send the GIDs you need connectivity info for to their owning processors
      MPI_Alltoallv(ghostGIDs,vertsendcounts,degreesdispls,MPI_INT,
                    recvGIDs, vertrecvcounts,degreerdispls,MPI_INT,MPI_COMM_WORLD);
      //printf("--Rank %d: computing degrees for received GIDs\n",comm->getRank());
      //replace entries in recvGIDs with their degrees
      //printf("--Rank %d: received GIDs: ",comm->getRank());
      int *sendDegrees = new int[recvvertcount];
      int adj_len = 0;
      int adjsendcounts[comm->getSize()];
      for(int i = 0; i < comm->getSize(); i++){
        adjsendcounts[i] = 0;
        for(int j = degreerdispls[i]; j < degreerdispls[i+1]; j++){
          lno_t lid = mapOwned->getLocalElement(recvGIDs[j]);
          //if(lid == -1) printf("--Rank %d doesn't own vertex %d, index %d sent by rank %d\n",comm->getRank(),recvGIDs[j],j-degreerdispls[i],i);
          int degree = offsets[lid+1] - offsets[lid];
          sendDegrees[j] = degree;
          adj_len += degree;
          adjsendcounts[i]+=degree;
        }
      }
      //If it comes out burnt, I send it back! -Grati
      //(send the degrees back to the procs that need them)
      int* ghost_degrees = new int[sendvertcount];
      int* send_adjs = new int[adj_len]; //the total adjacency list, still need to determine sdispls/rdispls
      //printf("--Rank %d: sending degrees to ghost owners\n",comm->getRank());
      MPI_Alltoallv(sendDegrees,vertrecvcounts,degreerdispls,MPI_INT,
                    ghost_degrees,vertsendcounts,degreesdispls,MPI_INT,MPI_COMM_WORLD);
      delete [] sendDegrees;
      
      for(int i = 0; i < sendvertcount; i++){
        //printf("--Rank %d: ghost vertex %d has degree %d\n",comm->getRank(),ghostGIDs[i],ghost_degrees[i]);
      }
      
      //construct offests with the received vertex degrees
      int * ghost_offsets = new int[sendvertcount+1];
      ghost_offsets[0] = 0;
      //printf("--Rank %d: Ghost offsets:\n\t",comm->getRank());
      for(int i = 1; i < sendvertcount+1; i++){
        ghost_offsets[i] = ghost_offsets[i-1] + ghost_degrees[i-1];
        //printf("%d ",ghost_offsets[i]);
      }
      //printf("\n");
      //build adjacency lists to send
      //use recvGIDs to look up adjacency information
      int adjidx = 0;
      for(int i = 0; i < recvvertcount; i++){
        lno_t lid = mapOwned->getLocalElement(recvGIDs[i]);
        for(offset_t j = offsets[lid]; j<offsets[lid+1]; j++){
          send_adjs[adjidx++] = adjs[j];
        }
      }
      delete [] recvGIDs;
      //use adjsendcounts to build adjsdispls
      int adjsdispls[comm->getSize()+1];
      adjsdispls[0]=0;
      for(int i = 1; i < comm->getSize()+1; i++){
        adjsdispls[i] = adjsdispls[i-1] + adjsendcounts[i-1];
      }
      //use ghost_degrees to build up adjrevcounts and adjrdispls
      int adjrecvcounts[comm->getSize()];
      int adjrecvcount = 0;
      for(int i = 0; i < comm->getSize(); i++){
        adjrecvcounts[i] = 0;
        for(int j = degreesdispls[i]; j < degreesdispls[i+1]; j++){
          adjrecvcounts[i] += ghost_degrees[j];
          adjrecvcount+=ghost_degrees[j];
        }
      }
      delete [] ghost_degrees;
      int adjrdispls[comm->getSize()+1];
      adjrdispls[0]=0;
      for(int i = 1; i < comm->getSize()+1; i++){
        adjrdispls[i] = adjrdispls[i-1] + adjrecvcounts[i-1];
      }
      int * ghost_adjs = new int[adjrecvcount];
      
      MPI_Alltoallv(send_adjs, adjsendcounts, adjsdispls, MPI_INT,
                    ghost_adjs,adjrecvcounts, adjrdispls, MPI_INT, MPI_COMM_WORLD);
      delete [] send_adjs;
      //recvbuf should be the adjacency information
     /* printf("--Rank %d received %d adjacencies:\n\t",comm->getRank(),adjrecvcount);
      for(int i = 0; i < adjrecvcount; i++){
        printf("%d ",ghost_adjs[i]);
      }
      printf("\n");
      printf("--Rank %d: ghost graph:\n",comm->getRank());
      for(int i = 0; i < sendvertcount; i++){
        printf("--Rank %d: global vertex %d is adjacent to: ",comm->getRank(),ownedPlusGhosts[i+offsets.size()-1]);
        for(int j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
          printf("%d ",ghost_adjs[j]);
        }
        printf("\n");
      }*/
      //from this point on, I only need ghostGIDs, ghost_adjs, and ghost_offsets.
      //need to append the ghost offsets onto the local offsets, and the ghost adjacencies onto the local adjacencies.
      //This works for now, but it is not good. I need to change the types and figure out how to alter
      //MPI calls to change their sizes with the templated types (offset_t and gno_t.
      std::vector<gno_t> ghost_adjs_vec;
      std::vector<offset_t> ghost_offset_vec;
      //printf("--Rank %d: Ghost adjs vec:\n\t",comm->getRank());
      for(int i = 0; i < adjrecvcount; i++){
        adjs_2GL.push_back(ghost_adjs[i]);
        //printf("%d ",adjs_2GL.back());
      }
      //printf("\n");
      //printf("--Rank %d: Ghost offsets vec:\n\t",comm->getRank());
      for(int i = 0; i < sendvertcount+1; i++){
        offsets_2GL.push_back(ghost_offsets[i]);
        //printf("%d ", offsets_2GL.back());
      }
      //printf("\n");
      
      //adjs_2GL = Teuchos::arrayView<const gno_t>((gno_t*)ghost_adjs, adjrecvcount);
      //offsets_2GL = Teuchos::arrayView<const offset_t>((offset_t*)ghost_offsets, sendvertcount+1);
          

  
    }
    
    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution){
      //Timers I want in the future:
      // Setup, timing how long it takes to construct the local graph
      // Color, how long KokkosKernels takes
      // Communicate, how long the communication takes
      // Conflict-Resolve, how long the final resolution takes
      Teuchos::RCP<Teuchos::Time>
             timeReorder(Teuchos::TimeMonitor::getNewTimer("00 REORDER"));
      //convert from global graph to local graph
      timeReorder->start();
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      // the weights are not used at this point.
      
      ArrayView<const gno_t> adjs;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);      
      

      std::unordered_map<gno_t,lno_t> globalToLocal;
      std::vector<gno_t> ownedPlusGhosts;
      std::vector<int> owners;
      for(int i = 0; i < vtxIDs.size(); i++){
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
        owners.push_back(comm->getRank());
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
      
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                           <Tpetra::global_size_t>::invalid();
      RCP<const map_t> mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));
      std::vector<gno_t> ghosts;
      std::vector<int> ghostowners;
      for(size_t i = nVtx; i < nVtx+nGhosts; i++){
        ghosts.push_back(ownedPlusGhosts[i]);
        ghostowners.push_back(-1);
      }
      
      //get the owners of the vertices
      ArrayView<int> owningProcs = Teuchos::arrayViewFromVector(ghostowners);
      ArrayView<const gno_t> gids = Teuchos::arrayViewFromVector(ghosts);
      Tpetra::LookupStatus ls = mapOwned->getRemoteIndexList(gids, owningProcs);
      
      for(size_t i = 0; i < ghostowners.size(); i++){
        owners.push_back(ghostowners[i]);
      }
      
      //use the mapOwned to find the owners of ghosts
      std::vector< gno_t> first_layer_ghost_adjs;
      std::vector< offset_t> first_layer_ghost_offsets;
      constructSecondGhostLayer(ownedPlusGhosts,owners, adjs, offsets, mapOwned,
                                first_layer_ghost_adjs, first_layer_ghost_offsets);
      //we potentially reordered the local IDs of the ghost vertices, so we need
      //to re-insert the GIDs into the global to local ID mapping.
      globalToLocal.clear();
      for(size_t i = 0; i < ownedPlusGhosts.size(); i++){
        globalToLocal[ownedPlusGhosts[i]] = i;
      }
      
      for(int i = 0 ; i < adjs.size(); i++){
        local_adjs[i] = globalToLocal[adjs[i]];
      }
      //at this point, we have ownedPlusGhosts with 1layer ghosts' GIDs.
      //need to add 2layer ghost GIDs, and add them to the map.
      size_t n2Ghosts = 0;
      std::vector<lno_t> local_ghost_adjs;
      for(size_t i = 0; i< first_layer_ghost_adjs.size(); i++ ){
        if(globalToLocal.count(first_layer_ghost_adjs[i]) == 0){
          ownedPlusGhosts.push_back(first_layer_ghost_adjs[i]);
          globalToLocal[first_layer_ghost_adjs[i]] = vtxIDs.size() + nGhosts + n2Ghosts;
          n2Ghosts++;
        }
        local_ghost_adjs.push_back(globalToLocal[first_layer_ghost_adjs[i]]);
      }

      dummy = Teuchos::OrdinalTraits <Tpetra::global_size_t>::invalid();
      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy,
                                           Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                           0, comm));
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned,
                                                    importer, 1, true));
      
      timeReorder->stop();
      //Create random numbers seeded on global IDs, as in AlgHybridGMB.
      //This may or may not have an effect on this algorithm, but we
      //might as well see.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(size_t i = 0; i < rand.size(); i++){
        Zoltan_Srand((unsigned int) ownedPlusGhosts[i], NULL);
        rand[i] = (int) (((double) Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
      }

      

      Teuchos::ArrayView<const lno_t> local_adjs_view = Teuchos::arrayViewFromVector(local_adjs);
      Teuchos::ArrayView<const offset_t> ghost_offsets = Teuchos::arrayViewFromVector(first_layer_ghost_offsets);
      Teuchos::ArrayView<const lno_t> ghost_adjacencies = Teuchos::arrayViewFromVector(local_ghost_adjs);
      //call the coloring algorithm
      twoGhostLayer(nVtx, nVtx+nGhosts, local_adjs_view, offsets, ghost_adjacencies, ghost_offsets,
                    femv, ownedPlusGhosts, globalToLocal, rand);
      
      //copy colors to the output array
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
         colors[i] = femv->getData(0)[i];
      }
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();  
    }

    void twoGhostLayer(const size_t n_local, const size_t n_total,
                       const Teuchos::ArrayView<const lno_t>& adjs,
                       const Teuchos::ArrayView<const offset_t>& offsets,
                       const Teuchos::ArrayView<const lno_t>& ghost_adjs,
                       const Teuchos::ArrayView<const offset_t>& ghost_offsets,
                       const Teuchos::RCP<femv_t>& femv,
                       const std::vector<gno_t>& gids, 
                       const std::unordered_map<gno_t,lno_t>& globalToLocal,
                       const std::vector<int>& rand){
      
      Teuchos::RCP<Teuchos::Time>
             timeInterior(Teuchos::TimeMonitor::getNewTimer("01 INTERIOR")),
             timeBoundary(Teuchos::TimeMonitor::getNewTimer("02 BOUNDARY")),
             timeBoundaryComm(Teuchos::TimeMonitor::getNewTimer("03 BOUNDARY-COMM")),
             timeBoundaryComp(Teuchos::TimeMonitor::getNewTimer("04 BOUNDARY-COMP"));
      timeBoundary->start();
      timeInterior->start();
      Kokkos::View<offset_t*, device_type> host_offsets("Host Offset View", offsets.size());
      Kokkos::View<lno_t*, device_type> host_adjs("Host Adjacencies View", adjs.size());
      for(Teuchos_Ordinal i = 0; i < offsets.size(); i++) host_offsets(i) = offsets[i];
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) host_adjs(i) = adjs[i];

      //give the entire local graph to KokkosKernels to color
      this->colorInterior(n_local, host_adjs, host_offsets, femv,false);
      timeInterior->stop();

      //communicate the initial coloring.
      timeBoundaryComm->start();
      femv->switchActiveMultiVector();
      femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
      femv->switchActiveMultiVector();
      timeBoundaryComm->stop();
      
      //create the graph structures which allow KokkosKernels to recolor the conflicting vertices
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> dist_degrees("Owned+Ghost degree view",rand.size());
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) dist_degrees(adjs[i])++;
      for(Teuchos_Ordinal i = 0; i < ghost_adjs.size(); i++) dist_degrees(ghost_adjs[i])++;
      for(Teuchos_Ordinal i = 0; i < offsets.size()-1; i++) dist_degrees(i) = offsets[i+1] - offsets[i];
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++) dist_degrees(i+n_local) = ghost_offsets[i+1] - ghost_offsets[i];
      
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> dist_offsets("Owned+Ghost Offset view", rand.size()+1);
      dist_offsets(0) = 0;
      offset_t total_adjs = 0;
      for(size_t i = 1; i < rand.size()+1; i++){
        dist_offsets(i) = dist_degrees(i-1) + dist_offsets(i-1);
        total_adjs += dist_degrees(i-1);
      }
      Kokkos::View<lno_t*, Tpetra::Map<>::device_type> dist_adjs("Owned+Ghost adjacency view", total_adjs);
      for(size_t i = 0; i < rand.size(); i++){
        dist_degrees(i) = 0;
      }
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) dist_adjs(i) = adjs[i];
      for(Teuchos_Ordinal i = adjs.size(); i < adjs.size() + ghost_adjs.size(); i++) dist_adjs(i) = ghost_adjs[i-adjs.size()];
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++){
        for(offset_t j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
          if((size_t)ghost_adjs[j] >= n_total){
            dist_adjs(dist_offsets(ghost_adjs[j]) + dist_degrees(ghost_adjs[j])) = i + n_local;
            dist_degrees(ghost_adjs[j])++;
          }
        }
      }

      
      //we can find all the conflicts with one loop through the ghost vertices.
      timeBoundaryComp->start();
      
      //this view represents how many conflicts were found
      Kokkos::View<int[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      //keep an atomic version so that we can increment from multiple threads
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize;

      //create views for the ghost adjacencies, as they can detect all conflicts.
      Kokkos::View<offset_t*, device_type> ghost_offset_view("Ghost Offsets", ghost_offsets.size());
      Kokkos::View<lno_t*, device_type> ghost_adjs_view("Ghost Adjacencies", ghost_adjs.size());
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size(); i++) ghost_offset_view(i) = ghost_offsets[i];
      for(Teuchos_Ordinal i = 0; i < ghost_adjs.size(); i++) ghost_adjs_view(i) = ghost_adjs[i];

      //get the color view from the FEMultiVector
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);

      //create a view for the tie-breaking numbers.
      Kokkos::View<int*, device_type> rand_view("Random View", rand.size());
      for(size_t i = 0; i < rand.size(); i ++) rand_view(i) = rand[i];
      
      //detect conflicts only for ghost vertices
      Kokkos::parallel_for(ghost_offsets.size()-1, KOKKOS_LAMBDA (const int& i){
        lno_t localIdx = i + n_local;
        for(offset_t j = ghost_offset_view(i); j < ghost_offset_view(i+1); j++){
          int currColor = femv_colors(localIdx);
          int nborColor = femv_colors(ghost_adjs_view(j));
          if(currColor == nborColor ){
            if(rand_view(localIdx) > rand_view(ghost_adjs_view(j))){
              recoloringSize_atomic(0)++;
              femv_colors(localIdx) = 0;
            }
            if(rand_view(ghost_adjs_view(j)) > rand_view(localIdx)){
              recoloringSize_atomic(0)++;
              femv_colors(ghost_adjs_view(j)) = 0;
            }
          }
        }
      });
      //ensure that the parallel_for finishes before continuing
      Kokkos::fence();
      //all conflicts detected!

      //variables for statistics
      int vertsPerRound[100];
      int distributedRounds = 0;
      
      //see if recoloring is necessary.
      bool done = (recoloringSize(0) == 0);
      if(comm->getSize()==1) done = true;
      
      //recolor until no conflicts are left
      while(!done){
        if(distributedRounds < 100) vertsPerRound[distributedRounds++] = recoloringSize(0);
       
        //recolor using KokkosKernels' coloring function 
        this->colorInterior(femv_colors.size(), dist_adjs, dist_offsets, femv, true);
        recoloringSize(0) = 0;
        
        //communicate the new colors
        timeBoundaryComp->stop();
        timeBoundaryComm->start();
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();
        timeBoundaryComm->stop();
        timeBoundaryComp->start();

        //check for further conflicts
        Kokkos::parallel_for(ghost_offsets.size()-1, KOKKOS_LAMBDA (const int& i){
          lno_t localIdx = i + n_local;
          for(offset_t j = ghost_offset_view(i); j < ghost_offset_view(i+1); j++){
            int currColor = femv_colors(localIdx);
            int nborColor = femv_colors(ghost_adjs_view(j));
            if(currColor == nborColor ){
              if(rand_view(localIdx) > rand_view(ghost_adjs_view(j))){
                recoloringSize_atomic(0)++;
                femv_colors(localIdx) = 0;
              }
              if(rand_view(ghost_adjs_view(j)) > rand_view(localIdx)){
                recoloringSize_atomic(0)++;
                femv_colors(ghost_adjs_view(j)) = 0;
              }
            }
          }
        });
        
        //ensure the parallel_for finishes before continuing
        Kokkos::fence(); 
        int localDone = recoloringSize(0);
        int globalDone = 0;
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
        done = !globalDone;
      
      }//end coloring
      timeBoundary->stop();
      timeBoundaryComp->stop();
      
      if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n",distributedRounds);
      int totalVertsPerRound[100];
      for(int i= 0; i < 100; i++) totalVertsPerRound[i]=0;
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 100, vertsPerRound, totalVertsPerRound);
      for(int i = 0; i < std::min(distributedRounds, 100); i++){
        if(comm->getRank()==0) printf("recolored %d vertices in round %d\n", totalVertsPerRound[i],i);
      }
    }
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
                       Teuchos::RCP<femv_t> femv,
                       bool recolor=false) {
  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <size_t, lno_t, lno_t, execution_space, memory_space, memory_space> KernelHandle;
  KernelHandle kh;
  
  kh.set_team_work_size(-1);
  kh.set_shmem_size(16128);
  kh.set_suggested_team_size(-1);
  kh.set_suggested_vector_size(-1);
  kh.set_dynamic_scheduling(0);
  kh.set_verbose(0);
  if(recolor){
    kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
  } else {
    kh.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
  }
  Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
  Kokkos::View<int*, device_type> sv = subview(femvColors, Kokkos::ALL, 0);
  kh.get_graph_coloring_handle()->set_vertex_colors(sv);
  
  KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx,nVtx, offset_view, adjs_view);
  
  Kokkos::fence();

  numColors = kh.get_graph_coloring_handle()->get_num_colors();
  
}//end colorInterior

}//end namespace Zoltan2

#endif
