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
      for(int i = 0; i < owners.size(); i++){
        if(owners[i] != comm->getRank()) vertsendcounts[owners[i]]++;
      }
      printf("--Rank %d: sending vert counts\n",comm->getRank());
      //alltoall to send counts out
      MPI_Alltoall(vertsendcounts,1,MPI_INT,vertrecvcounts,1,MPI_INT,MPI_COMM_WORLD);
      //construct sdispls and rdispls
      printf("--Rank %d: computing GIDs to send\n",comm->getRank());
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
      for(int i = offsets.size()-1; i < owners.size(); i++){
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
      
      printf("--Rank %d: sending GIDs to owning procs\n",comm->getRank());
      //Send the GIDs you need connectivity info for to their owning processors
      MPI_Alltoallv(ghostGIDs,vertsendcounts,degreesdispls,MPI_INT,
                    recvGIDs, vertrecvcounts,degreerdispls,MPI_INT,MPI_COMM_WORLD);
      printf("--Rank %d: computing degrees for received GIDs\n",comm->getRank());
      //replace entries in recvGIDs with their degrees
      //printf("--Rank %d: received GIDs: ",comm->getRank());
      int *sendDegrees = new int[recvvertcount];
      int adj_len = 0;
      int adjsendcounts[comm->getSize()];
      for(int i = 0; i < comm->getSize(); i++){
        adjsendcounts[i] = 0;
        for(int j = degreerdispls[i]; j < degreerdispls[i+1]; j++){
          lno_t lid = mapOwned->getLocalElement(recvGIDs[j]);
          if(lid == -1) printf("--Rank %d doesn't own vertex %d, index %d sent by rank %d\n",comm->getRank(),recvGIDs[j],j-degreerdispls[i],i);
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
      printf("--Rank %d: sending degrees to ghost owners\n",comm->getRank());
      MPI_Alltoallv(sendDegrees,vertrecvcounts,degreerdispls,MPI_INT,
                    ghost_degrees,vertsendcounts,degreesdispls,MPI_INT,MPI_COMM_WORLD);
      delete [] sendDegrees;
      
      for(int i = 0; i < sendvertcount; i++){
        //printf("--Rank %d: ghost vertex %d has degree %d\n",comm->getRank(),ghostGIDs[i],ghost_degrees[i]);
      }
      
      //construct offests with the received vertex degrees
      int * ghost_offsets = new int[sendvertcount+1];
      ghost_offsets[0] = 0;
      printf("--Rank %d: Ghost offsets:\n\t",comm->getRank());
      for(int i = 1; i < sendvertcount+1; i++){
        ghost_offsets[i] = ghost_offsets[i-1] + ghost_degrees[i-1];
        //printf("%d ",ghost_offsets[i]);
      }
      printf("\n");
      //build adjacency lists to send
      //use recvGIDs to look up adjacency information
      int adjidx = 0;
      for(int i = 0; i < recvvertcount; i++){
        lno_t lid = mapOwned->getLocalElement(recvGIDs[i]);
        for(int j = offsets[lid]; j<offsets[lid+1]; j++){
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
      printf("--Rank %d: Ghost adjs vec:\n\t",comm->getRank());
      for(int i = 0; i < adjrecvcount; i++){
        adjs_2GL.push_back(ghost_adjs[i]);
        //printf("%d ",adjs_2GL.back());
      }
      printf("\n");
      printf("--Rank %d: Ghost offsets vec:\n\t",comm->getRank());
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
      std::vector<int> owners;
      for(int i = 0; i < vtxIDs.size(); i++){
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
        owners.push_back(comm->getRank());
      }
      printf("--Rank %d: OwnedPlusGhosts before ghosts:\n\t",comm->getRank());
      for(int i = 0; i< ownedPlusGhosts.size(); i++){
       // printf("%d ",ownedPlusGhosts[i]);
      }
      printf("\n");

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
      printf("--Rank %d: OwnedPlusGhosts after ghosts:\n\t",comm->getRank());
      for(int i = 0; i< ownedPlusGhosts.size(); i++){
        //printf("%d ",ownedPlusGhosts[i]);
      }
      printf("\n");
      
      RCP<const map_t> mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));
      std::vector<gno_t> ghosts;
      std::vector<gno_t> ghostowners;
      for(int i = nVtx; i < nVtx+nGhosts; i++){
        ghosts.push_back(ownedPlusGhosts[i]);
        ghostowners.push_back(-1);
      }

      //Using this for now
      //getGhostOwners(ghosts,mapOwned,ghostowners);
      //printf("--Rank %d: numGlobalElements: %d\n",comm->getRank(), mapOwned->getGlobalNumElements());
      //printf("--Rank %d: description: %s\n", comm->getRank(), mapOwned->description().c_str());
      ArrayView<int> owningProcs = Teuchos::arrayViewFromVector(ghostowners);
      ArrayView<const gno_t> gids = Teuchos::arrayViewFromVector(ghosts);
      //for( int i = 0; i < owningProcs.size(); i++){
      //  printf("--Rank %d, vertex %d is owned by proc %d\n",comm->getRank(),gids[i],owningProcs[i]);
      //}
      Tpetra::LookupStatus ls = mapOwned->getRemoteIndexList(gids, owningProcs());
      
      printf("--Rank %d owns %d vertices\n",comm->getRank(),nVtx);
      for( int i = 0; i < ghostowners.size(); i++){
        owners.push_back(ghostowners[i]);
        //printf("--Rank %d, vertex %d is owned by proc %d\n",comm->getRank(),ghosts[i],ghostowners[i]);
      }
      
      //use the mapOwned to find the owners of ghosts
      std::vector< gno_t> first_layer_ghost_adjs;
      std::vector< offset_t> first_layer_ghost_offsets;
      constructSecondGhostLayer(ownedPlusGhosts,owners, adjs, offsets, mapOwned,
                                first_layer_ghost_adjs, first_layer_ghost_offsets);
      //we potentially reordered the local IDs of the ghost vertices, so we need
      //to re-insert the GIDs into the global to local ID mapping.
      globalToLocal.clear();
      printf("--Rank %d: ownedPlusGhosts before adding second layer:\n\t",comm->getRank());
      for(int i = 0; i < ownedPlusGhosts.size(); i++){
        //printf("%d ",ownedPlusGhosts[i]);
        globalToLocal[ownedPlusGhosts[i]] = i;
      }
      printf("\n");      
      
      for(int i = 0 ; i < adjs.size(); i++){
        local_adjs[i] = globalToLocal[adjs[i]];
      }
      //at this point, we have ownedPlusGhosts with 1layer ghosts' GIDs.
      //need to add 2layer ghost GIDs, and add them to the map.
      int n2Ghosts = 0;
      std::vector<lno_t> local_ghost_adjs;
      for(int i = 0; i< first_layer_ghost_adjs.size(); i++ ){
        if(globalToLocal.count(first_layer_ghost_adjs[i]) == 0){
          ownedPlusGhosts.push_back(first_layer_ghost_adjs[i]);
          globalToLocal[first_layer_ghost_adjs[i]] = vtxIDs.size() + nGhosts + n2Ghosts;
          n2Ghosts++;
        }
        local_ghost_adjs.push_back(globalToLocal[first_layer_ghost_adjs[i]]);
      }
      printf("--Rank %d: ownedPlusGhosts:\n\t",comm->getRank());
      for(int i = 0; i < ownedPlusGhosts.size(); i++){
        //printf("%d ",ownedPlusGhosts[i]);
      }
      printf("\n");      

      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy,
                                           Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                           0, comm));
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned,
                                                    importer, 1, true));
      
      
      //Create random numbers seeded on global IDs, as in AlgHybridGMB.
      //This may or may not have an effect on this algorithm, but we
      //might as well see.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(int i = 0; i < rand.size(); i++){
        Zoltan_Srand((unsigned int) ownedPlusGhosts[i], NULL);
        rand[i] = (int) (((double) Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
      }

      //debug: print out local graph with GIDs
      /*for(int i = 0 ; i < nVtx; i ++){
        printf("--Rank %d: global vertex %d is adjacent to\n\t",comm->getRank(),ownedPlusGhosts[i]);
        for(int j = offsets[i]; j < offsets[i+1]; j++){
          printf("%d ",ownedPlusGhosts[local_adjs[j]]);
        }
        printf("\n");
      }*/
      //debug: print out ghost layer with GIDs
      /*printf("--Rank %d: ghost layer offsets:\n\t",comm->getRank());
      for(int i = 0; i < first_layer_ghost_offsets.size(); i++){
        printf("%d ",first_layer_ghost_offsets[i]);
      }
      printf("\n");
      printf("--Rank %d: ghost layer global adjs:\n\t",comm->getRank());
      for(int i = 0; i < first_layer_ghost_adjs.size(); i++){
        printf("%d ",first_layer_ghost_adjs[i]);
      }
      printf("\n");*/
      /*for(int i = 0; i < nGhosts; i++){
        printf("--Rank %d: global ghost vertex %d is adjacent to\n\t",comm->getRank(),
                                                   ownedPlusGhosts[i+nVtx]);
        for(int j = first_layer_ghost_offsets[i]; j< first_layer_ghost_offsets[i+1]; j++){
          printf("%d ",ownedPlusGhosts[local_ghost_adjs[j]]);
        }
        printf("\n");
      }*/
      

      Teuchos::ArrayView<const lno_t> local_adjs_view = Teuchos::arrayViewFromVector(local_adjs);
      Teuchos::ArrayView<const offset_t> ghost_offsets = Teuchos::arrayViewFromVector(first_layer_ghost_offsets);
      Teuchos::ArrayView<const lno_t> ghost_adjacencies = Teuchos::arrayViewFromVector(local_ghost_adjs);
      //call actual coloring algorithm
      twoGhostLayer(nVtx, nVtx+nGhosts, local_adjs_view, offsets, ghost_adjacencies, ghost_offsets,
                    femv, ownedPlusGhosts, globalToLocal, rand);
      /*printf("--Rank %d coloring: \n\t",comm->getRank());
      printf("\n\t(vertex, color) : ");
      for(int i = 0; i < femv->getData(0).size(); i++){
        printf("(%d, %d) ",ownedPlusGhosts[i],femv->getData(0)[i]);
      }
      printf("\n");*/
      //copy colors to the output array
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
         colors[i] = femv->getData(0)[i];
      }
      
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
      
      Kokkos::View<offset_t*, device_type> host_offsets("Host Offset View", offsets.size());
      Kokkos::View<lno_t*, device_type> host_adjs("Host Adjacencies View", adjs.size());
      //for(int i = 0; i < offsets.size(); i++) host_offsets(i) = offsets[i];
      //for(int i = 0; i < adjs.size(); i++) host_adjs(i) = adjs[i];
      Kokkos::parallel_for(offsets.size(),KOKKOS_LAMBDA (const int& i){
        host_offsets(i) = offsets[i];
      });
      Kokkos::parallel_for(adjs.size(), KOKKOS_LAMBDA (const int& i){
        host_adjs(i) = adjs[i];
      });

      //give the entire local graph to KokkosKernels to color
      this->colorInterior(n_local, host_adjs, host_offsets, femv);

      //communicate!
      femv->switchActiveMultiVector();
      femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
      femv->switchActiveMultiVector();

      /*printf("--Rank %d coloring (BEFORE RESOLUTION): \n\t",comm->getRank());
      printf("\n\t(vertex, color) : ");
      for(int i = 0; i < femv->getData(0).size(); i++){
        printf("(%d, %d) ",gids[i],femv->getData(0)[i]);
      }
      printf("\n");*/
      //resolve conflicts on the ghosts first, then the boundary vertices
      //(can keep a queue of the boundary vertices when we see them from the ghosts)
      
      //we can find all the conflicts with one loop through the ghost vertices.
      std::priority_queue<gno_t> resolve;
      //std::queue<lno_t> resolve_ghosts;
      //std::queue<lno_t> resolve_local;
      printf("Finding initial conflicts\n");
      Kokkos::View<lno_t*, device_type> conflicts("Conflict Array",n_total);
      Kokkos::parallel_for(n_total, KOKKOS_LAMBDA(const int& i){
        conflicts(i) = -1;
      });
      Kokkos::View<lno_t*,device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > conflicts_atomic = conflicts; 
      Kokkos::View<int[1], device_type> conflict_count("Number of conflicts");
      conflict_count(0) = 0; 
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > conflict_count_atomic = conflict_count;
      //Kokkos::View<offset_t*, device_type> ghost_offset_view("Ghost Offsets", ghost_offsets.size());
      //Kokkos::View<lno_t*, device_type> ghost_adjs_view("Ghost Adjacencies", ghost_adjs.size());
      //Kokkos::parallel_for(ghost_offsets.size(), KOKKOS_LAMBDA (const int& i){
      //  ghost_offset_view(i) = ghost_offsets[i];
      //});
      //Kokkos::parallel_for(ghost_adjs.size(), KOKKOS_LAMBDA (const int& i){
      //  ghost_adjs_view(i) = ghost_adjs[i];
      //});
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
      //Kokkos::View<int*, device_type> rand_view("Random View", rand.size());
      //Kokkos::parallel_for(rand.size(), KOKKOS_LAMBDA(const int& i){
      //  rand_view(i) = rand[i];
      //});
      Kokkos::parallel_for(ghost_offsets.size()-1, KOKKOS_LAMBDA (const int& i){
        lno_t localIdx = i + n_local;
        for( offset_t j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
          int currColor = femv_colors(localIdx);
          int nborColor = femv_colors(ghost_adjs[j]);
          if(currColor == nborColor){
            if(rand[localIdx] < rand[ghost_adjs[j]]){
              femv_colors(localIdx) = 0;
              conflict_count_atomic(0)++;
              //conflicts_atomic(conflict_count_atomic(0)++) = localIdx;
            } else if(ghost_adjs[j] < n_total){
              femv_colors(ghost_adjs[j]) = 0;
              conflict_count_atomic(0)++;
              //conflicts_atomic(conflict_count_atomic(0)++) = ghost_adjs[j];
            }
          }
        }  
      });
      /*for(int i = n_local; i < n_total; i++){
        offset_t ghost_idx = i - n_local;
        for(offset_t j = ghost_offsets[ghost_idx]; j < ghost_offsets[ghost_idx+1]; j++){
          int currColor = femv->getData(0)[i];
          int nborColor = femv->getData(0)[ghost_adjs[j]];
          if(currColor == nborColor){
            //printf("--Rank %d: Detected conflict between global %d and %d\n",comm->getRank(),gids[i],gids[ghost_adjs[j]]);
            if(gids[i] < gids[ghost_adjs[j]]){
              resolve.push(gids[i]);
              //resolve_ghosts.push(ghost_idx);
            } else { //We're resolving the neighbor, need to figure out what the neighbor is.
              if(ghost_adjs[j] < n_total){ //if the neighbor is a local or ghost vertex
                //resolve_local.push(ghost_adjs[j]);
                resolve.push(gids[ghost_adjs[j]]);     
              } else {
                //the conflict is between a double-ghost and a ghost. Just assume the double ghost will change.
                //resolve_ghosts.push(ghost_idx);
              }
            }
          }
        }
      }*/

      //printf("--Rank %d: conflicts found:\n\t",comm->getRank());
      //for(int i = 0; i< conflict_count(0); i++){
       // printf("%d ",conflicts(i));
      //}
      printf("\n");
      Teuchos::ArrayRCP<bool> forbiddenColors;
      forbiddenColors.resize(numColors);
      for(int i = 0; i < numColors; i++) forbiddenColors[i] = false;
      
      bool done = (conflict_count(0) == 0);
      printf("Starting to resolve conflicts\n");
      while(!done){
      //while(resolve.size()){
        //while(conflict_count(0) > 0) {  
          //gno_t globalVtx = resolve.top();
          //printf("--Rank %d: recoloring vertex %d\n",comm->getRank(), globalVtx);
          //lno_t currVtx = globalToLocal.at(globalVtx);
          //resolve.pop(); 
          /*lno_t currVtx = conflicts(conflict_count(0)-1);
          conflict_count(0) -= 1;
          if(currVtx >= n_local){ // this is a ghost
            for(offset_t nborIdx = ghost_offsets[currVtx-n_local]; nborIdx < ghost_offsets[currVtx+1-n_local]; nborIdx++){
              lno_t nbor = ghost_adjs[nborIdx];
              int nborColor = femv_colors(nbor);

              if(nborColor > 0){
                if(nborColor > numColors){
                  forbiddenColors.resize(nborColor);
                  numColors = nborColor;
                }
                forbiddenColors[nborColor-1] = true;
              }
            }
          } else {
            for(offset_t nborIdx = offsets[currVtx]; nborIdx < offsets[currVtx+1]; nborIdx++){
              lno_t nbor = adjs[nborIdx];
              int nborColor = femv->getData(0)[nbor];
              if(nborColor > 0){
                if(nborColor > numColors){
                  forbiddenColors.resize(nborColor);
                  numColors=nborColor;
                }
                forbiddenColors[nborColor-1] = true;
              }
            }
          }
          bool colored = false;
          for(int i = 0; i< numColors; i++){
            if(!forbiddenColors[i]){
              femv->replaceLocalValue(currVtx,0,i+1);
              colored=true;
              break;
            }
          }
          if(!colored){
            femv->replaceLocalValue(currVtx,0,numColors+1);
            numColors++;
            forbiddenColors.resize(numColors);
          }
          for(int i = 0; i < numColors; i++){
            forbiddenColors[i] = false;
          }
          
        }*///end recolor
        this->colorInterior(n_local, host_adjs, host_offsets, femv);
        conflict_count_atomic(0) = 0;
        printf("recoloring done, communicating changes\n");
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();
        printf("checking for further conflicts\n");
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        Kokkos::parallel_for(ghost_offsets.size()-1, KOKKOS_LAMBDA (const int& i){
          lno_t localIdx = i + n_local;
          for( offset_t j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
            int currColor = femv_colors(localIdx);
            int nborColor = femv_colors(ghost_adjs[j]);
            if(currColor == nborColor){
              if(rand[localIdx] < rand[ghost_adjs[j]]){
                femv_colors(localIdx) = 0;
                conflict_count_atomic(0)++;
                //conflicts_atomic(conflict_count_atomic(0)++) = localIdx;
              } else if(ghost_adjs[j] < n_total){
                femv_colors(ghost_adjs[j]) = 0;
                conflict_count_atomic(0)++;
                //conflicts_atomic(conflict_count_atomic(0)++) = ghost_adjs[j];
              }
            }
          }  
        });
          /*for(int i = n_local; i < n_total; i++){
            offset_t ghost_idx = i - n_local;
            for(offset_t j = ghost_offsets[ghost_idx]; j < ghost_offsets[ghost_idx+1]; j++){
              int currColor = femv->getData(0)[i];
              int nborColor = femv->getData(0)[ghost_adjs[j]];
              if(currColor == nborColor){
                //printf("--Rank %d: Detected conflict between global %d and %d\n",comm->getRank(),gids[i],gids[ghost_adjs[j]]);
                if(gids[i] < gids[ghost_adjs[j]]){
                  resolve.push(gids[i]);
                  //resolve_ghosts.push(ghost_idx);
                } else { //We're resolving the neighbor, need to figure out what the neighbor is.
                  if(ghost_adjs[j] < n_total){ //if the neighbor is a local or ghost vertex
                    //resolve_local.push(ghost_adjs[j]);
                    resolve.push(gids[ghost_adjs[j]]);     
                  } else {
                    //the conflict is between a double-ghost and a ghost. Just assume the double ghost will change.
                    //resolve_ghosts.push(ghost_idx);
                  }
                }
              }
            }
          }*/
          int localDone = conflict_count(0);
          int globalDone = 0;
          Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
          printf("global conflicts: %d\n",globalDone); 
          done = !globalDone;
      
      }//end coloring
      /*while(resolve_ghosts.size()){
        lno_t currVtx = resolve_ghosts.front();
        resolve_ghosts.pop();
        for(offset_t nborIdx = ghost_offsets[currVtx]; nborIdx < ghost_offsets[currVtx+1]; nborIdx++){
          lno_t nbor = ghost_adjs[nborIdx];
          int nborColor = femv->getData(0)[nbor];
          
          if(nborColor > 0){
            if(nborColor > numColors){
              forbiddenColors.resize(nborColor);
              numColors = nborColor;
            }
            forbiddenColors[nborColor-1] = true;
          }
        }
        //pick the color for currVtx based on the forbidden array
        //keep in mind currVtx is a ghost index; to get lid need to add n_local.
        bool colored = false;
        for(int i = 0; i < numColors; i++){
          if(!forbiddenColors[i]){
            femv->replaceLocalValue(currVtx+n_local,0,i+1);
            colored = true;
            break;
          }
        }
        if(!colored){
          femv->replaceLocalValue(currVtx+n_local,0,numColors+1);
          numColors++;
          forbiddenColors.resize(numColors);
        }
        for(int i = 0; i < numColors; i++){
          forbiddenColors[i] = false;
        }
      }
      
      while(resolve_local.size()){
        lno_t currVtx = resolve_local.front();
        resolve_local.pop();
        for(offset_t nborIdx = offsets[currVtx]; nborIdx< offsets[currVtx+1]; nborIdx++){
          lno_t nbor = adjs[nborIdx];
          int nborColor = femv->getData(0)[nbor];
          if(nborColor > 0){
            if(nborColor > numColors){
              forbiddenColors.resize(nborColor);
              numColors = nborColor;
            }
            forbiddenColors[nborColor-1] = true;
          }
        }
        bool colored = false;
        for(int i = 0; i < numColors; i++){
          if(!forbiddenColors[i]){
            femv->replaceLocalValue(currVtx,0,i+1);
            colored = true;
            break;
          }
        }
        if(!colored){
          femv->replaceLocalValue(currVtx,0,numColors+1);
          numColors++;
          forbiddenColors.resize(numColors);
        }
        for(int i = 0; i < numColors; i++){
          forbiddenColors[i] = false;
        }
      }*/
      
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
                       Teuchos::RCP<femv_t> femv) {
  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <size_t, lno_t, lno_t, execution_space, memory_space, memory_space> KernelHandle;
  KernelHandle kh;
  
  kh.set_team_work_size(-1);
  kh.set_shmem_size(16128);
  kh.set_suggested_team_size(-1);
  kh.set_suggested_vector_size(-1);
  kh.set_dynamic_scheduling(0);
  kh.set_verbose(0);
  kh.create_graph_coloring_handle(KokkosGraph::COLORING_VB);
  Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
  Kokkos::View<int*, device_type> sv = subview(femvColors, Kokkos::ALL, 0);
  kh.get_graph_coloring_handle()->set_vertex_colors(sv);
  
  KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx,nVtx, offset_view, adjs_view);
  
  Kokkos::fence();

  numColors = kh.get_graph_coloring_handle()->get_num_colors();
  
}//end colorInterior

}//end namespace Zoltan2

#endif
