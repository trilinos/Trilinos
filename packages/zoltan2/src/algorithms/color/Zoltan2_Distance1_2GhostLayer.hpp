#ifndef _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_
#define _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>

#include "Zoltan2_Algorithm.hpp"
#include "Zoltan2_GraphModel.hpp"
#include "Zoltan2_ColoringSolution.hpp"
#include "Zoltan2_Util.hpp"
#include "Zoltan2_TPLTraits.hpp"
#include "Zoltan2_AlltoAll.hpp"


#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"

/////////////////////////////////////////////////
//! \file Zoltan2_Distance1_2GhostLayer.hpp
//! \brief A Communication Avoidant Distance-1 Coloring Algorithm


namespace Zoltan2{

template <typename Adapter>
class AlgDistance1TwoGhostLayer : public Algorithm<Adapter> {

  public:
    
    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t,gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>; 
    using device_type = Tpetra::Map<>::device_type;
    using execution_space = Tpetra::Map<>::execution_space;
    using memory_space = Tpetra::Map<>::memory_space;
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
    

    void constructSecondGhostLayer(std::vector<gno_t>& ownedPlusGhosts, //this argument changes
                                   const std::vector<int>& owners,
                                   ArrayView<const gno_t> adjs,
                                   ArrayView<const offset_t> offsets,
                                   RCP<const map_t> mapOwned,
                                   std::vector< gno_t>& adjs_2GL,
                                   std::vector< offset_t>& offsets_2GL) {
      Teuchos::Array<int> sendcounts(comm->getSize(),0);
      Teuchos::Array<gno_t> sdispls(comm->getSize()+1,0);
      //loop through owners, count how many vertices we'll send to each processor
      for(size_t i = 0; i < owners.size(); i++){
        if(owners[i] != comm->getRank()&& owners[i] !=-1) sendcounts[owners[i]]++;
      }
      //construct sdispls (for building sendbuf), and sum the total sendcount
      gno_t sendcount = 0;
      for(int i = 1; i < comm->getSize()+1; i++){
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
        sendcount += sendcounts[i-1];
      }
      std::vector<gno_t> idx(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++){
        idx[i]=sdispls[i];
      }
      //construct sendbuf to send GIDs to owning processes
      std::vector<gno_t> sendbuf(sendcount,0);
      for(size_t i = offsets.size()-1; i < owners.size(); i++){
        if(owners[i] != comm->getRank() && owners[i] != -1){
          sendbuf[idx[owners[i]]++] = ownedPlusGhosts[i];
        }
      }
      //remap the GIDs so we receive the adjacencies in the same order as the current processes LIDs
      for(gno_t i = 0; i < sendcount; i++){
        ownedPlusGhosts[i+offsets.size()-1] = sendbuf[i];
      }
      
      //communicate GIDs to owners
      //Teuchos::ArrayView<int> sendcounts_view = Teuchos::arrayViewFromVector(sendcounts);
      //Teuchos::ArrayView<gno_t> sendbuf_view = Teuchos::arrayViewFromVector(sendbuf);
      Teuchos::ArrayRCP<gno_t>  recvbuf;
      Teuchos::Array<int> recvcounts(comm->getSize(),0);
      //Teuchos::ArrayView<int> recvcounts_view = Teuchos::arrayViewFromVector(recvcounts);
      Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf, sendcounts, recvbuf, recvcounts);
      
      //replace entries in recvGIDs with their degrees
      gno_t recvcounttotal = 0;
      Teuchos::Array<int> rdispls(comm->getSize()+1,0);
      for(size_t i = 1; i<recvcounts.size()+1; i++){
        rdispls[i] = rdispls[i-1] + recvcounts[i-1];
        recvcounttotal += recvcounts[i-1];
      }
      //send back the degrees to the requesting processes,
      //build the adjacency counts
      Teuchos::Array<offset_t> sendDegrees(recvcounttotal,0);
      gno_t adj_len = 0;
      Teuchos::Array<int> adjsendcounts(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++){
        adjsendcounts[i] = 0;
        for(int j = rdispls[i]; j < rdispls[i+1]; j++){
          lno_t lid = mapOwned->getLocalElement(recvbuf[j]);
          offset_t degree = offsets[lid+1] - offsets[lid];
          sendDegrees[j] = degree;
          adj_len +=degree;
          adjsendcounts[i] += degree; 
        }
      }
      //communicate the degrees back to the requesting processes
      //Teuchos::ArrayView<offset_t> sendDegrees_view = Teuchos::arrayViewFromVector(sendDegrees);
      Teuchos::ArrayRCP<offset_t> recvDegrees;
      Teuchos::Array<int> recvDegreesCount(comm->getSize(),0);
      //Teuchos::ArrayView<int> recvDegreesCount_view = Teuchos::arrayViewFromVector(recvDegreesCount);
      Zoltan2::AlltoAllv<offset_t>(*comm, *env, sendDegrees, recvcounts, recvDegrees, recvDegreesCount);
     
      //calculate number of rounds of AlltoAllv's that are necessary on this process
      int rounds = 1;
      for(int i = 0; i < comm->getSize(); i++){
        if(adjsendcounts[i]*sizeof(gno_t)/ INT_MAX > rounds){
	  rounds = (adjsendcounts[i]*sizeof(gno_t)/INT_MAX)+1;
	}
      }

      //see what the max number of rounds will be globally
      int max_rounds = 0;
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_MAX, 1, &rounds, &max_rounds);
      
      Teuchos::Array<Teuchos::Array<gno_t> > per_proc_round_adj_sums;
      Teuchos::Array<Teuchos::Array<gno_t> > per_proc_round_vtx_sums;
      per_proc_round_adj_sums.resize(max_rounds+1,Teuchos::Array<gno_t>(comm->getSize(),0));
      per_proc_round_vtx_sums.resize(max_rounds+1,Teuchos::Array<gno_t>(comm->getSize(),0));
      
      for(int proc_to_send = 0; proc_to_send < comm->getSize(); proc_to_send++){
        int curr_round = 0;
	for(int j = sdispls[proc_to_send]; j < sdispls[proc_to_send+1]; j++){
	  if((per_proc_round_adj_sums[curr_round][proc_to_send] + recvDegrees[j])*sizeof(gno_t) > INT_MAX){
	    curr_round++;
	  }
	  per_proc_round_adj_sums[curr_round][proc_to_send] += recvDegrees[j];
	  per_proc_round_vtx_sums[curr_round][proc_to_send]++;
	}
      }

      Teuchos::Array<Teuchos::Array<Teuchos::Array<gno_t> > > recv_GID_per_proc_per_round;
      recv_GID_per_proc_per_round.resize(max_rounds+1);
      for(int i = 0; i < max_rounds; i++){
	recv_GID_per_proc_per_round[i].resize(comm->getSize());
        for(int j = 0; j < comm->getSize(); j++){
          recv_GID_per_proc_per_round[i][j].resize(sendcounts[j],0);
        }
      }
      
      std::cout<<comm->getRank()<<" finished constructing recv_GID_per_proc_per_round\n";

      for(int i = 0; i < comm->getSize(); i++){
        int curr_round = 0;
	int curr_idx = 0;
	for(int j = sdispls[i]; j < sdispls[i+1]; j++){
	  if(curr_idx > per_proc_round_vtx_sums[curr_round][i]){
	    curr_round++;
	    curr_idx = 0;
	  }
	  if(comm->getRank()==0) std::cout<<comm->getRank()<<" curr_round("<<curr_round<<") >= max_rounds("<<max_rounds<<"\n";
	  recv_GID_per_proc_per_round[curr_round][i][curr_idx++] = j;
	}
      }
      
      std::cout<<comm->getRank()<<" finished assigning recv_GID_per_proc_per_round\n";

      Teuchos::Array<gno_t> final_gid_arr(sendcount,0);
      Teuchos::Array<offset_t> final_degree_arr(sendcount,0);
      gno_t reorder_idx = 0;
      for(int i = 0; i < max_rounds; i++){
        for(int j = 0; j < comm->getSize(); j++){
	  for(int k = 0; k < per_proc_round_vtx_sums[i][j]; k++){
	    final_gid_arr[reorder_idx] = sendbuf[recv_GID_per_proc_per_round[i][j][k]];
	    final_degree_arr[reorder_idx++] = recvDegrees[recv_GID_per_proc_per_round[i][j][k]];
	  }
	}
      }
     
      std::cout<<comm->getRank()<<" finished reorganizing GIDs\n";

      //remap the GIDs so we receive the adjacencies in the same order as the current process's LIDs
      for(gno_t i = 0; i < sendcount; i++){
        ownedPlusGhosts[i+offsets.size()-1] = final_gid_arr[i];
      }

      //construct offsets with the received vertex degrees
      Teuchos::Array<offset_t> ghost_offsets(sendcount+1,0);
      Teuchos::Array<lno_t> send_adjs(adj_len,0);
      for(int i = 1; i < sendcount+1; i++){
        ghost_offsets[i] = ghost_offsets[i-1] + recvDegrees[i-1];
      }
      //build adjacency lists to send back
      //use recvGIDs to look up adjacency information
      offset_t adjidx = 0;
      for (int i = 0; i < recvcounttotal; i++){
        lno_t lid = mapOwned->getLocalElement(recvbuf[i]);
        for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
          send_adjs[adjidx++] = adjs[j];
        }
      }
      //count how many adjacencies we will receive
      offset_t recvadjscount = 0;
      for(int i = 0; i < recvDegrees.size(); i++){
        recvadjscount+= recvDegrees[i];
      }

      Teuchos::Array<gno_t> curr_idx_per_proc(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++) curr_idx_per_proc[i] = rdispls[i];
      for(int round = 0; round < max_rounds; round++){
        Teuchos::Array<gno_t> send_adj;
	Teuchos::Array<int> send_adj_counts(comm->getSize(),0);
	for(int curr_proc=0; curr_proc < comm->getSize(); curr_proc++){
	  gno_t curr_adj_sum = 0;
	  while( curr_idx_per_proc[curr_proc] < rdispls[curr_proc+1]){
	    lno_t lid = mapOwned->getLocalElement(recvbuf[curr_idx_per_proc[curr_proc]++]);
	    if((curr_adj_sum + (offsets[lid+1] - offsets[lid]))*sizeof(gno_t) >= INT_MAX){
	      break;
	    }
	    curr_adj_sum += (offsets[lid+1] - offsets[lid]);
	    for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
	      send_adj.push_back(adjs[j]);
	    }
	  }
	  send_adj_counts[curr_proc] = curr_adj_sum;
	}

        //communicate adjacencies back to requesting processes
        //Teuchos::ArrayView<lno_t> send_adjs_view = Teuchos::arrayViewFromVector(send_adjs);
        //Teuchos::ArrayView<int> adjsendcounts_view = Teuchos::arrayViewFromVector(adjsendcounts);
        Teuchos::ArrayRCP<gno_t> ghost_adjs;
	gno_t recv_adjs_count = 0;
	for(int i = 0; i < comm->getSize(); i++){
	  recv_adjs_count += per_proc_round_adj_sums[round][i];
	}
        Teuchos::Array<int> adjrecvcounts(comm->getSize(),0);
        //Teuchos::ArrayView<int> adjsrecvcounts_view = Teuchos::arrayViewFromVector(adjrecvcounts);
        Zoltan2::AlltoAllv<gno_t>(*comm, *env, send_adj, send_adj_counts , ghost_adjs, adjrecvcounts);

        //build the adjacencies and offsets for the second ghost layer
        for(offset_t i = 0; i < ghost_adjs.size(); i ++){
          adjs_2GL.push_back(ghost_adjs[i]);
        }
      }
      for(int i = 0; i < sendcount+1; i++){
        offsets_2GL.push_back(ghost_offsets[i]);
      }
    }
    
    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution){
      //convert from global graph to local graph
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
      using import_t = Tpetra::Import<lno_t, gno_t>;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned,
                                                    importer, 1, true));
      //Create random numbers seeded on global IDs, as in AlgHybridGMB.
      //This may or may not have an effect on this algorithm, but we
      //might as well see.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(size_t i = 0; i < rand.size(); i++){
        std::srand(ownedPlusGhosts[i]);
        rand[i] = std::rand();
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
      for(Teuchos_Ordinal i = 0; i < offsets.size(); i++) host_offsets(i) = offsets[i];
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) host_adjs(i) = adjs[i];

      //give the entire local graph to KokkosKernels to color
      this->colorInterior(n_local, host_adjs, host_offsets, femv);

      //communicate the initial coloring.
      femv->switchActiveMultiVector();
      femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
      femv->switchActiveMultiVector();
      
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
      
      //this view represents how many conflicts were found
      Kokkos::View<gno_t[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      //keep an atomic version so that we can increment from multiple threads
      Kokkos::View<gno_t[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize;

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
      
      Kokkos::View<gno_t*, device_type> gid_view("GIDs", gids.size());
      for(size_t i = 0; i < gids.size(); i++) gid_view(i) = gids[i];
 
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
            }else if(rand_view(ghost_adjs_view(j)) > rand_view(localIdx)){
              recoloringSize_atomic(0)++;
              femv_colors(ghost_adjs_view(j)) = 0;
            } else {
              if (gid_view(localIdx) >= gid_view(ghost_adjs_view(j))){
                femv_colors(localIdx) = 0;
                recoloringSize_atomic(0)++;
              } else {
                femv_colors(ghost_adjs_view(j)) = 0;
                recoloringSize_atomic(0)++;
              }
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
      gno_t totalConflicts = 0;
      gno_t localConflicts = recoloringSize(0);
      Teuchos::reduceAll<int,gno_t>(*comm, Teuchos::REDUCE_SUM, 1, &localConflicts, &totalConflicts);
      bool done = !totalConflicts;
      if(comm->getSize()==1) done = true;
      
      //recolor until no conflicts are left
      while(!done){
        if(distributedRounds < 100) vertsPerRound[distributedRounds++] = recoloringSize(0);
       
        //recolor using KokkosKernels' coloring function 
        this->colorInterior(femv_colors.size(), dist_adjs, dist_offsets, femv);
        recoloringSize(0) = 0;
        
        //communicate the new colors
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();

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
              }else if(rand_view(ghost_adjs_view(j)) > rand_view(localIdx)){
                recoloringSize_atomic(0)++;
                femv_colors(ghost_adjs_view(j)) = 0;
              } else {
                if (gid_view(localIdx) >= gid_view(ghost_adjs_view(j))){
                  femv_colors(localIdx) = 0;
                  recoloringSize_atomic(0)++;
                } else {
                  femv_colors(ghost_adjs_view(j)) = 0;
                  recoloringSize_atomic(0)++;
                }
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
      
      /*RESULT REPORTING INSTRUMENTATION
      if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n",distributedRounds);
      int totalVertsPerRound[100];
      for(int i= 0; i < 100; i++) totalVertsPerRound[i]=0;
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 100, vertsPerRound, totalVertsPerRound);
      for(int i = 0; i < std::min(distributedRounds, 100); i++){
        if(comm->getRank()==0) printf("recolored %d vertices in round %d\n", totalVertsPerRound[i],i);
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
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
      <size_t, lno_t, lno_t, execution_space, memory_space, memory_space>;
  KernelHandle kh;
  
  kh.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
  
  Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
  Kokkos::View<int*, device_type> sv = subview(femvColors, Kokkos::ALL, 0);
  kh.get_graph_coloring_handle()->set_vertex_colors(sv);
  
  KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx,nVtx, offset_view, adjs_view);
  
  numColors = kh.get_graph_coloring_handle()->get_num_colors();
  
}//end colorInterior

}//end namespace Zoltan2

#endif
