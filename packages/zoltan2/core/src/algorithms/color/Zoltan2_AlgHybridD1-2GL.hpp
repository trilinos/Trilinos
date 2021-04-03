#ifndef _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_
#define _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <sys/time.h>

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
    double timer(){
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
    }
  private:
    
    void buildModel(modelFlag_t &flags);
    template <class ExecutionSpace, typename TempMemorySpace,
	      typename MemorySpace>
    void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,Kokkos::Device<ExecutionSpace, MemorySpace> > adjs_view,
                       Kokkos::View<offset_t*,Kokkos::Device<ExecutionSpace, MemorySpace> > offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){

      //setup types to be used by KokkosKernels
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, ExecutionSpace, MemorySpace, MemorySpace>;
      using lno_row_view_t = Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;
      using lno_nnz_view_t = Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;
      KernelHandle kh;
  
      //pick which KokkosKernels algorithm to use, based on max graph degree
      if(recolor){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
      }
      //vertex_list_size indicates whether we have provided a list of vertices to recolor.
      //Only makes a difference if the algorithm to be used is VBBIT
      if(vertex_list_size != 0){
        kh.get_graph_coloring_handle()->set_vertex_list(vertex_list,vertex_list_size);
      }
  
      kh.set_verbose(verbose);
  
      //set the initial coloring to be the colors from femv.
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<MemorySpace>();
      Kokkos::View<int*, Kokkos::Device<ExecutionSpace, MemorySpace>> sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      kh.get_graph_coloring_handle()->set_tictoc(verbose);
  
      KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t> (&kh, nVtx,nVtx, offset_view, adjs_view);

      numColors = kh.get_graph_coloring_handle()->get_num_colors();

      if(verbose){
        std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
      }
    }
    
    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    bool verbose;
    int numColors;
    std::unordered_map<lno_t, std::vector<int>> procs_to_send;

  public:
    AlgDistance1TwoGhostLayer(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_){
      verbose = pl->get<bool>("verbose",false);
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
      std::vector<int> sendcounts(comm->getSize(),0);
      std::vector<gno_t> sdispls(comm->getSize()+1,0);
      //loop through owners, count how many vertices we'll send to each processor
      if(verbose) std::cout<<comm->getRank()<<": building sendcounts\n";
      for(size_t i = 0; i < owners.size(); i++){
        if(owners[i] != comm->getRank()&& owners[i] !=-1) sendcounts[owners[i]]++;
      }
      //construct sdispls (for building sendbuf), and sum the total sendcount
      if(verbose) std::cout<<comm->getRank()<<": building sdispls\n";
      gno_t sendcount = 0;
      for(int i = 1; i < comm->getSize()+1; i++){
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
        sendcount += sendcounts[i-1];
      }

      if(verbose) std::cout<<comm->getRank()<<": building idx\n";
      std::vector<gno_t> idx(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++){
        idx[i]=sdispls[i];
      }
      //construct sendbuf to send GIDs to owning processes
      if(verbose) std::cout<<comm->getRank()<<": building sendbuf\n";

      std::vector<gno_t> sendbuf(sendcount,0);
      for(size_t i = offsets.size()-1; i < owners.size(); i++){
        if(owners[i] != comm->getRank() && owners[i] != -1){
          sendbuf[idx[owners[i]]++] = ownedPlusGhosts[i];
        }
      }
      //remap the GIDs so we receive the adjacencies in the same order as the current processes LIDs
      if(verbose) std::cout<<comm->getRank()<<": updating ownedPlusGhosts\n";
      for(gno_t i = 0; i < sendcount; i++){
        ownedPlusGhosts[i+offsets.size()-1] = sendbuf[i];
      }
      
      //communicate GIDs to owners
      if(verbose) std::cout<<comm->getRank()<<": requesting GIDs from owners\n";
      Teuchos::ArrayView<int> sendcounts_view = Teuchos::arrayViewFromVector(sendcounts);
      Teuchos::ArrayView<gno_t> sendbuf_view = Teuchos::arrayViewFromVector(sendbuf);
      Teuchos::ArrayRCP<gno_t>  recvbuf;
      std::vector<int> recvcounts(comm->getSize(),0);
      Teuchos::ArrayView<int> recvcounts_view = Teuchos::arrayViewFromVector(recvcounts);
      Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcounts_view, recvbuf, recvcounts_view);
      
      if(verbose) std::cout<<comm->getRank()<<": done communicating\n";
      //replace entries in recvGIDs with their degrees

      if(verbose) std::cout<<comm->getRank()<<": building rdispls\n";
      gno_t recvcounttotal = 0;
      std::vector<int> rdispls(comm->getSize()+1,0);
      for(size_t i = 1; i<recvcounts.size()+1; i++){
        rdispls[i] = rdispls[i-1] + recvcounts[i-1];
        recvcounttotal += recvcounts[i-1];
      }
      //send back the degrees to the requesting processes,
      //build the adjacency counts
      std::vector<offset_t> sendDegrees(recvcounttotal,0);
      gno_t adj_len = 0;
      std::vector<int> adjsendcounts(comm->getSize(),0);
      if(verbose) std::cout<<comm->getRank()<<": building adjacency counts\n";
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
      if(verbose) std::cout<<comm->getRank()<<": sending degrees back to requestors\n";
      Teuchos::ArrayView<offset_t> sendDegrees_view = Teuchos::arrayViewFromVector(sendDegrees);
      Teuchos::ArrayRCP<offset_t> recvDegrees;
      std::vector<int> recvDegreesCount(comm->getSize(),0);
      Teuchos::ArrayView<int> recvDegreesCount_view = Teuchos::arrayViewFromVector(recvDegreesCount);
      Zoltan2::AlltoAllv<offset_t>(*comm, *env, sendDegrees_view, recvcounts_view, recvDegrees, recvDegreesCount_view);
     
      //calculate number of rounds of AlltoAllv's that are necessary on this process

      if(verbose) std::cout<<comm->getRank()<<": determining number of rounds necessary\n";
      int rounds = 1;
      for(int i = 0; i < comm->getSize(); i++){
        if(adjsendcounts[i]*sizeof(gno_t)/ INT_MAX > (size_t)rounds){
          rounds = (adjsendcounts[i]*sizeof(gno_t)/INT_MAX)+1;
        }
      }

      //see what the max number of rounds will be globally
      int max_rounds = 0;
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_MAX, 1, &rounds, &max_rounds);
      
      if(verbose) std::cout<<comm->getRank()<<": building per_proc sums\n";

      std::vector<std::vector<uint64_t>> per_proc_round_adj_sums(max_rounds+1,std::vector<uint64_t>(comm->getSize(),0));
      std::vector<std::vector<uint64_t>> per_proc_round_vtx_sums(max_rounds+1,std::vector<uint64_t>(comm->getSize(),0));

      if(verbose) std::cout<<comm->getRank()<<": filling per_proc sums\n";
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
      
      if(verbose) std::cout<<comm->getRank()<<": building recv GID schedule\n";
      std::vector<std::vector<std::vector<gno_t>>> recv_GID_per_proc_per_round(
		      max_rounds+1,std::vector<std::vector<gno_t>>(
			      comm->getSize(),std::vector<gno_t>(0)));
      for(int i = 0; i < max_rounds; i++){
        for(int j = 0; j < comm->getSize(); j++){
	  recv_GID_per_proc_per_round[i][j] = std::vector<gno_t>(sendcounts[j],0);
	}
      }

      if(verbose) std::cout<<comm->getRank()<<": filling out recv GID schedule\n";
      for(int i = 0; i < comm->getSize(); i++){
        int curr_round = 0;
        size_t curr_idx = 0;
        for(int j = sdispls[i]; j < sdispls[i+1]; j++){
          if(curr_idx > per_proc_round_vtx_sums[curr_round][i]){
            curr_round++;
            curr_idx = 0;
          }
          recv_GID_per_proc_per_round[curr_round][i][curr_idx++] = j;
        }
      }
      
      if(verbose) std::cout<<comm->getRank()<<": reordering gids and degrees in the order they'll be received\n";
      std::vector<gno_t> final_gid_vec(sendcount, 0);
      std::vector<offset_t> final_degree_vec(sendcount,0);
      gno_t reorder_idx = 0;
      for(int i = 0; i < max_rounds; i++){
        for(int j = 0; j < comm->getSize(); j++){
          for(size_t k = 0; k < per_proc_round_vtx_sums[i][j]; k++){
            final_gid_vec[reorder_idx] = sendbuf[recv_GID_per_proc_per_round[i][j][k]];
            final_degree_vec[reorder_idx++] = recvDegrees[recv_GID_per_proc_per_round[i][j][k]];
          }
        }
      }
      bool reorganized = false;
      for(int i = 0; i < sendcount; i++){
        if(final_gid_vec[i] != sendbuf[i]) reorganized = true;
      }
      if(verbose){
        if(!reorganized && (max_rounds > 1)) std::cout<<comm->getRank()<<": did not reorgainze GIDs, but probably should have\n";
        if(reorganized && (max_rounds == 1)) std::cout<<comm->getRank()<<": reorganized GIDs, but probably should not have\n";
      }
      //remap the GIDs so we receive the adjacencies in the same order as the current processes LIDs
      for (gno_t i = 0; i < sendcount; i++){
        ownedPlusGhosts[i+offsets.size()-1] = final_gid_vec[i];
      }
      if(verbose) {
        std::cout<<comm->getRank()<<": done remapping\n";
        std::cout<<comm->getRank()<<": building ghost offsets\n";
      }
      std::vector<offset_t> ghost_offsets(sendcount+1,0);
      std::vector<lno_t> send_adjs(adj_len,0);
      for(int i = 1; i < sendcount+1; i++){
        ghost_offsets[i] = ghost_offsets[i-1] + final_degree_vec[i-1];
      }

      if(verbose) std::cout<<comm->getRank()<<": building send_adjs\n";
      offset_t adjidx = 0;
      for(gno_t i = 0; i < recvcounttotal; i++){
        lno_t lid = mapOwned->getLocalElement(recvbuf[i]);
        for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
          send_adjs[adjidx++] = adjs[j];
        }
      }

      if(verbose) std::cout<<comm->getRank()<<": building recvadjscount\n";
      offset_t recvadjscount = 0;
      for(int i = 0; i < recvDegrees.size(); i++){
        recvadjscount += final_degree_vec[i];
      }
      if(verbose) std::cout<<comm->getRank()<<": going through the sending rounds\n";

      std::vector<uint64_t> curr_idx_per_proc(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++) curr_idx_per_proc[i] = rdispls[i];
      for(int round = 0; round < max_rounds; round++){
        std::vector<gno_t> send_adj;
        std::vector<int> send_adj_counts(comm->getSize(),0);
        if(verbose) std::cout<<comm->getRank()<<": round "<<round<<", constructing send_adj\n";
        for(int curr_proc = 0; curr_proc < comm->getSize(); curr_proc++){
          uint64_t curr_adj_sum = 0;
          while( curr_idx_per_proc[curr_proc] < (size_t)rdispls[curr_proc+1]){
            lno_t lid = mapOwned->getLocalElement(recvbuf[curr_idx_per_proc[curr_proc]++]);
            if((curr_adj_sum + (offsets[lid+1]-offsets[lid]))*sizeof(gno_t) >= INT_MAX){
              break;
            }
            curr_adj_sum += (offsets[lid+1] - offsets[lid]);
            for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
              send_adj.push_back(adjs[j]);
            }
          }
          send_adj_counts[curr_proc] = curr_adj_sum;
        }
        if(verbose) std::cout<<comm->getRank()<<": round "<<round<<", sending...\n";
        Teuchos::ArrayView<gno_t> send_adjs_view = Teuchos::arrayViewFromVector(send_adj);
        Teuchos::ArrayView<int> adjsendcounts_view = Teuchos::arrayViewFromVector(send_adj_counts);
        Teuchos::ArrayRCP<gno_t> ghost_adjs;
        uint64_t recv_adjs_count = 0;
        for(int i = 0; i < comm->getSize(); i++){
          recv_adjs_count += per_proc_round_adj_sums[round][i];
        }
        std::vector<int> adjrecvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> adjsrecvcounts_view = Teuchos::arrayViewFromVector(adjrecvcounts);
        Zoltan2::AlltoAllv<gno_t>(*comm, *env, send_adjs_view, adjsendcounts_view, ghost_adjs, adjsrecvcounts_view);
        for(offset_t i = 0; i< (offset_t)ghost_adjs.size(); i++){
          adjs_2GL.push_back(ghost_adjs[i]);
        }
      }
      if(verbose) std::cout<<comm->getRank()<<": constructing offsets\n";
      for(int i = 0; i < sendcount+1; i++){
        offsets_2GL.push_back(ghost_offsets[i]);
      }
      if(verbose) std::cout<<comm->getRank()<<": done building 2nd ghost layer\n";
    }
    
    double doOwnedToGhosts(RCP<const map_t> mapOwnedPlusGhosts,
                           size_t nVtx,
			   Kokkos::View<lno_t*,device_type> verts_to_send,
			   Kokkos::View<size_t[1],device_type> verts_to_send_size,
                           Kokkos::View<int*, device_type>& colors,
                           gno_t& total_sent, gno_t& total_recvd){
      int nprocs = comm->getSize();
      std::vector<int> sendcnts(comm->getSize(), 0);
      std::vector<gno_t> sdispls(comm->getSize()+1, 0);
   
      for(size_t i = 0; i < verts_to_send_size(0); i++){
        for(size_t j = 0; j < procs_to_send[verts_to_send(i)].size(); j++){
	  sendcnts[procs_to_send[verts_to_send(i)][j]] += 2;
	}
      }
      
      sdispls[0] = 0;
      gno_t sendsize = 0;
      std::vector<int> sentcount(nprocs, 0);

      for(int i = 1; i < comm->getSize()+1; i++){
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
	sendsize += sendcnts[i-1];
      }
      total_sent = sendsize;
      std::vector<gno_t> sendbuf(sendsize,0);

      for(size_t i = 0; i < verts_to_send_size(0); i++){
        std::vector<int> procs = procs_to_send[verts_to_send(i)];
	for(size_t j = 0; j < procs.size(); j++){
	  size_t idx = sdispls[procs[j]] + sentcount[procs[j]];
	  sentcount[procs[j]] += 2;
	  sendbuf[idx++] = mapOwnedPlusGhosts->getGlobalElement(verts_to_send(i));
	  sendbuf[idx] = colors(verts_to_send(i));
	}
      }

      Teuchos::ArrayView<int> sendcnts_view = Teuchos::arrayViewFromVector(sendcnts);
      Teuchos::ArrayView<gno_t> sendbuf_view = Teuchos::arrayViewFromVector(sendbuf);
      Teuchos::ArrayRCP<gno_t> recvbuf;
      std::vector<int> recvcnts(comm->getSize(), 0);
      Teuchos::ArrayView<int> recvcnts_view = Teuchos::arrayViewFromVector(recvcnts);

      //if we're computing statistics, remove the computation imbalance from the comm timer
      if(verbose) comm->barrier();
      double comm_total = 0.0;
      double comm_temp = timer();

      Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcnts_view, recvbuf, recvcnts_view);
      comm_total += timer() - comm_temp;

      gno_t recvsize = 0;
      for(int i = 0; i < recvcnts_view.size(); i++){
        recvsize += recvcnts_view[i];
      }
      total_recvd = recvsize;
      for(int i = 0; i < recvsize; i+=2){
        size_t lid = mapOwnedPlusGhosts->getLocalElement(recvbuf[i]);
	if(lid < nVtx && verbose) std::cout<<comm->getRank()<<": received a locally owned vertex, somehow\n";
	colors(lid) = recvbuf[i+1];
      }

      return comm_total;
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
      model->getEdgeList(adjs, offsets, ewgts);      
      

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
      mapOwned->getRemoteIndexList(gids, owningProcs);
      
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
      if(verbose) std::cout<<comm->getRank()<<": constructing Tpetra stuff\n";
      dummy = Teuchos::OrdinalTraits <Tpetra::global_size_t>::invalid();
      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy,
                                           Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                           0, comm));
      if(verbose) std::cout<<comm->getRank()<<": done constructing map with copies\n";

      using import_t = Tpetra::Import<lno_t, gno_t>;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                            mapWithCopies));
      if(verbose) std::cout<<comm->getRank()<<": done constructing importer\n";
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned,
                                                    importer, 1, true));
      if(verbose) std::cout<<comm->getRank()<<": done constructing femv\n";

      //Create random numbers seeded on global IDs, as in AlgHybridGMB.
      //This may or may not have an effect on this algorithm, but we
      //might as well see.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(size_t i = 0; i < rand.size(); i++){
        std::srand(ownedPlusGhosts[i]);
        rand[i] = std::rand();
      }

      
      std::vector<int> ghostOwners2(ownedPlusGhosts.size() -nVtx);
      std::vector<gno_t> ghosts2(ownedPlusGhosts.size() - nVtx);
      for(size_t i = nVtx; i < ownedPlusGhosts.size(); i++) ghosts2[i-nVtx] = ownedPlusGhosts[i];
      Teuchos::ArrayView<int> owners2 = Teuchos::arrayViewFromVector(ghostOwners2);
      Teuchos::ArrayView<const gno_t> ghostGIDs = Teuchos::arrayViewFromVector(ghosts2);
      mapOwned->getRemoteIndexList(ghostGIDs,owners2);
      if(verbose) std::cout<<comm->getRank()<<": done getting ghost owners\n";
      
      //calculations for how many 2GL verts are in the boundary of another process, only
      //do this if it's going to be displayed.
      if(verbose) {
        std::cout<<comm->getRank()<<": calculating 2GL stats...\n";
        
        std::vector<int> sendcounts(comm->getSize(),0);
        std::vector<gno_t> sdispls(comm->getSize()+1,0);
        //loop through owners, count how many vertices we'll send to each processor
        for(int i = nGhosts; i < ghostGIDs.size(); i++){
          if(owners2[i] != comm->getRank()&& owners2[i] !=-1) sendcounts[owners2[i]]++;
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
        for(size_t i = nGhosts; i < (size_t)ghostGIDs.size(); i++){
          if(owners2[i] != comm->getRank() && owners2[i] != -1){
            sendbuf[idx[owners2[i]]++] = ghostGIDs[i];
          }
        }
        Teuchos::ArrayView<int> sendcounts_view = Teuchos::arrayViewFromVector(sendcounts);
        Teuchos::ArrayView<gno_t> sendbuf_view = Teuchos::arrayViewFromVector(sendbuf);
        Teuchos::ArrayRCP<gno_t>  recvbuf;
        std::vector<int> recvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> recvcounts_view = Teuchos::arrayViewFromVector(recvcounts);
        Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcounts_view, recvbuf, recvcounts_view);
        std::vector<int> is_bndry_send(recvbuf.size(),0);
        for(int i = 0; i < recvbuf.size(); i++){
          size_t lid = mapWithCopies->getLocalElement(recvbuf[i]);
          is_bndry_send[i] = 0;
          if(lid < nVtx){
            for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
              if((size_t)local_adjs[j] >= nVtx) is_bndry_send[i] = 1;
            }
          } else{
            for(offset_t j = first_layer_ghost_offsets[lid]; j < first_layer_ghost_offsets[lid+1]; j++){
              if((size_t)local_ghost_adjs[j] >= nVtx) is_bndry_send[i] = 1;
            }
          }
        }
        
        Teuchos::ArrayView<int> is_bndry_send_view = Teuchos::arrayViewFromVector(is_bndry_send);
        Teuchos::ArrayRCP<int> is_bndry_recv;
        std::vector<int> bndry_recvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> bndry_recvcounts_view = Teuchos::arrayViewFromVector(bndry_recvcounts);
        Zoltan2::AlltoAllv<int> (*comm, *env, is_bndry_send_view, recvcounts_view, is_bndry_recv, bndry_recvcounts_view);

        int boundaryverts = 0;
        for(int i = 0; i < is_bndry_recv.size(); i++){
          boundaryverts+= is_bndry_recv[i];
        }
        std::cout<<comm->getRank()<<": "<<boundaryverts<<" boundary verts out of "<<n2Ghosts<<" verts in 2GL\n";
      }

      Teuchos::ArrayView<const lno_t> local_adjs_view = Teuchos::arrayViewFromVector(local_adjs);
      Teuchos::ArrayView<const offset_t> ghost_offsets = Teuchos::arrayViewFromVector(first_layer_ghost_offsets);
      Teuchos::ArrayView<const lno_t> ghost_adjacencies = Teuchos::arrayViewFromVector(local_ghost_adjs);

      Teuchos::ArrayView<const lno_t> exportLIDs = importer->getExportLIDs();
      Teuchos::ArrayView<const int> exportPIDs = importer->getExportPIDs();

      for(int i = 0; i < exportLIDs.size(); i++){
        procs_to_send[exportLIDs[i]].push_back(exportPIDs[i]);
      }

      //call the coloring algorithm
      twoGhostLayer(nVtx, nVtx+nGhosts, local_adjs_view, offsets, ghost_adjacencies, ghost_offsets,
                    femv, ownedPlusGhosts, globalToLocal, rand, owners2, mapWithCopies);
      
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
                       const std::vector<int>& rand,
                       ArrayView<int> owners,
                       RCP<const map_t> mapOwnedPlusGhosts){
      double total_time = 0.0;
      double interior_time = 0.0;
      double comm_time = 0.0;
      double comp_time = 0.0;
      double recoloring_time=0.0;
      double conflict_detection = 0.0;

      //get the degrees of all ghost vertices
      std::vector<int> deg_send_cnts(comm->getSize(),0);
      std::vector<gno_t> deg_sdispls(comm->getSize()+1,0);
      for(int i = 0; i < owners.size(); i++){
        deg_send_cnts[owners[i]]++;
      }
      deg_sdispls[0] = 0;
      gno_t deg_sendsize = 0;
      std::vector<int> deg_sentcount(comm->getSize(),0);
      for(int i = 1; i < comm->getSize()+1; i++){
        deg_sdispls[i] = deg_sdispls[i-1] + deg_send_cnts[i-1];
	deg_sendsize += deg_send_cnts[i-1];
      }
      std::vector<gno_t> deg_sendbuf(deg_sendsize,0);
      for(int i = 0; i < owners.size(); i++){
        size_t idx = deg_sdispls[owners[i]] + deg_sentcount[owners[i]];
	deg_sentcount[owners[i]]++;
	deg_sendbuf[idx] = mapOwnedPlusGhosts->getGlobalElement(i+n_local);
      }
      Teuchos::ArrayView<int> deg_send_cnts_view = Teuchos::arrayViewFromVector(deg_send_cnts);
      Teuchos::ArrayView<gno_t> deg_sendbuf_view = Teuchos::arrayViewFromVector(deg_sendbuf);
      Teuchos::ArrayRCP<gno_t> deg_recvbuf;
      std::vector<int> deg_recvcnts(comm->getSize(),0);
      Teuchos::ArrayView<int> deg_recvcnts_view = Teuchos::arrayViewFromVector(deg_recvcnts);
      AlltoAllv<gno_t>(*comm, *env, deg_sendbuf_view, deg_send_cnts_view, deg_recvbuf, deg_recvcnts_view);

      //replace GID with local degree
      for(int i = 0; i < deg_recvbuf.size(); i++){
        lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_recvbuf[i]);
	deg_recvbuf[i] = offsets[lid+1] - offsets[lid];
      }
      //send modified buffer back
      ArrayRCP<gno_t> ghost_degrees;
      AlltoAllv<gno_t>(*comm, *env, deg_recvbuf(), deg_recvcnts_view, ghost_degrees, deg_send_cnts_view);

      //create the ghost degree views
      Kokkos::View<gno_t*, device_type> ghost_degrees_dev("ghost degree view",ghost_degrees.size());
      typename Kokkos::View<gno_t*, device_type>::HostMirror ghost_degrees_host = Kokkos::create_mirror(ghost_degrees_dev);
      for(int i = 0; i < ghost_degrees.size(); i++){
        lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_sendbuf[i]);
	ghost_degrees_host(lid-n_local) = ghost_degrees[i];
      }
      Kokkos::deep_copy(ghost_degrees_dev, ghost_degrees_host);

      gno_t recvPerRound[100];
      gno_t sentPerRound[100];
      //find global max degree
      offset_t local_max_degree = 0;
      offset_t global_max_degree = 0;
      for(size_t i = 0; i < n_local; i++){
	offset_t curr_degree = offsets[i+1] - offsets[i];
	if(curr_degree > local_max_degree){
	  local_max_degree = curr_degree;
	}
      }
      Teuchos::reduceAll<int, offset_t>(*comm, Teuchos::REDUCE_MAX,1, &local_max_degree, &global_max_degree);
      if(comm->getRank() == 0 && verbose) std::cout<<"Input has max degree "<<global_max_degree<<"\n";
      
      if(verbose) std::cout<<comm->getRank()<<": constructing Kokkos Views for initial coloring\n";
      Kokkos::View<offset_t*, device_type> offsets_dev("Host Offset View", offsets.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror offsets_host = Kokkos::create_mirror(offsets_dev);
      Kokkos::View<lno_t*, device_type> adjs_dev("Host Adjacencies View", adjs.size());
      typename Kokkos::View<lno_t*, device_type>::HostMirror adjs_host = Kokkos::create_mirror(adjs_dev);
      for(Teuchos_Ordinal i = 0; i < offsets.size(); i++) offsets_host(i) = offsets[i];
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) adjs_host(i) = adjs[i]; 
      Kokkos::deep_copy(offsets_dev,offsets_host);
      Kokkos::deep_copy(adjs_dev, adjs_host);

      
      //create the graph structures which allow KokkosKernels to recolor the conflicting vertices
      if(verbose) std::cout<<comm->getRank()<<": constructing Kokkos Views for recoloring\n";
      Kokkos::View<offset_t*, device_type> dist_degrees_dev("Owned+Ghost degree view",rand.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_degrees_host = Kokkos::create_mirror(dist_degrees_dev);
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) dist_degrees_host(adjs[i])++;
      for(Teuchos_Ordinal i = 0; i < ghost_adjs.size(); i++) dist_degrees_host(ghost_adjs[i])++;
      for(Teuchos_Ordinal i = 0; i < offsets.size()-1; i++) dist_degrees_host(i) = offsets[i+1] - offsets[i];
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++) dist_degrees_host(i+n_local) = ghost_offsets[i+1] - ghost_offsets[i];
      
      Kokkos::View<offset_t*, device_type> dist_offsets_dev("Owned+Ghost Offset view", rand.size()+1);
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host = Kokkos::create_mirror(dist_offsets_dev);
      dist_offsets_host(0) = 0;
      offset_t total_adjs = 0;
      for(size_t i = 1; i < rand.size()+1; i++){
        dist_offsets_host(i) = dist_degrees_host(i-1) + dist_offsets_host(i-1);
        total_adjs += dist_degrees_host(i-1);
      }
      Kokkos::View<lno_t*, device_type> dist_adjs_dev("Owned+Ghost adjacency view", total_adjs);
      typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host = Kokkos::create_mirror(dist_adjs_dev);
      for(size_t i = 0; i < rand.size(); i++){
        dist_degrees_host(i) = 0;
      }
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) dist_adjs_host(i) = adjs[i];
      for(Teuchos_Ordinal i = adjs.size(); i < adjs.size() + ghost_adjs.size(); i++) dist_adjs_host(i) = ghost_adjs[i-adjs.size()];
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++){
        for(offset_t j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
          if((size_t)ghost_adjs[j] >= n_total){
            dist_adjs_host(dist_offsets_host(ghost_adjs[j]) + dist_degrees_host(ghost_adjs[j])) = i + n_local;
            dist_degrees_host(ghost_adjs[j])++;
          }
        }
      }

      Kokkos::deep_copy(dist_degrees_dev,dist_degrees_host);
      Kokkos::deep_copy(dist_offsets_dev,dist_offsets_host);
      Kokkos::deep_copy(dist_adjs_dev, dist_adjs_host);
      //we can find all the conflicts with one loop through the ghost vertices.
      
      //this view represents how many conflicts were found
      Kokkos::View<gno_t[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      //keep an atomic version so that we can increment from multiple threads
      Kokkos::View<gno_t[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize;

      //create views for the ghost adjacencies, as they can detect all conflicts.
      if(verbose) std::cout<<comm->getRank()<<": constructing ghost connectivity views\n";
      Kokkos::View<offset_t*, device_type> ghost_offset_dev("Ghost Offsets", ghost_offsets.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror ghost_offset_host = Kokkos::create_mirror(ghost_offset_dev);
      Kokkos::View<lno_t*, device_type> ghost_adjs_dev("Ghost Adjacencies", ghost_adjs.size());
      typename Kokkos::View<lno_t*, device_type>::HostMirror ghost_adjs_host = Kokkos::create_mirror(ghost_adjs_dev);
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size(); i++) ghost_offset_host(i) = ghost_offsets[i];
      for(Teuchos_Ordinal i = 0; i < ghost_adjs.size(); i++) ghost_adjs_host(i) = ghost_adjs[i];

      Kokkos::deep_copy(ghost_offset_dev, ghost_offset_host);
      Kokkos::deep_copy(ghost_adjs_dev, ghost_adjs_host);

      //create a view for the tie-breaking numbers.
      if(verbose) std::cout<<comm->getRank()<<": constructing rand and GIDs views\n";
      Kokkos::View<int*, device_type> rand_dev("Random View", rand.size());
      typename Kokkos::View<int*, device_type>::HostMirror rand_host = Kokkos::create_mirror(rand_dev);
      for(size_t i = 0; i < rand.size(); i ++) rand_host(i) = rand[i];
      Kokkos::deep_copy(rand_dev,rand_host);
      
      Kokkos::View<gno_t*, device_type> gid_dev("GIDs", gids.size());
      typename Kokkos::View<gno_t*, device_type>::HostMirror gid_host = Kokkos::create_mirror(gid_dev);
      for(size_t i = 0; i < gids.size(); i++) gid_host(i) = gids[i];
      Kokkos::deep_copy(gid_dev,gid_host);
       
      if(verbose) std::cout<<comm->getRank()<<": counting boundary\n";
      offset_t boundary_size = 0;
      for(size_t i = 0; i < n_local; i++){
        for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
	  if((size_t)dist_adjs_host(j) >= n_local){
	    boundary_size++;
	    break;
	  }
	  bool found = false;
	  for(offset_t k = dist_offsets_host(dist_adjs_host(j)); k < dist_offsets_host(dist_adjs_host(j)+1); k++){
	    if((size_t)dist_adjs_host(k) >= n_local){
	      boundary_size++;
	      found = true;
	      break;
	    }
	  }
	  if(found) break;
	}
      }
      Kokkos::View<bool*, device_type> send_flags("send flags",n_local);
      Kokkos::parallel_for(n_local,KOKKOS_LAMBDA(const int& i){
	send_flags(i) = false;	      
      });

      if(verbose) std::cout<<comm->getRank()<<": constructing communication and recoloring lists\n";
      Kokkos::View<lno_t*, device_type> verts_to_recolor_view("verts to recolor", rand.size());
      Kokkos::parallel_for(boundary_size, KOKKOS_LAMBDA(const int& i){
	verts_to_recolor_view(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_atomic = verts_to_recolor_view;
      
      Kokkos::View<int[1], device_type> verts_to_recolor_size("verts to recolor size");
      verts_to_recolor_size(0) = 0;
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_size_atomic = verts_to_recolor_size;

      Kokkos::View<lno_t*, device_type> verts_to_send_view("verts to send", boundary_size);
      Kokkos::parallel_for(boundary_size, KOKKOS_LAMBDA(const int& i){
        verts_to_send_view(i) = -1;		      
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_atomic = verts_to_send_view;
      

      Kokkos::View<size_t[1], device_type> verts_to_send_size("verts to send size");
      verts_to_send_size(0) = 0;
      Kokkos::View<size_t[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic = verts_to_send_size;
      if(verbose) std::cout<<comm->getRank()<<": initializing the list of vertices to send\n";
      Kokkos::parallel_for(n_local, KOKKOS_LAMBDA(const int& i){
        for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
	  if((size_t)dist_adjs_dev(j) >= n_local){
	    verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	    send_flags(i) = true;
	    break;
	  }
	  bool found = false;
	  for(offset_t k = dist_offsets_dev(dist_adjs_dev(j)); k < dist_offsets_dev(dist_adjs_dev(j)+1); k++){
	    if((size_t)dist_adjs_dev(k) >= n_local){
	      verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	      found = true;
	      send_flags(i) = true;
	      break;
	    }
	  }
	  if(found) break;
	}		      
      });
      Kokkos::fence();
      
      bool use_vbbit = (global_max_degree < 6000);
      //Done initializing, start coloring!
      comm->barrier();
      interior_time = timer();
      total_time = timer();
      //give the entire local graph to KokkosKernels to color
      //this->colorInterior(n_local, adjs_dev, offsets_dev, femv,verts_to_recolor_view,verts_to_recolor_size(0),use_vbbit);
      this->colorInterior<execution_space, memory_space, memory_space>
	      (n_local, adjs_dev, offsets_dev, femv,adjs_dev,0,use_vbbit);
      interior_time = timer() - interior_time;
      comp_time = interior_time;
      //get the color view from the FEMultiVector
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
      Kokkos::View<int*,device_type> ghost_colors("ghost color backups", rand.size() - n_local);
      comm->barrier();

      //communicate the initial coloring.
      if(verbose) std::cout<<comm->getRank()<<": communicating\n";
      gno_t recv,sent;
      comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,n_local,verts_to_send_view,verts_to_send_size,femv_colors,sent,recv);
      sentPerRound[0] = sent;
      recvPerRound[0] = recv; 

      Kokkos::parallel_for(rand.size() - n_local, KOKKOS_LAMBDA(const int& i){
        ghost_colors(i) = femv_colors(i+n_local);
      });
      Kokkos::fence();
      
      verts_to_send_size(0) = 0;
      double temp = timer();
      //detect conflicts only for ghost vertices
      bool recolor_degrees = this->pl->template get<bool>("recolor_degrees",false);
      if(verbose) std::cout<<comm->getRank()<<": detecting conflicts\n";
      Kokkos::RangePolicy<execution_space> policy(n_local,rand.size());
      Kokkos::parallel_reduce("conflict detection",policy, KOKKOS_LAMBDA (const int& i, gno_t& recoloring_size){
        lno_t localIdx = i;
        int currColor = femv_colors(localIdx);
	int currDegree = ghost_degrees_dev(i-n_local);
        for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
          int nborColor = femv_colors(dist_adjs_dev(j));
	  int nborDegree = 0;
	  if((size_t)dist_adjs_dev(j) < n_local) nborDegree = dist_offsets_dev(dist_adjs_dev(j)+1) - dist_offsets_dev(dist_adjs_dev(j));
	  else nborDegree = ghost_degrees_dev(dist_adjs_dev(j) - n_local);
          if(currColor == nborColor ){
	    if(currDegree < nborDegree && recolor_degrees){
	      femv_colors(localIdx) = 0;
	      recoloring_size++;
	      break;
	    }else if(nborDegree < currDegree && recolor_degrees){
	      femv_colors(dist_adjs_dev(j)) = 0;
              recoloring_size++;
	    }else if(rand_dev(localIdx) > rand_dev(dist_adjs_dev(j))){
              recoloring_size++;
              femv_colors(localIdx) = 0;
              break;
            }else if(rand_dev(dist_adjs_dev(j)) > rand_dev(localIdx)){
                recoloring_size++;
                femv_colors(dist_adjs_dev(j)) = 0;
            } else {
              if (gid_dev(localIdx) >= gid_dev(dist_adjs_dev(j))){
                femv_colors(localIdx) = 0;
                recoloring_size++;
                break;
              } else {
                  femv_colors(dist_adjs_dev(j)) = 0;
                  recoloring_size++;
              }
            }
          }
        }
      },recoloringSize(0));
      Kokkos::fence();
      if(comm->getSize() > 1){
        conflict_detection = timer() - temp;
        comp_time += conflict_detection;
      }
      Kokkos::parallel_for(femv_colors.size(), KOKKOS_LAMBDA (const size_t& i){
        if(femv_colors(i) == 0){
          if(i < n_local){
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
          } 
          verts_to_recolor_atomic(verts_to_recolor_size_atomic(0)++) = i;
        }
      });

      //ensure that the parallel_for finishes before continuing
      Kokkos::fence();
      
      //all conflicts detected!
      if(verbose) std::cout<<comm->getRank()<<": starting to recolor\n";
      //variables for statistics
      double totalPerRound[100];
      double commPerRound[100];
      double compPerRound[100];
      double recoloringPerRound[100];
      double conflictDetectionPerRound[100];
      uint64_t vertsPerRound[100];
      uint64_t incorrectGhostsPerRound[100];
      int distributedRounds = 1;
      totalPerRound[0] = interior_time + comm_time + conflict_detection;
      recoloringPerRound[0] = 0;
      commPerRound[0] = comm_time;
      compPerRound[0] = interior_time + conflict_detection;
      conflictDetectionPerRound[0] = conflict_detection;
      recoloringPerRound[0] = 0;
      vertsPerRound[0] = 0;
      incorrectGhostsPerRound[0]=0;
      typename Kokkos::View<int*, device_type>::HostMirror colors_host;
      typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_recolor_host;
      typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send_host;
      typename Kokkos::View<int*, device_type>::HostMirror ghost_colors_host;
      int serial_threshold = this->pl->template get<int>("serial_threshold",0);
      //see if recoloring is necessary.
      gno_t totalConflicts = 0;
      gno_t localConflicts = recoloringSize(0);
      Teuchos::reduceAll<int,gno_t>(*comm, Teuchos::REDUCE_SUM, 1, &localConflicts, &totalConflicts);
      bool done = !totalConflicts;
      if(comm->getSize()==1) done = true;
      
      //recolor until no conflicts are left
      while(!done){
	if(recoloringSize(0) < serial_threshold) break;
        if(distributedRounds < 100) {
          vertsPerRound[distributedRounds] = verts_to_recolor_size(0);
        }
        
        comm->barrier();
        double recolor_temp = timer();
        //recolor using KokkosKernels' coloring function 
        if(verts_to_recolor_size(0) > 0){
	  this->colorInterior<execution_space, memory_space, memory_space>
		  (femv_colors.size(), dist_adjs_dev, dist_offsets_dev,femv,verts_to_recolor_view,verts_to_recolor_size(0),use_vbbit);
	}
        recoloringPerRound[distributedRounds] = timer() - recolor_temp;
        recoloring_time += recoloringPerRound[distributedRounds];
        total_time += recoloringPerRound[distributedRounds];
        comp_time += recoloringPerRound[distributedRounds];
        compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        recoloringSize(0) = 0;
        verts_to_recolor_size(0) = 0;
        Kokkos::parallel_for(rand.size() - n_local, KOKKOS_LAMBDA(const int& i){
          femv_colors(i+n_local) = ghost_colors(i);
        });
        Kokkos::fence();
        /*std::vector<int> old_colors;
	//only loop through these if we're going to use them
	if(verbose){
          for(size_t i = 0; i < femv_colors.size(); i++){
            old_colors.push_back(femv_colors(i));
          }
	}*/
        comm->barrier();
        //communicate the new colors
        commPerRound[distributedRounds] = doOwnedToGhosts(mapOwnedPlusGhosts,n_local,verts_to_send_view,verts_to_send_size,femv_colors,sent,recv);
        recvPerRound[distributedRounds] = recv;
        sentPerRound[distributedRounds] = sent;
	if(verbose) {
          std::cout<<comm->getRank()<<": total sent in round "<<distributedRounds<<" = "<<sent<<"\n";
          std::cout<<comm->getRank()<<": total recv in round "<<distributedRounds<<" = "<<recv<<"\n";
	}
        //incorrectGhostsPerRound[distributedRounds] = 0;
        //only compute this metric if it will be displayed
	/*if(verbose){
          for(size_t i = 0; i < femv_colors.size(); i++){
            if(i >= n_local && femv_colors(i) != old_colors[i]) incorrectGhostsPerRound[distributedRounds]++;
          }
	}*/
        comm_time += commPerRound[distributedRounds];
        totalPerRound[distributedRounds] += commPerRound[distributedRounds];
        
	verts_to_send_size(0) = 0;
        Kokkos::parallel_for(rand.size() - n_local, KOKKOS_LAMBDA(const int& i){
          ghost_colors(i) = femv_colors(i+n_local);
        });
        Kokkos::fence();
        comm->barrier();
        //check for further conflicts
        double detection_temp = timer();
        Kokkos::parallel_reduce("conflict detection2",policy, KOKKOS_LAMBDA (const int& i, gno_t& recoloring_size){
          lno_t localIdx = i;
          int currColor = femv_colors(localIdx);
	  int currDegree = ghost_degrees_dev(i-n_local);
          for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
            int nborColor = femv_colors(dist_adjs_dev(j));
	    int nborDegree = 0;
	    if((size_t)dist_adjs_dev(j) < n_local) nborDegree = dist_offsets_dev(dist_adjs_dev(j)+1) - dist_offsets_dev(dist_adjs_dev(j));
	    else nborDegree = ghost_degrees_dev(dist_adjs_dev(j) - n_local);
            if(currColor == nborColor ){
	      if(currDegree < nborDegree && recolor_degrees){
	        femv_colors(localIdx) = 0;
	        recoloring_size++;
	        break;
	      }else if(nborDegree < currDegree && recolor_degrees){
	        femv_colors(dist_adjs_dev(j)) = 0;
                recoloring_size++;
	      }else if(rand_dev(localIdx) > rand_dev(dist_adjs_dev(j))){
                recoloring_size++;
                femv_colors(localIdx) = 0;
                break;
              }else if(rand_dev(dist_adjs_dev(j)) > rand_dev(localIdx)){
                  recoloring_size++;
                  femv_colors(dist_adjs_dev(j)) = 0;
              } else {
                if (gid_dev(localIdx) >= gid_dev(dist_adjs_dev(j))){
                  femv_colors(localIdx) = 0;
                  recoloring_size++;
                  break;
                } else {
                    femv_colors(dist_adjs_dev(j)) = 0;
                    recoloring_size++;
                }
              }
            }
          }
        },recoloringSize(0));
        Kokkos::fence();
        Kokkos::parallel_for(femv_colors.size(), KOKKOS_LAMBDA (const size_t& i){
          if(femv_colors(i) == 0){
            if(i < n_local){
              verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
            }
            verts_to_recolor_atomic(verts_to_recolor_size_atomic(0)++) = i;
          }
        });

        //ensure the parallel_for finishes before continuing
        Kokkos::fence(); 
        
	conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
        conflict_detection += conflictDetectionPerRound[distributedRounds];
        compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        comp_time += conflictDetectionPerRound[distributedRounds];

        distributedRounds++;
        int localDone = recoloringSize(0);
        int globalDone = 0;
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
        done = !globalDone;
        
      }//end device coloring

      if(recoloringSize(0) > 0 || !done){
	colors_host = Kokkos::create_mirror_view(femv_colors);
	deep_copy(colors_host, femv_colors);
	ghost_colors_host = Kokkos::create_mirror_view(ghost_colors);
	deep_copy(ghost_colors_host, ghost_colors);
        verts_to_recolor_host = Kokkos::create_mirror_view(verts_to_recolor_view);
	deep_copy(verts_to_recolor_host, verts_to_recolor_view);
	verts_to_send_host = Kokkos::create_mirror_view(verts_to_send_view);
	deep_copy(verts_to_send_host, verts_to_send_view);
      }

      while(recoloringSize(0) > 0 || !done){
	if(distributedRounds < 100){
	  vertsPerRound[distributedRounds] = recoloringSize(0);
	}
	if(verbose) std::cout<<comm->getRank()<<": starting to recolor, serial\n";
        comm->barrier();
	double recolor_temp = timer();
	if(verts_to_recolor_size(0) > 0){
	  this->colorInterior<Kokkos::DefaultHostExecutionSpace,
		              memory_space,
			      memory_space>
		  (femv_colors.size(), dist_adjs_host, dist_offsets_host, femv, verts_to_recolor_host, verts_to_recolor_size(0), true);
	}
	recoloringPerRound[distributedRounds] = timer() - recolor_temp;
	recoloring_time += recoloringPerRound[distributedRounds];
	total_time += recoloringPerRound[distributedRounds];
	comp_time += recoloringPerRound[distributedRounds];
	compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        recoloringSize(0) = 0;
	verts_to_recolor_size(0) = 0;
	for(size_t i = 0; i < rand.size() - n_local; i++){
	  femv_colors(i+n_local) = ghost_colors_host(i);
	}
        comm->barrier();
	commPerRound[distributedRounds] = doOwnedToGhosts(mapOwnedPlusGhosts, n_local,verts_to_send_view, verts_to_send_size, femv_colors, sent,recv);
	recvPerRound[distributedRounds] = recv;
	sentPerRound[distributedRounds] = sent;
	if(verbose) {
	  std::cout<<comm->getRank()<<": total sent in round "<<distributedRounds<<" = "<<sent<<"\n";
	  std::cout<<comm->getRank()<<": total recv in round "<<distributedRounds<<" = "<<recv<<"\n";
	}
        comm_time += commPerRound[distributedRounds];
	totalPerRound[distributedRounds] += commPerRound[distributedRounds];
	verts_to_send_size(0) = 0;
	for(size_t i = 0; i < rand.size() - n_local; i++){
	  ghost_colors(i) = femv_colors(i+n_local);
	}
	comm->barrier();
	double detection_temp = timer();
	for(size_t i = n_local; i < rand.size(); i++){
	  lno_t localIdx = i;
	  int currColor = femv_colors(localIdx);
	  int currDegree = ghost_degrees_host(i-n_local);
	  for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
	    int nborColor = femv_colors(dist_adjs_host(j));
	    int nborDegree = 0;
	    if((size_t)dist_adjs_host(j) < n_local){
	      nborDegree = dist_offsets_host(dist_adjs_host(j)+1) - dist_offsets_host(dist_adjs_host(j));
	    } else {
	      nborDegree = ghost_degrees_host(dist_adjs_host(j) - n_local);
	    }
	    if(currColor == nborColor) {  
	      if(currDegree < nborDegree && recolor_degrees){
	        femv_colors(localIdx) = 0;
	        recoloringSize(0)++;
		//verts_to_recolor_host(verts_to_recolor_size(0)++) = localIdx;
	        break;
	      }else if(nborDegree < currDegree && recolor_degrees){
	        femv_colors(dist_adjs_host(j)) = 0;
                recoloringSize(0)++;
		//verts_to_recolor_host(verts_to_recolor_size(0)++) = dist_adjs_host(j);
		//if(dist_adjs_host(j) < n_local){
		//  verts_to_send_host(verts_to_send_size(0)++) = dist_adjs_host(j);
		//}
	      }else if(rand_host(localIdx) > rand_host(dist_adjs_host(j))){
	        femv_colors(localIdx) = 0;
		recoloringSize(0)++;
		//verts_to_recolor_host(verts_to_recolor_size(0)++) = localIdx;
		break;
	      } else if (rand_host(dist_adjs_host(j)) > rand_host(localIdx)) {
	        femv_colors(dist_adjs_host(j)) = 0;
		recoloringSize(0)++;
		//verts_to_recolor_host(verts_to_recolor_size(0)++) = dist_adjs_host(j);
		//if(dist_adjs_host(j) < n_local){
		//  verts_to_send_host(verts_to_send_size(0)++) = dist_adjs_host(j);
		//}
	      } else {
	        if(gid_host(localIdx) >= gid_host(dist_adjs_host(j))){
		  femv_colors(localIdx) = 0;
		  recoloringSize(0)++;
		  //verts_to_recolor_host(verts_to_recolor_size(0)++) = localIdx;
		  break;
		} else {
		  femv_colors(dist_adjs_host(j)) = 0;
		  recoloringSize(0)++;
		  //verts_to_recolor_host(verts_to_recolor_size(0)++) = dist_adjs_host(j);
		  //if(dist_adjs_host(j) < n_local){
		  //  verts_to_send_host(verts_to_send_size(0)++) = dist_adjs_host(j);
		  //}
		}
	      }
	    } 
	  }
	}
        for(size_t i = 0; i < femv_colors.size(); i++){
	  if(femv_colors(i) == 0){
	    if(i < n_local){
	      verts_to_send_host(verts_to_send_size(0)++) = i;
	    }
	    verts_to_recolor_host(verts_to_recolor_size(0)++) = i;
	  }
	} 
	conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
	conflict_detection += conflictDetectionPerRound[distributedRounds];
	compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
	totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
	comp_time += conflictDetectionPerRound[distributedRounds];

	int globalDone = 0;
	int localDone = recoloringSize(0);
	Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
	distributedRounds++;
	done = !globalDone;
      }

      total_time = timer() - total_time;
      //only compute stats if they will be displayed
      if(verbose){
        uint64_t localBoundaryVertices = 0;
        for(size_t i = 0; i < n_local; i++){
          for(offset_t j = offsets[i]; j < offsets[i+1]; j++){
            if((size_t)adjs[j] >= n_local){
              localBoundaryVertices++;
              break;
            }
          }
        }      
        
        if(comm->getRank() == 0) printf("did %d rounds of distributed coloring\n", distributedRounds);
        uint64_t totalVertsPerRound[100];
        uint64_t totalBoundarySize = 0;
        uint64_t totalIncorrectGhostsPerRound[100];
        double finalTotalPerRound[100];
        double maxRecoloringPerRound[100];
        double minRecoloringPerRound[100];
        double finalCommPerRound[100];
        double finalCompPerRound[100];
        double finalConflictDetectionPerRound[100];
        gno_t finalRecvPerRound[100];
        gno_t finalSentPerRound[100];
        for(int i = 0; i < 100; i++){
          totalVertsPerRound[i] = 0;
          finalTotalPerRound[i] = 0.0;
          maxRecoloringPerRound[i] = 0.0;
          minRecoloringPerRound[i] = 0.0;
          finalCommPerRound[i] = 0.0;
          finalCompPerRound[i] = 0.0;
          finalConflictDetectionPerRound[i] = 0.0;
          finalRecvPerRound[i] = 0;
          finalSentPerRound[i] = 0;
        }
        Teuchos::reduceAll<int,uint64_t>(*comm, Teuchos::REDUCE_SUM,1,&localBoundaryVertices, &totalBoundarySize);
        Teuchos::reduceAll<int,uint64_t>(*comm, Teuchos::REDUCE_SUM,100,vertsPerRound,totalVertsPerRound);
        Teuchos::reduceAll<int,uint64_t>(*comm, Teuchos::REDUCE_SUM,100,incorrectGhostsPerRound,totalIncorrectGhostsPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,totalPerRound, finalTotalPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,recoloringPerRound,maxRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,100,recoloringPerRound,minRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,commPerRound,finalCommPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,compPerRound,finalCompPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,conflictDetectionPerRound, finalConflictDetectionPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,recvPerRound,finalRecvPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,sentPerRound,finalSentPerRound);
        printf("Rank %d: boundary size: %ld\n",comm->getRank(),localBoundaryVertices);
        for(int i = 0; i < std::min((int)distributedRounds,100); i++){
          printf("Rank %d: recolor %ld vertices in round %d\n",comm->getRank(), vertsPerRound[i],i);
          printf("Rank %d: sentbuf had %lld entries in round %d\n", comm->getRank(), sentPerRound[i],i);
          if(comm->getRank()==0){
            printf("recolored %ld vertices in round %d\n",totalVertsPerRound[i], i);
            printf("%ld inconsistent ghosts in round %d\n",totalIncorrectGhostsPerRound[i],i);
            printf("total time in round %d: %f\n",i,finalTotalPerRound[i]);
            printf("recoloring time in round %d: %f\n",i,maxRecoloringPerRound[i]);
            printf("min recoloring time in round %d: %f\n",i,minRecoloringPerRound[i]);
            printf("conflict detection time in round %d: %f\n",i,finalConflictDetectionPerRound[i]);
            printf("comm time in round %d: %f\n",i,finalCommPerRound[i]);
            printf("recvbuf size in round %d: %lld\n",i,finalRecvPerRound[i]);
            printf("sendbuf size in round %d: %lld\n",i,finalSentPerRound[i]);
            printf("comp time in round %d: %f\n",i,finalCompPerRound[i]);
          }
        }
        double global_total_time = 0.0;
        double global_recoloring_time = 0.0;
        double global_min_recoloring_time = 0.0;
        double global_conflict_detection=0.0;
        double global_comm_time=0.0;
        double global_comp_time=0.0;
        double global_interior_time=0.0;
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&total_time,&global_total_time);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&recoloring_time,&global_recoloring_time);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,1,&recoloring_time,&global_min_recoloring_time);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&conflict_detection,&global_conflict_detection);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&comm_time,&global_comm_time);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&comp_time,&global_comp_time);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&interior_time,&global_interior_time);
        comm->barrier();
        fflush(stdout);
        if(comm->getRank()==0){
          printf("Boundary size: %ld\n",totalBoundarySize);
          printf("Total Time: %f\n",global_total_time);
          printf("Interior Time: %f\n",global_interior_time);
          printf("Recoloring Time: %f\n",global_recoloring_time);
          printf("Min Recoloring Time: %f\n",global_min_recoloring_time);
          printf("Conflict Detection Time: %f\n",global_conflict_detection);
          printf("Comm Time: %f\n",global_comm_time);
          printf("Comp Time: %f\n",global_comp_time);
        }
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

}//end namespace Zoltan2

#endif
