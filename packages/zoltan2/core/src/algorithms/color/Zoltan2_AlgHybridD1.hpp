#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <queue>
#include <sys/time.h>
#include <algorithm>

#include "Zoltan2_Algorithm.hpp"
#include "Zoltan2_GraphModel.hpp"
#include "Zoltan2_ColoringSolution.hpp"
#include "Zoltan2_Util.hpp"
#include "Zoltan2_TPLTraits.hpp"

#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"
#include <stdlib.h>

//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin and Manne

namespace Zoltan2 {

template <typename Adapter>
class AlgHybridGMB : public Algorithm<Adapter>
{
  public:
  
    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t, gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>;
    using device_type = Tpetra::Map<>::device_type;
    using execution_space = Tpetra::Map<>::execution_space;
    using memory_space = Tpetra::Map<>::memory_space;
    double timer() {
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
    }
    
  private:

    void buildModel(modelFlag_t &flags); 

    //function to invoke KokkosKernels distance-1 coloring    
    template <class ExecutionSpace, typename TempMemorySpace, 
              typename MemorySpace>
    void colorInterior(const size_t nVtx, 
                       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > adjs_view,
                       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > offset_view, 
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){
      
      //set up types to be used by KokkosKernels
      using KernelHandle =  KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, ExecutionSpace, MemorySpace, 
           MemorySpace>;
      using lno_row_view_t = Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;
      using lno_nnz_view_t = Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;

      KernelHandle kh;
      
      //pick which KokkosKernels algorithm to use. Typically based on max graph degree.
      if(recolor){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);  
      }
      
      //vertex_list_size indicates whether we have provided a list of vertices to recolor.
      //Only makes a difference if the algorithm to be used is VBBIT.
      if(vertex_list_size != 0){
        kh.get_graph_coloring_handle()->set_vertex_list(vertex_list,vertex_list_size);
      }

      kh.set_verbose(verbose);

      //set the initial coloring of the kh.get_graph_coloring_handle() to be
      //the data view from the femv.
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<MemorySpace>();
      Kokkos::View<int*, Tpetra::Map<>::device_type >  sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      kh.get_graph_coloring_handle()->set_tictoc(verbose);
      KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t>
                                                     (&kh, nVtx, nVtx, offset_view, adjs_view);
      numColors = kh.get_graph_coloring_handle()->get_num_colors();

      if(verbose){
        std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
      }
    }
    
    double doOwnedToGhosts(RCP<const map_t> mapOwnedPlusGhosts,
                         size_t nVtx,
			 typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send,
			 Kokkos::View<size_t[1], device_type>& verts_to_send_size,
                         Kokkos::View<int*, device_type>& colors,
                         gno_t& recv, gno_t& send){
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
      send = sendsize;
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

      //if we're computing statistics, remove the computation imbalance from the comm timer.
      if(verbose) comm->barrier();
      double comm_total = 0.0;
      double comm_temp = timer();

      Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcnts_view, recvbuf, recvcnts_view);
      comm_total += timer() - comm_temp;
      
      gno_t recvsize = 0;

      for(int i = 0; i < recvcnts_view.size(); i++){
        recvsize += recvcnts_view[i];
      }
      recv = recvsize;
      for(int i = 0; i < recvsize; i+=2){
        size_t lid = mapOwnedPlusGhosts->getLocalElement(recvbuf[i]);
	if(lid<nVtx && verbose) std::cout<<comm->getRank()<<": received a locally owned vertex, somehow\n";
	colors(lid) = recvbuf[i+1];
      }
       
      return comm_total;
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
    //constructor for the  hybrid distributed distance-1 algorithm
    AlgHybridGMB(
      const RCP<const base_adapter_t> &adapter_, 
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_) {
      verbose = pl->get<bool>("verbose",false);
      if(verbose) std::cout<<comm->getRank()<<": inside coloring constructor\n";
      numColors = 4;
      modelFlag_t flags;
      flags.reset();
      buildModel(flags);
      if(verbose) std::cout<<comm->getRank()<<": done constructing coloring class\n";
    }


    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution ) {
      if(verbose) std::cout<<comm->getRank()<<": inside coloring function\n"; 
      //this will color the global graph in a manner similar to Zoltan
      
      //get vertex GIDs in a locally indexed array
      if(verbose) std::cout<<comm->getRank()<<": getting owned vtxIDs\n";
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      //we do not use weights at this point
      if(verbose) std::cout<<comm->getRank()<<": getting edge list\n";
      //get edge information from the model
      ArrayView<const gno_t> adjs;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      model->getEdgeList(adjs, offsets, ewgts);
      //again, weights are not used
      
      RCP<const map_t> mapOwned;
      RCP<const map_t> mapWithCopies;
      
      std::vector<gno_t> finalGIDs;
      std::vector<offset_t> finalOffset_vec;
      std::vector<lno_t> finalAdjs_vec;      

      lno_t nInterior = 0;
      std::vector<lno_t> reorderToLocal;
      for(size_t i = 0;  i< nVtx; i++) reorderToLocal.push_back(i);
      if(verbose) std::cout<<comm->getRank()<<": Setting up local datastructures\n";
      //Set up a typical local mapping here.
      std::unordered_map<gno_t,lno_t> globalToLocal;
      std::vector<gno_t> ownedPlusGhosts;
      for (gno_t i = 0; i < vtxIDs.size(); i++){
        if(vtxIDs[i] < 0 && verbose) std::cout<<comm->getRank()<<": found a negative GID\n";
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
      }
      gno_t nghosts = 0;
      for (int i = 0; i < adjs.size(); i++){
        if(globalToLocal.count(adjs[i]) == 0){
          //new unique ghost found
          if(adjs[i] < 0 && verbose) std::cout<<comm->getRank()<<": found a negative adjacency\n";
          ownedPlusGhosts.push_back(adjs[i]);
          globalToLocal[adjs[i]] = vtxIDs.size() + nghosts;
          nghosts++;
            
        }
      }
      if(verbose) std::cout<<comm->getRank()<<": vec.max_size() = "<<finalAdjs_vec.max_size()<<", adjs.size() = "<<adjs.size()<<"\n";  
      finalAdjs_vec.resize(adjs.size()); 
      for(size_t i = 0; i < finalAdjs_vec.size();i++){
        finalAdjs_vec[i] = globalToLocal[adjs[i]];
      }
      for(int i = 0; i < offsets.size(); i++) finalOffset_vec.push_back(offsets[i]);
      finalGIDs = ownedPlusGhosts;
      

      if(verbose) std::cout<<comm->getRank()<<": creating Tpetra Maps\n";
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                             <Tpetra::global_size_t>::invalid();
      mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));

      dummy = Teuchos::OrdinalTraits <Tpetra::global_size_t>::invalid();
      mapWithCopies = rcp(new map_t(dummy, 
                                Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                0, comm)); 
                                      
      //create the FEMultiVector for the distributed communication.
      //We also use the views from this datastructure as arguments to
      //KokkosKernels coloring functions.
      if(verbose)std::cout<<comm->getRank()<<": creating FEMultiVector\n";
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, 
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, 
                                                    importer, 1, true));
      //Get color array to fill
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
        colors[i] = 0;
      } 
      
      //Create random numbers seeded on global IDs so that we don't
      //need to communicate for consistency. These numbers determine
      //which vertex gets recolored in the event of a conflict.
      //taken directly from the Zoltan coloring implementation 
      std::vector<int> rand(finalGIDs.size());
      for(size_t i = 0; i < finalGIDs.size(); i++){
        std::srand(finalGIDs[i]);
        rand[i] = std::rand();
      }

      std::vector<int> ghostOwners(finalGIDs.size() - nVtx);
      std::vector<gno_t> ghosts(finalGIDs.size() - nVtx);
      for(size_t i = nVtx; i < finalGIDs.size(); i++) ghosts[i-nVtx] = finalGIDs[i];
      ArrayView<int> owners = Teuchos::arrayViewFromVector(ghostOwners);
      ArrayView<const gno_t> ghostGIDs = Teuchos::arrayViewFromVector(ghosts);
      mapOwned->getRemoteIndexList(ghostGIDs, owners);
      
      for(size_t i = 0; i < finalOffset_vec.size()-1; i++){
        std::sort(finalAdjs_vec.begin()+finalOffset_vec[i],finalAdjs_vec.begin()+finalOffset_vec[i+1]);
      }
      
      ArrayView<const offset_t> finalOffsets = Teuchos::arrayViewFromVector(finalOffset_vec);
      ArrayView<const lno_t> finalAdjs = Teuchos::arrayViewFromVector(finalAdjs_vec);
      
      Teuchos::ArrayView<const lno_t> exportLIDs = importer->getExportLIDs();
      Teuchos::ArrayView<const int> exportPIDs = importer->getExportPIDs();
      
      for(int i = 0; i < exportLIDs.size(); i++){
        procs_to_send[exportLIDs[i]].push_back(exportPIDs[i]);  
      }

      // call coloring function
      hybridGMB(nVtx, nInterior, finalAdjs, finalOffsets,colors,femv,finalGIDs,rand,owners,mapWithCopies);
      
      //copy colors to the output array.
      for(int i = 0; i < colors.size(); i++){
        colors[reorderToLocal[i]] = femv->getData(0)[i];
      }
    
    }
     
    void hybridGMB(const size_t nVtx,lno_t nInterior, Teuchos::ArrayView<const lno_t> adjs, 
                   Teuchos::ArrayView<const offset_t> offsets, 
                   Teuchos::ArrayRCP<int> colors, Teuchos::RCP<femv_t> femv,
                   std::vector<gno_t> reorderGIDs,
                   std::vector<int> rand,
                   ArrayView<int> owners,
                   RCP<const map_t> mapOwnedPlusGhosts){
      if(verbose) std::cout<<comm->getRank()<<": inside coloring algorithm\n";
      
      //initialize timers
      double total_time = 0.0;
      double interior_time = 0.0;
      double comm_time = 0.0;
      double comp_time = 0.0;
      double recoloring_time = 0.0;
      double conflict_detection = 0.0;
      
      //determine max degree locally and globally
      offset_t local_max_degree = 0;
      offset_t global_max_degree = 0;
      for(size_t i = 0; i < nVtx; i++){
        offset_t curr_degree = offsets[i+1] - offsets[i];
        if(curr_degree > local_max_degree){
          local_max_degree = curr_degree;
        }
      }
      Teuchos::reduceAll<int, offset_t>(*comm,Teuchos::REDUCE_MAX,1, &local_max_degree, &global_max_degree);
      if(comm->getRank() == 0 && verbose) std::cout<<"Input has max degree "<<global_max_degree<<"\n";
      if(verbose)std::cout<<comm->getRank()<<": creating Kokkos Views\n"; 
     
      Kokkos::View<offset_t*, device_type> dist_degrees("Owned+Ghost degree view",rand.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_degrees_host = Kokkos::create_mirror(dist_degrees);
      //set degree counts for ghosts
      for(int i = 0; i < adjs.size(); i++){
        if((size_t)adjs[i] < nVtx) continue;
        dist_degrees_host(adjs[i])++;
      }
      //set degree counts for owned verts
      for(int i = 0; i < offsets.size()-1; i++){
        dist_degrees_host(i) = offsets[i+1] - offsets[i];
      }
      
      Kokkos::View<offset_t*, device_type> dist_offsets("Owned+Ghost Offset view", rand.size()+1);
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host = Kokkos::create_mirror(dist_offsets);

      //set offsets and total # of adjacencies
      dist_offsets_host(0) = 0;
      uint64_t total_adjs = 0;
      for(size_t i = 1; i < rand.size()+1; i++){
        dist_offsets_host(i) = dist_degrees_host(i-1) + dist_offsets_host(i-1);
        total_adjs+= dist_degrees_host(i-1);
      }

      Kokkos::View<lno_t*, device_type> dist_adjs("Owned+Ghost adjacency view", total_adjs);
      typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host = Kokkos::create_mirror(dist_adjs);
      //now, use the degree view as a counter
      for(size_t i = 0; i < rand.size(); i++){
        dist_degrees_host(i) = 0;
      }
      for(int i = 0; i < adjs.size(); i++) dist_adjs_host(i) = adjs[i];
      if(comm->getSize() > 1){
        for(size_t i = 0; i < nVtx; i++){
          for(size_t j = offsets[i]; j < offsets[i+1]; j++){
            //if the adjacency is a ghost
            if( (size_t)adjs[j] >= nVtx){
              //add the symmetric edge to its adjacency list (already accounted for by offsets)
              dist_adjs_host(dist_offsets_host(adjs[j]) + dist_degrees_host(adjs[j])) = i;
              dist_degrees_host(adjs[j])++;
            }
          }
      	}
      }
       
      if(verbose) std::cout<<comm->getRank()<<": copying host mirrors to device views\n";
      //copy the host data to device memory
      Kokkos::deep_copy(dist_degrees, dist_degrees_host); //may be unnecessary
      Kokkos::deep_copy(dist_offsets, dist_offsets_host);
      Kokkos::deep_copy(dist_adjs, dist_adjs_host);
      if(verbose) std::cout<<comm->getRank()<<": done copying to device\n";
      
      //counter in UVM memory for how many vertices need recoloring.
      Kokkos::View<gno_t[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;

      //device copy of the random tie-breakers
      Kokkos::View<int*,device_type> rand_dev("randVec",rand.size());
      typename Kokkos::View<int*, device_type>::HostMirror rand_host = Kokkos::create_mirror(rand_dev);
      for(size_t i = 0; i < rand.size(); i++){
        rand_host(i) = rand[i];
      }

      //device copy of global IDs
      Kokkos::View<gno_t*, device_type> gid_dev("GIDs",reorderGIDs.size());
      typename Kokkos::View<gno_t*,device_type>::HostMirror gid_host = Kokkos::create_mirror(gid_dev);
      for(size_t i = 0; i < reorderGIDs.size(); i++){
        gid_host(i) = reorderGIDs[i];
      }

      //copy host views to device memory
      Kokkos::deep_copy(rand_dev,rand_host);
      Kokkos::deep_copy(gid_dev, gid_host);

      if(verbose)std::cout<<comm->getRank()<<": done creating recoloring datastructures\n";
      //count boundary size to allocate list of vertices to recolor.
      offset_t boundary_size = 0;
      for(offset_t i = 0; i < nVtx; i++){
        for(offset_t j = offsets[i]; j < offsets[i+1]; j++){
          if((size_t)adjs[j] >= nVtx) {
            boundary_size++;
            break;
          }
        }
      }
      if(verbose)std::cout<<comm->getRank()<<": creating send views\n";
      
      //list of vertices to send to remote processes
      Kokkos::View<lno_t*, device_type> verts_to_send_view("verts to send",boundary_size);
      Kokkos::parallel_for(boundary_size, KOKKOS_LAMBDA(const int& i){
        verts_to_send_view(i) = -1;
      });
      //atomic copy
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_atomic = verts_to_send_view;
      
      //size information for the list of vertices to send. Also includes an atomic copy
      Kokkos::View<size_t[1], device_type> verts_to_send_size("verts to send size");
      verts_to_send_size(0) = 0;
      Kokkos::View<size_t[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_size_atomic = verts_to_send_size;

      if(verbose)std::cout<<comm->getRank()<<": Done creating send views, initializing...\n";
      if(verbose)std::cout<<comm->getRank()<<": boundary_size = "<<boundary_size<<" verts_to_send_size_atomic(0) = "<<verts_to_send_size_atomic(0)<<"\n";
      //initially the verts to send include all boundary vertices.
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int&i){
        for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
	  if((size_t)dist_adjs(j) >= nVtx){
	    verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	    break;
	  }
	} 
      });
      Kokkos::fence();
      
      
      //used to replace the incorrectly recolored ghosts
      Kokkos::View<int*, device_type> ghost_colors("ghost color backups", rand.size()-nVtx);
      if(verbose)std::cout<<comm->getRank()<<": Done initializing\n";
      gno_t sentPerRound[100];
      gno_t recvPerRound[100];

      if(verbose) std::cout<<comm->getRank()<<": Coloring interior\n";
      //initialize interior and total timers, barrier to prevent any imbalance from setup.
      //Only use a barrier if timing is happening.
      if(verbose) comm->barrier();
      interior_time = timer();
      total_time = timer();
      //call the KokkosKernels coloring function with the Tpetra default spaces.
      bool use_vbbit = (global_max_degree < 6000);
      this->colorInterior<execution_space, memory_space,memory_space>
                 (nVtx, dist_adjs, dist_offsets, femv,dist_adjs,0,use_vbbit);
      if(verbose){
        interior_time = timer() - interior_time;
        comp_time = interior_time;
        std::cout<<comm->getRank()<<": Going to recolor\n";
      }

      //if there is more than a single process, check distributed conflicts and recolor
      if(comm->getSize() > 1){
	//get the colors from the femv
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
	
        if(verbose)std::cout<<comm->getRank()<<": going to communicate\n";

	//communicate
        gno_t recv, sent;
        comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,nVtx,verts_to_send_view,verts_to_send_size,femv_colors,recv,sent);
        sentPerRound[0] = sent;
        recvPerRound[0] = recv;
        if(verbose) std::cout<<comm->getRank()<<": done communicating\n";
	verts_to_send_size(0) = 0;

        //set the old ghost colors
        Kokkos::parallel_for(rand.size()-nVtx,KOKKOS_LAMBDA(const int& i){
          ghost_colors(i) = femv_colors(i+nVtx);
        });
	Kokkos::fence();
	//detect conflicts on the device, uncolor conflicts consistently.
        double temp = timer();
        Kokkos::parallel_reduce(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i,gno_t& recoloring_size){
          lno_t lid = i+nVtx;
          int currColor = femv_colors(lid);
          for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            if(currColor == nborColor ){
              if(rand_dev(lid) > rand_dev(dist_adjs(j))){
                femv_colors(lid) = 0;
                recoloring_size++;
                break;
              }else if (rand_dev(dist_adjs(j)) > rand_dev(lid)){
                femv_colors(dist_adjs(j)) = 0;
                recoloring_size++;
              } else {
                if (gid_dev(lid) >= gid_dev(dist_adjs(j))){
                  femv_colors(lid) = 0;
                  recoloring_size++;
                  break;
                } else {
                  femv_colors(dist_adjs(j)) = 0;
                  recoloring_size++;
                }
              }
            }
          }
        },recoloringSize(0));
        Kokkos::fence();
	//build the verts to send list.
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
          if(femv_colors(i) == 0){
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	  }
        });
        //ensure the parallel_for has completed before continuing.
        Kokkos::fence();
        conflict_detection += timer() - temp;
        //total_time += conflict_detection;
        comp_time += conflict_detection;
      }
      //end the initial recoloring
      if(verbose)std::cout<<comm->getRank()<<": done initial recoloring, begin recoloring loop\n";
      double totalPerRound[100];
      double commPerRound[100];
      double compPerRound[100];
      double recoloringPerRound[100];
      double conflictDetectionPerRound[100];
      double serialRecoloringPerRound[100];
      int vertsPerRound[100];
      bool done = false; //We're only done when all processors are done
      if(comm->getSize() == 1) done = true;
      totalPerRound[0] = interior_time + comm_time + conflict_detection;
      recoloringPerRound[0] = 0;
      commPerRound[0] = comm_time;
      compPerRound[0] = interior_time + conflict_detection;
      conflictDetectionPerRound[0] = conflict_detection;
      recoloringPerRound[0] = 0;
      vertsPerRound[0] = 0;
      int distributedRounds = 1; //this is the same across all processors
      Kokkos::View<lno_t*, device_type> verts_to_recolor("verts_to_recolor", boundary_size);
      //while the queue is not empty
      while(recoloringSize(0) > 0 || !done){
        //get a subview of the colors:
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        //color everything in the recoloring queue, put everything on conflict queue
        if(distributedRounds < 100) {
          vertsPerRound[distributedRounds] = recoloringSize(0);
        }
        if(verbose) std::cout<<comm->getRank()<<": starting to recolor\n";
	
	Kokkos::deep_copy(verts_to_recolor, verts_to_send_view);

        double recolor_temp = timer();
        //use KokkosKernels to recolor the conflicting vertices.  
	if(verts_to_send_size(0) > 0){
            this->colorInterior<execution_space,
                                memory_space,
                                memory_space>
                               (femv_colors.size(),dist_adjs,dist_offsets,femv,verts_to_recolor,verts_to_send_size(0),use_vbbit);
	}
        recoloringPerRound[distributedRounds] = timer() - recolor_temp;
        recoloring_time += recoloringPerRound[distributedRounds];
        comp_time += recoloringPerRound[distributedRounds];
        compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	
        
        if(verbose) std::cout<<comm->getRank()<<": done recoloring\n";
        recoloringSize(0) = 0;

        Kokkos::parallel_for(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i){
          femv_colors(i+nVtx) = ghost_colors(i);
        });
        Kokkos::fence();
        //communicate
        typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send_host = Kokkos::create_mirror(verts_to_send_view);
        Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
	gno_t sent,recv;
        commPerRound[distributedRounds] = doOwnedToGhosts(mapOwnedPlusGhosts,nVtx,verts_to_send_host, verts_to_send_size,femv_colors,recv,sent); 	
        if(verbose){
          std::cout<<comm->getRank()<<": total sent in round "<<distributedRounds<<" = "<<sent<<"\n";
          std::cout<<comm->getRank()<<": total recv in round "<<distributedRounds<<" = "<<recv<<"\n";
	}
        sentPerRound[distributedRounds] = sent;
        recvPerRound[distributedRounds] = recv;
        comm_time += commPerRound[distributedRounds];
        totalPerRound[distributedRounds] += commPerRound[distributedRounds];

        //detect conflicts in parallel. For a detected conflict,
        //reset the vertex-to-be-recolored's color to 0, in order to
        //allow KokkosKernels to recolor correctly.
        Kokkos::parallel_for(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i){
          ghost_colors(i) = femv_colors(i+nVtx);
        });
        Kokkos::fence();
	verts_to_send_size(0) = 0;
        double detection_temp = timer();
        Kokkos::parallel_reduce(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i,gno_t& recoloring_size){
          lno_t lid = i+nVtx;
          int currColor = femv_colors(lid);
          for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            if(currColor == nborColor ){
              if(rand_dev(lid) > rand_dev(dist_adjs(j))){
                femv_colors(lid) = 0;
                recoloring_size++;
                break;
              }else if (rand_dev(dist_adjs(j)) > rand_dev(lid)){
                femv_colors(dist_adjs(j)) = 0;
                recoloring_size++;
              } else {
                if (gid_dev(lid) >= gid_dev(dist_adjs(j))){
                  femv_colors(lid) = 0;
                  recoloring_size++;
                  break;
                } else {
                  femv_colors(dist_adjs(j)) = 0;
                  recoloring_size++;
                }
              }
            }
          }
        },recoloringSize(0));
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
          if(femv_colors(i) == 0){
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
          }
        });
        
        //For Cuda, this fence is necessary to ensure the Kokkos::parallel_for is finished
        //before continuing with the coloring. 
        Kokkos::fence();
        conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
        conflict_detection += conflictDetectionPerRound[distributedRounds];
        compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        comp_time += conflictDetectionPerRound[distributedRounds];
        
        //do a reduction to determine if we're done
        int globalDone = 0;
        int localDone = recoloringSize(0);
        Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
        //We're only allowed to stop once everyone has no work to do.
        //collectives will hang if one process exits. 
        distributedRounds++;
        done = !globalDone;
      }
      total_time = timer() - total_time;
      
      //only take time to compute statistics if the user wants them 
      if(verbose){
        std::cout<<comm->getRank()<<": done recoloring loop, computing statistics\n";
        int localBoundaryVertices = 0;
        for(size_t i = 0; i < nVtx; i++){
          for(offset_t j = offsets[i]; j < offsets[i+1]; j++){
            if((size_t)adjs[j] >= nVtx){
              localBoundaryVertices++;
              break;
            }
          }
        }
        //print how many rounds of speculating/correcting happened (this should be the same for all ranks):
        if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n", distributedRounds);
        int totalVertsPerRound[100];
        int totalBoundarySize = 0;
        double finalTotalPerRound[100];
        double maxRecoloringPerRound[100];
        double finalSerialRecoloringPerRound[100];
        double minRecoloringPerRound[100];
        double finalCommPerRound[100];
        double finalCompPerRound[100];
        double finalConflictDetectionPerRound[100];
        gno_t finalRecvPerRound[100];
        gno_t finalSentPerRound[100];
        for(int i = 0; i < 100; i++) {
          totalVertsPerRound[i] = 0;
          finalTotalPerRound[i] = 0.0;
          maxRecoloringPerRound[i] = 0.0;
          minRecoloringPerRound[i] = 0.0;
          finalCommPerRound[i] = 0.0;
          finalCompPerRound[i] = 0.0;
          finalConflictDetectionPerRound[i] = 0.0;
          finalSentPerRound[i] = 0;
          finalRecvPerRound[i] = 0;
        }
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,1, &localBoundaryVertices,&totalBoundarySize);
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,100,vertsPerRound,totalVertsPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,totalPerRound,finalTotalPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,recoloringPerRound,maxRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,100,recoloringPerRound,minRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,serialRecoloringPerRound,finalSerialRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,commPerRound,finalCommPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,compPerRound,finalCompPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,conflictDetectionPerRound, finalConflictDetectionPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,recvPerRound, finalRecvPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,sentPerRound, finalSentPerRound);
        
        printf("Rank %d: boundary size: %d\n",comm->getRank(),localBoundaryVertices);
        for(int i = 0; i < std::min(distributedRounds,100); i++){
          printf("Rank %d: recolor %d vertices in round %d\n",comm->getRank(),vertsPerRound[i],i);
          if(comm->getRank()==0) printf("recolored %d vertices in round %d\n",totalVertsPerRound[i],i);
          if(comm->getRank()==0) printf("total time in round %d: %f\n",i,finalTotalPerRound[i]);
          if(comm->getRank()==0) printf("recoloring time in round %d: %f\n",i,maxRecoloringPerRound[i]);
          if(comm->getRank()==0) printf("serial recoloring time in round %d: %f\n",i,finalSerialRecoloringPerRound[i]);
          if(comm->getRank()==0) printf("min recoloring time in round %d: %f\n",i,minRecoloringPerRound[i]);
          if(comm->getRank()==0) printf("conflict detection time in round %d: %f\n",i,finalConflictDetectionPerRound[i]);
          if(comm->getRank()==0) printf("comm time in round %d: %f\n",i,finalCommPerRound[i]);
          if(comm->getRank()==0) printf("total sent in round %d: %lld\n",i,finalSentPerRound[i]);
          if(comm->getRank()==0) printf("total recv in round %d: %lld\n",i,finalRecvPerRound[i]);
          if(comm->getRank()==0) printf("comp time in round %d: %f\n",i,finalCompPerRound[i]);
        }
        
        double global_total_time = 0.0;
        double global_recoloring_time=0.0;
        double global_min_recoloring_time=0.0;
        double global_conflict_detection=0.0;
        double global_comm_time=0.0;
        double global_comp_time=0.0;
        double global_interior_time = 0.0;
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
          printf("Boundary size: %d\n",totalBoundarySize);
          printf("Total Time: %f\n",global_total_time);
          printf("Interior Time: %f\n",global_interior_time);
          printf("Recoloring Time: %f\n",global_recoloring_time);
          printf("Min Recoloring Time: %f\n",global_min_recoloring_time);
          printf("Conflict Detection Time: %f\n",global_conflict_detection);
          printf("Comm Time: %f\n",global_comm_time);
          printf("Comp Time: %f\n",global_comp_time);
        }
        std::cout<<comm->getRank()<<": exiting coloring\n";
      }
    }
};

template <typename Adapter>
void AlgHybridGMB<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
    
  this->env->debug(DETAILED_STATUS, "   building graph model");
  std::cout<<comm->getRank()<<": starting to construct graph model\n";
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
  std::cout<<comm->getRank()<<": done constructing graph model\n";
  this->env->debug(DETAILED_STATUS, "   graph model built");
}


}//end namespace Zoltan2
#endif
