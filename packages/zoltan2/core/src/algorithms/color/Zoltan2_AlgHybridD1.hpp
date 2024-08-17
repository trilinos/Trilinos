// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGHYBRIDD1_HPP_
#define _ZOLTAN2_ALGHYBRIDD1_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <queue>
#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <algorithm>

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
#include <stdlib.h>

//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridD1.hpp
//! \brief A hybrid MPI+Kokkos version of the framework proposed by Gebremedhin and Manne

namespace Zoltan2 {

template <typename Adapter>
class AlgDistance1 : public Algorithm<Adapter>
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
    using device_type = typename femv_t::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using host_exec = typename femv_t::host_view_type::device_type::execution_space;
    using host_mem = typename femv_t::host_view_type::device_type::memory_space;
    double timer() {
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
    }
    
  private:

    void buildModel(modelFlag_t &flags); 

    //function to invoke KokkosKernels distance-1 coloring
    //
    //OUTPUT ARG:
    //  
    //  femv: FEMultiVector containing dual-view of colors.
    //        After this call each vertex will have a locally-valid color.
    //
    //INPUT ARGS:
    //  
    //  nVtx: number of local owned vertices
    //
    //  adjs_view: CSR adjacencies for the local graph
    //
    //  offset_view: CSR offsets for indexing adjs_view
    //
    //  vertex_list: list of local IDs of vertices that 
    //               need to be recolored
    //
    //  vertex_list_size: 0 means no vertex list given,
    //                    otherwise it simply is the size
    //                    of the vertex_list
    //
    //  recolor: switches between VB_BIT and EB KokkosKernels
    //           algorithms
    //
    template <class ExecutionSpace, typename MemorySpace>
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
      
      //pick which KokkosKernels algorithm to use.
      //VBBIT is more efficient for inputs with max degree < ~6000
      //EB is more efficient for inputs with max degree > ~6000
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
      auto femvColors = femv->template getLocalView<Kokkos::Device<ExecutionSpace,MemorySpace> >(Tpetra::Access::ReadWrite);
      auto  sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      kh.get_graph_coloring_handle()->set_tictoc(verbose);
      KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t>
                                                     (&kh, nVtx, nVtx, offset_view, adjs_view);
      //numColors = kh.get_graph_coloring_handle()->get_num_colors();

      if(verbose){
        std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
      }
    }
    
  public:
    //contains all distance-1 conflict detection logic
    //
    //OUTPUT ARGS:
    //  
    //  femv_colors: This function uncolors vertices that have
    //               distance-1 conflicts, so these colors will
    //               change if there are any conflicts present
    //
    //INPUT ARGS:
    //
    //  n_local: number of locally owned vertices
    //
    //  n_ghosts: number of ghosts on this process
    //
    //  dist_offsets: CSR offsets of the local graph
    //
    //  dist_adjs: CSR adjacencies of the local graph
    //
    //  verts_to_send_view: Used to construct a list of verts to send on device.
    //
    //  verts_to_send_size_atomic: atomic version of the verts_to_send_size view
    //                             Used to construct a list on device,
    //                             the atomic is necessary for correctness.
    //
    //  recoloringSize: holds the total amount of recoloring that will be done
    //                  locally. Does not need to be atomic.
    //
    //  rand: view that holds the tie-breaking random numbers indexed by LID.
    //
    //  gid: view that holds GIDs, for tie breaking in the case that rand
    //       numbers are the same for two vertices. 
    //
    //  ghost_degrees: view that holds degrees only for ghost vertices.
    //                 LID n_local corresponds to ghost_degrees(0).
    //
    //  recolor_degrees: if true, we factor degrees into the conflict detection
    //                   if false, we resolve using only consistent random numbers.
    //
    template <class ExecutionSpace, typename MemorySpace>
    void detectConflicts(const size_t n_local, const size_t n_ghosts,
		         Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> dist_offsets,
			 Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> dist_adjs,
			 Kokkos::View<int*, Kokkos::Device<ExecutionSpace, MemorySpace>> femv_colors,
			 Kokkos::View<lno_t*,
			              Kokkos::Device<ExecutionSpace, MemorySpace>> verts_to_send_view,
		         Kokkos::View<size_t*,
			              Kokkos::Device<ExecutionSpace, MemorySpace>,
				      Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic,
			 Kokkos::View<gno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> recoloringSize,
			 Kokkos::View<int*,
			              Kokkos::Device<ExecutionSpace, MemorySpace>> rand,
			 Kokkos::View<gno_t*,
			              Kokkos::Device<ExecutionSpace, MemorySpace>> gid,
                         Kokkos::View<gno_t*,
			              Kokkos::Device<ExecutionSpace, MemorySpace>> ghost_degrees,
			 bool recolor_degrees){
      gno_t local_recoloring_size;
      Kokkos::parallel_reduce("Conflict Detection",
		              Kokkos::RangePolicy<ExecutionSpace>(0,n_ghosts), 
		              KOKKOS_LAMBDA(const int& i,
				            gno_t& recoloring_size){
        lno_t lid = i+n_local;
        int currColor = femv_colors(lid);
	int currDegree = ghost_degrees(i);
        for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
          int nborColor = femv_colors(dist_adjs(j));
	  int nborDegree = 0;
	  if((size_t)dist_adjs(j) < n_local) nborDegree = dist_offsets(dist_adjs(j)+1) - dist_offsets(dist_adjs(j));
	  else nborDegree = ghost_degrees(dist_adjs(j) - n_local); 
          if(currColor == nborColor ){
	    if(currDegree < nborDegree && recolor_degrees){
	      femv_colors(lid) = 0;
	      recoloring_size++;
	    }else if (nborDegree < currDegree && recolor_degrees){
	      femv_colors(dist_adjs(j)) = 0;
	      recoloring_size++;
            }else if(rand(lid) > rand(dist_adjs(j))){
              femv_colors(lid) = 0;
              recoloring_size++;
              break;
            }else if (rand(dist_adjs(j)) > rand(lid)){
              femv_colors(dist_adjs(j)) = 0;
              recoloring_size++;
            } else {
              if (gid(lid) >= gid(dist_adjs(j))){
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
      },local_recoloring_size);
      Kokkos::deep_copy(recoloringSize, local_recoloring_size);
      Kokkos::fence();
      Kokkos::parallel_for("Rebuild verts_to_send_view",
		           Kokkos::RangePolicy<ExecutionSpace>(0,n_local), 
			   KOKKOS_LAMBDA(const int& i){
        if(femv_colors(i) == 0){
          verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
        }
      });
      Kokkos::fence();
    }

  private:
    //Communicates owned vertex colors to ghost copies.
    //
    //RETURN VALUE:
    //
    //  returns the time it took for the communication call to complete
    //
    //OUTPUT ARGS:
    //  
    //  femv: FEMultivector that holds the dual-view of the colors
    //        After this call, ghost colors will be up-to-date.
    //
    //  recv: returns the size of the recv buffer
    //
    //  send: returns the size of the send buffer
    //  
    //INPUT ARGS:
    //
    //  mapOwnedPlusGhosts: a Tpetra map that translates between
    //                      LID and GID for any vertex on this process.
    //
    //  nVtx: the number of locally owned vertices.
    //
    //  verts_to_send: hostmirror of verts to send. This function sends
    //                 all vertices in this list to their ghost copies.
    //
    //  verts_to_send_size: hostmirror of verts_to_send_size, holds the
    //                      size of verts_to_send.
    //
    //  procs_to_send: map that translates LID into a list of processes that
    //                 have a ghost copy of the vertex. 
    //
    double doOwnedToGhosts(RCP<const map_t> mapOwnedPlusGhosts,
                         size_t nVtx,
			 typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send,
			 typename Kokkos::View<size_t*, device_type>::HostMirror& verts_to_send_size,
                         Teuchos::RCP<femv_t> femv,
			 const std::unordered_map<lno_t, std::vector<int>>& procs_to_send,
                         gno_t& recv, gno_t& send){ 
      auto femvColors = femv->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto colors = subview(femvColors, Kokkos::ALL, 0);
      int nprocs = comm->getSize();
      
      std::vector<int> sendcnts(comm->getSize(), 0);
      std::vector<gno_t> sdispls(comm->getSize()+1, 0);
      for(size_t i = 0; i < verts_to_send_size(0); i++){
	for(size_t j = 0; j < procs_to_send.at(verts_to_send(i)).size(); j++){
	  sendcnts[procs_to_send.at(verts_to_send(i))[j]] += 2;
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
	std::vector<int> procs = procs_to_send.at(verts_to_send(i));
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

      //if we're reporting timings, remove the computation imbalance from the comm timer.
      if(timing) comm->barrier();
      double comm_total = 0.0;
      double comm_temp = timer();

      AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcnts_view, recvbuf, recvcnts_view);
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
    bool timing;
  public:
    //constructor for the  hybrid distributed distance-1 algorithm
    AlgDistance1(
      const RCP<const base_adapter_t> &adapter_, 
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_) {
      verbose = pl->get<bool>("verbose",false);
      timing = pl->get<bool>("timing", false);
      if(verbose) std::cout<<comm->getRank()<<": inside coloring constructor\n";
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
      
      //find out who owns the ghosts on this process
      std::vector<int> ghostOwners(finalGIDs.size() - nVtx);
      std::vector<gno_t> ghosts(finalGIDs.size() - nVtx);
      for(size_t i = nVtx; i < finalGIDs.size(); i++) ghosts[i-nVtx] = finalGIDs[i];
      ArrayView<int> owners = Teuchos::arrayViewFromVector(ghostOwners);
      ArrayView<const gno_t> ghostGIDs = Teuchos::arrayViewFromVector(ghosts);
      mapOwned->getRemoteIndexList(ghostGIDs, owners);
      
     //create const views of the CSR representation of the local graph 
      ArrayView<const offset_t> finalOffsets = Teuchos::arrayViewFromVector(finalOffset_vec);
      ArrayView<const lno_t> finalAdjs = Teuchos::arrayViewFromVector(finalAdjs_vec);
      ArrayView<const int> rand_view = Teuchos::arrayViewFromVector(rand);
      ArrayView<const gno_t> gid_view = Teuchos::arrayViewFromVector(finalGIDs);
      
      //find out which remote processes have a ghost copy of a local owned vertex.
      Teuchos::ArrayView<const lno_t> exportLIDs = importer->getExportLIDs();
      Teuchos::ArrayView<const int> exportPIDs = importer->getExportPIDs();
      
      //create a quick lookup from LID -> remote processes with a ghost copy
      std::unordered_map<lno_t, std::vector<int>> procs_to_send; 
      for(int i = 0; i < exportLIDs.size(); i++){
        procs_to_send[exportLIDs[i]].push_back(exportPIDs[i]);  
      }

      // call coloring function
      hybridGMB(nVtx, finalAdjs, finalOffsets,femv,gid_view,rand_view,owners,mapWithCopies, procs_to_send);
      
      //copy colors to the output array.
      auto femvdata = femv->getData(0);
      for(int i = 0; i < colors.size(); i++){
        colors[reorderToLocal[i]] = femvdata[i];
      }
    
    }
    
    //This function contains all coloring logic.
    //
    //OUTPUT ARGS:
    //
    //  femv: FEMultiVector with a dual-view of the vertex colors.
    //        The colors will be a valid distance-1 coloring after this
    //        function returns.
    //
    //INPUT ARGS: 
    //  
    //  nVtx: number of locally owned vertices
    //
    //  adjs: CSR adjacencies for the local graph
    //
    //  offsets: CSR offsets into adjs for the local graph
    //
    //  gids: global IDs for every vertex on this process.
    //        indexed by local ID
    //
    //  rand: random number generated for every vertex on 
    //        this process. Indexed by local ID
    //
    //  owners: owners of each ghost vertex.
    //          owners[i] = the owning process for vertex with GID gids[i+nVtx]
    //
    //  mapOwnedPlusGhosts: Tpetra map that translates between LIDs and GIDs
    //
    //  procs_to_send: maps from LIDs to Process IDs.
    //                 procs_to_send[LID] gives a list of processes that have
    //                 a ghost copy of LID. 
    //
    void hybridGMB(const size_t nVtx,
		   const Teuchos::ArrayView<const lno_t>& adjs, 
                   const Teuchos::ArrayView<const offset_t>& offsets, 
                   const Teuchos::RCP<femv_t>& femv,
                   const Teuchos::ArrayView<const gno_t>& gids,
                   const Teuchos::ArrayView<const int>& rand,
                   const Teuchos::ArrayView<const int>& owners,
                   RCP<const map_t> mapOwnedPlusGhosts,
		   const std::unordered_map<lno_t, std::vector<int>>& procs_to_send){
      if(verbose) std::cout<<comm->getRank()<<": inside coloring algorithm\n";
      
      //initialize timers
      double total_time = 0.0;
      double interior_time = 0.0;
      double comm_time = 0.0;
      double comp_time = 0.0;
      double recoloring_time = 0.0;
      double conflict_detection = 0.0;
      
      const int numStatisticRecordingRounds = 100;

      //Put together local and remote degrees
      //1. Communicate remote GIDs to owning processes
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
	 deg_sendbuf[idx] = mapOwnedPlusGhosts->getGlobalElement(i+nVtx);
      }
      Teuchos::ArrayView<int> deg_send_cnts_view = Teuchos::arrayViewFromVector(deg_send_cnts);
      Teuchos::ArrayView<gno_t> deg_sendbuf_view = Teuchos::arrayViewFromVector(deg_sendbuf);
      Teuchos::ArrayRCP<gno_t> deg_recvbuf;
      std::vector<int> deg_recvcnts(comm->getSize(),0);
      Teuchos::ArrayView<int> deg_recvcnts_view = Teuchos::arrayViewFromVector(deg_recvcnts);
      AlltoAllv<gno_t>(*comm, *env, deg_sendbuf_view, deg_send_cnts_view, deg_recvbuf, deg_recvcnts_view);

      //2. replace GID with local degree
      for(int i = 0; i < deg_recvbuf.size(); i++){
        lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_recvbuf[i]);
	deg_recvbuf[i] = offsets[lid+1] - offsets[lid];
      }
      //3. send modified buffer back
      ArrayRCP<gno_t> ghost_degrees;
      AlltoAllv<gno_t>(*comm, *env, deg_recvbuf(), deg_recvcnts_view, ghost_degrees, deg_send_cnts_view);
      //determine max degree locally and globally
      Kokkos::View<gno_t*, device_type> ghost_degrees_dev("ghost degree view",ghost_degrees.size());
      typename Kokkos::View<gno_t*, device_type>::HostMirror ghost_degrees_host = Kokkos::create_mirror(ghost_degrees_dev);
      for(int i = 0; i < ghost_degrees.size(); i++){
	lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_sendbuf[i]);
        ghost_degrees_host(lid-nVtx) = ghost_degrees[i];
      }
      Kokkos::deep_copy(ghost_degrees_dev, ghost_degrees_host);
 
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
      for(Teuchos_Ordinal i = 1; i < rand.size()+1; i++){
        dist_offsets_host(i) = dist_degrees_host(i-1) + dist_offsets_host(i-1);
        total_adjs+= dist_degrees_host(i-1);
      }

      Kokkos::View<lno_t*, device_type> dist_adjs("Owned+Ghost adjacency view", total_adjs);
      typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host = Kokkos::create_mirror(dist_adjs);
      //now, use the degree view as a counter
      for(Teuchos_Ordinal i = 0; i < rand.size(); i++){
        dist_degrees_host(i) = 0;
      }
      for(int i = 0; i < adjs.size(); i++) dist_adjs_host(i) = adjs[i];
      if(comm->getSize() > 1){
        for(size_t i = 0; i < nVtx; i++){
          for(offset_t j = offsets[i]; j < offsets[i+1]; j++){
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
      Kokkos::View<gno_t*, device_type> recoloringSize("Recoloring Queue Size",1);
      typename Kokkos::View<gno_t*, device_type>::HostMirror recoloringSize_host = Kokkos::create_mirror(recoloringSize);
      recoloringSize_host(0) = 0;
      Kokkos::deep_copy(recoloringSize, recoloringSize_host);

      //device copy of the random tie-breakers
      Kokkos::View<int*,device_type> rand_dev("randVec",rand.size());
      typename Kokkos::View<int*, device_type>::HostMirror rand_host = Kokkos::create_mirror(rand_dev);
      for(Teuchos_Ordinal i = 0; i < rand.size(); i++){
        rand_host(i) = rand[i];
      }

      //device copy of global IDs
      Kokkos::View<gno_t*, device_type> gid_dev("GIDs",gids.size());
      typename Kokkos::View<gno_t*,device_type>::HostMirror gid_host = Kokkos::create_mirror(gid_dev);
      for(Teuchos_Ordinal i = 0; i < gids.size(); i++){
        gid_host(i) = gids[i];
      }

      //copy host views to device memory
      Kokkos::deep_copy(rand_dev,rand_host);
      Kokkos::deep_copy(gid_dev, gid_host);

      if(verbose)std::cout<<comm->getRank()<<": done creating recoloring datastructures\n";
      //count boundary size to allocate list of vertices to recolor.
      offset_t boundary_size = 0;
      for(size_t i = 0; i < nVtx; i++){
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
      Kokkos::parallel_for("init verts_to_send_view", 
        Kokkos::RangePolicy<execution_space, int>(0,boundary_size),
        KOKKOS_LAMBDA(const int& i){
          verts_to_send_view(i) = -1;
        });
      
      //size information for the list of vertices to send. Also includes an atomic copy
      Kokkos::View<size_t*, device_type> verts_to_send_size("verts to send size",1);
      Kokkos::View<size_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_size_atomic = verts_to_send_size;
      typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send_host = create_mirror(verts_to_send_view);
      typename Kokkos::View<size_t*,device_type>::HostMirror verts_to_send_size_host = create_mirror(verts_to_send_size);
      //initialize the device view with a value of zero
      verts_to_send_size_host(0) = 0;
      deep_copy(verts_to_send_size, verts_to_send_size_host);

      if(verbose)std::cout<<comm->getRank()<<": Done creating send views, initializing...\n";
      if(verbose)std::cout<<comm->getRank()<<": boundary_size = "<<boundary_size<<" verts_to_send_size_atomic(0) = "<<verts_to_send_size_atomic(0)<<"\n";
      //initially the verts to send include all boundary vertices.
      Kokkos::parallel_for("Initialize verts_to_send",
        Kokkos::RangePolicy<execution_space, int>(0,nVtx),
        KOKKOS_LAMBDA(const int&i){
          for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
	    if((size_t)dist_adjs(j) >= nVtx){
	      verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
	      break;
	    }
	  } 
        });
      Kokkos::fence();
      
      
      //used to replace the incorrectly recolored ghosts
      Kokkos::View<int*, device_type> ghost_colors("ghost color backups", rand.size()-nVtx);
      if(verbose)std::cout<<comm->getRank()<<": Done initializing\n";
      gno_t sentPerRound[numStatisticRecordingRounds];
      gno_t recvPerRound[numStatisticRecordingRounds];

      if(verbose) std::cout<<comm->getRank()<<": Coloring interior\n";
      //initialize interior and total timers, barrier to prevent any imbalance from setup.
      //Only use a barrier if timing is happening.
      if(timing) comm->barrier();
      interior_time = timer();
      total_time = timer();
      //call the KokkosKernels coloring function with the Tpetra default spaces.
      bool use_vbbit = (global_max_degree < 6000);
      this->colorInterior<execution_space,memory_space>
                 (nVtx, dist_adjs, dist_offsets, femv,dist_adjs,0,use_vbbit);
      if(timing){
        interior_time = timer() - interior_time;
        comp_time = interior_time;
      }
      if(verbose) std::cout<<comm->getRank()<<": Going to recolor\n";
      bool recolor_degrees = this->pl->template get<bool>("recolor_degrees", true);

      //if there is more than a single process, check distributed conflicts and recolor
      if(comm->getSize() > 1){
	
        if(verbose)std::cout<<comm->getRank()<<": going to communicate\n";

	//communicate
        Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
	Kokkos::deep_copy(verts_to_send_size_host, verts_to_send_size);
        gno_t recv, sent;
        comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,
			            nVtx,
				    verts_to_send_host,
				    verts_to_send_size_host,
				    femv,
				    procs_to_send,
				    recv,
				    sent);
        sentPerRound[0] = sent;
        recvPerRound[0] = recv;
        if(verbose) std::cout<<comm->getRank()<<": done communicating\n";
	verts_to_send_size_host(0) = 0;
        deep_copy(verts_to_send_size, verts_to_send_size_host);
        //set the old ghost colors
	//get the colors from the femv
        Kokkos::View<int**, Kokkos::LayoutLeft, device_type> femvColors =
	  femv->template getLocalView<device_type>(Tpetra::Access::ReadWrite); // Partial write
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        Kokkos::parallel_for("get colors from femv",
          Kokkos::RangePolicy<execution_space, int>(0,rand.size()-nVtx),
          KOKKOS_LAMBDA(const int& i){
            ghost_colors(i) = femv_colors(i+nVtx);
          });
	Kokkos::fence();
	//detect conflicts on the device, uncolor conflicts consistently.
        double temp = timer();
	detectConflicts<execution_space, memory_space>(nVtx,
			                               rand.size()-nVtx,
						       dist_offsets,
						       dist_adjs,
						       femv_colors,
						       verts_to_send_view,
						       verts_to_send_size_atomic,
						       recoloringSize,
						       rand_dev,
						       gid_dev,
						       ghost_degrees_dev,
						       recolor_degrees);
        deep_copy(recoloringSize_host, recoloringSize);
        conflict_detection += timer() - temp;
        comp_time += conflict_detection;
      }
      //end the initial recoloring
      if(verbose)std::cout<<comm->getRank()<<": done initial recoloring, begin recoloring loop\n";
      double totalPerRound[numStatisticRecordingRounds];
      double commPerRound[numStatisticRecordingRounds];
      double compPerRound[numStatisticRecordingRounds];
      double recoloringPerRound[numStatisticRecordingRounds];
      double conflictDetectionPerRound[numStatisticRecordingRounds];
      double serialRecoloringPerRound[numStatisticRecordingRounds];
      int vertsPerRound[numStatisticRecordingRounds];
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
      int serial_threshold = this->pl->template get<int>("serial_threshold",0);
      
      Kokkos::View<lno_t*, device_type> verts_to_recolor("verts_to_recolor", boundary_size);
      typename Kokkos::View<int*, device_type>::HostMirror ghost_colors_host;
      //while the queue is not empty
      while(recoloringSize_host(0) > 0 || !done){
	if(recoloringSize_host(0) < serial_threshold) break;
        //get a subview of the colors:
        auto femvColors = femv->getLocalViewDevice(Tpetra::Access::ReadWrite);
        auto femv_colors = subview(femvColors, Kokkos::ALL, 0);
        //color everything in the recoloring queue, put everything on conflict queue
        if(distributedRounds < numStatisticRecordingRounds) {
          vertsPerRound[distributedRounds] = recoloringSize_host(0);
        }
	
	//copying the send view to the recolor view is necessary because
	//KokkosKernels can change the view passed in, and we need the send view
	//intact for communication.
	Kokkos::deep_copy(verts_to_recolor, verts_to_send_view); 

        double recolor_temp = timer();
        //use KokkosKernels to recolor the conflicting vertices.  
        deep_copy(verts_to_send_size_host, verts_to_send_size);
	if(verts_to_send_size_host(0) > 0){
            this->colorInterior<execution_space,
                                memory_space>(femv_colors.size(),
					      dist_adjs,dist_offsets,
					      femv,
					      verts_to_recolor,
					      verts_to_send_size_host(0),
					      use_vbbit);
	}

	if(distributedRounds < numStatisticRecordingRounds){
          recoloringPerRound[distributedRounds] = timer() - recolor_temp;
          recoloring_time += recoloringPerRound[distributedRounds];
          comp_time += recoloringPerRound[distributedRounds];
          compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
          totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	} else if(timing) {
	  double recolor_round_time = timer() - recolor_temp;
	  recoloring_time += recolor_round_time;
	  comp_time += recolor_round_time;
	}
        
	//reset the recoloringSize device host and device views
	//to zero
        recoloringSize_host(0) = 0;
        Kokkos::deep_copy(recoloringSize,recoloringSize_host);

        Kokkos::parallel_for("set femv colors",
          Kokkos::RangePolicy<execution_space, int>(0,rand.size()-nVtx), 
          KOKKOS_LAMBDA(const int& i){
            femv_colors(i+nVtx) = ghost_colors(i);
          });
        Kokkos::fence();
        //communicate
        Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
	Kokkos::deep_copy(verts_to_send_size_host, verts_to_send_size);
	gno_t sent,recv;
        // Reset device views
        femvColors = decltype(femvColors)();
        femv_colors = decltype(femv_colors)();

	double curr_comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,
			                                  nVtx, 
							  verts_to_send_host, 
							  verts_to_send_size_host,
							  femv,
							  procs_to_send,
							  recv,
							  sent); 	
        comm_time += curr_comm_time;
        if(distributedRounds < numStatisticRecordingRounds){
	  commPerRound[distributedRounds] = curr_comm_time;
          sentPerRound[distributedRounds] = sent;
          recvPerRound[distributedRounds] = recv;
          totalPerRound[distributedRounds] += commPerRound[distributedRounds];
        }
        //detect conflicts in parallel. For a detected conflict,
        //reset the vertex-to-be-recolored's color to 0, in order to
        //allow KokkosKernels to recolor correctly.
	
        femvColors = femv->getLocalViewDevice(Tpetra::Access::ReadWrite);
        femv_colors = subview(femvColors, Kokkos::ALL, 0);
        Kokkos::parallel_for("get femv colors 2",
          Kokkos::RangePolicy<execution_space, int>(0,rand.size()-nVtx),
          KOKKOS_LAMBDA(const int& i){
            ghost_colors(i) = femv_colors(i+nVtx);
          });
        Kokkos::fence();
	verts_to_send_size_host(0) = 0;
	deep_copy(verts_to_send_size, verts_to_send_size_host);
        double detection_temp = timer();
	detectConflicts<execution_space, memory_space>(nVtx,
			                               rand.size()-nVtx,
						       dist_offsets,
						       dist_adjs,
						       femv_colors,
						       verts_to_send_view,
						       verts_to_send_size_atomic,
						       recoloringSize,
						       rand_dev,
						       gid_dev,
						       ghost_degrees_dev,
						       recolor_degrees);

	Kokkos::deep_copy(recoloringSize_host, recoloringSize);
        
	if(distributedRounds < numStatisticRecordingRounds){
          conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
          conflict_detection += conflictDetectionPerRound[distributedRounds];
          compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
          totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
          comp_time += conflictDetectionPerRound[distributedRounds];
	} else if(timing){
	  double conflict_detection_round_time = timer()- detection_temp;
	  conflict_detection += conflict_detection_round_time;
	  comp_time += conflict_detection_round_time;
	}
        //do a reduction to determine if we're done
        int globalDone = 0;
        int localDone = recoloringSize_host(0);
        Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
        //We're only allowed to stop once everyone has no work to do.
        //collectives will hang if one process exits. 
        distributedRounds++;
        done = !globalDone;
      }

      if(recoloringSize_host(0) > 0 || !done){
	ghost_colors_host = Kokkos::create_mirror_view(ghost_colors);
	deep_copy(ghost_colors_host, ghost_colors);
	deep_copy(verts_to_send_host, verts_to_send_view);
	deep_copy(verts_to_send_size_host, verts_to_send_size);
      }

      
      //finish the local coloring in serial on the host
      while(recoloringSize_host(0) > 0 || !done){
	//Use non-templated call to get the Host view
	auto femvColors = femv->getLocalViewHost(Tpetra::Access::ReadWrite);
	auto femv_colors = subview(femvColors, Kokkos::ALL, 0);
	/*Kokkos::View<int*, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>>femv_colors_host = create_mirror(femv_colors);
	Kokkos::deep_copy(femv_colors_host, femv_colors);*/
	if(distributedRounds < 100){
	  vertsPerRound[distributedRounds] = recoloringSize_host(0);
	}

	double recolor_temp = timer();
	//use KokkosKernels to recolor the conflicting vertices
	if(verts_to_send_size_host(0) > 0){
	  this->colorInterior<host_exec,
			      host_mem>
			      (femv_colors.size(), dist_adjs_host, dist_offsets_host, femv, verts_to_send_host, verts_to_send_size_host(0), true);
	}

	if(distributedRounds < numStatisticRecordingRounds){
	  recoloringPerRound[distributedRounds] = timer() - recolor_temp;
	  recoloring_time += recoloringPerRound[distributedRounds];
	  comp_time += recoloringPerRound[distributedRounds];
	  compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	  totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        } else if(timing){
	  double recolor_serial_round_time = timer() - recolor_temp;
	  recoloring_time += recolor_serial_round_time;
	  comp_time += recolor_serial_round_time;
	}

	recoloringSize_host(0) = 0;

	for(size_t i = 0; i < rand.size() -nVtx; i++){
	  femv_colors(i+nVtx) = ghost_colors_host(i);
	}

	gno_t sent,recv;
	double curr_comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,
			                                  nVtx, 
							  verts_to_send_host, 
							  verts_to_send_size_host, 
							  femv,
							  procs_to_send, 
							  recv, 
							  sent);
	comm_time += curr_comm_time;

	if(distributedRounds < numStatisticRecordingRounds){
	  commPerRound[distributedRounds] = curr_comm_time;
	  sentPerRound[distributedRounds] = sent;
	  recvPerRound[distributedRounds] = recv;
	  totalPerRound[distributedRounds] += commPerRound[distributedRounds];
	}
	for(size_t i = 0; i < rand.size()-nVtx; i++){
	  ghost_colors_host(i) = femv_colors(i+nVtx);
	}

	verts_to_send_size_host(0) = 0;
	double detection_temp = timer();
	detectConflicts<host_exec, host_mem>(nVtx,
			                     rand.size()-nVtx,
					     dist_offsets_host,
					     dist_adjs_host,
					     femv_colors,
					     verts_to_send_host,
					     verts_to_send_size_host,
					     recoloringSize_host,
					     rand_host,
					     gid_host,
					     ghost_degrees_host,
					     recolor_degrees);
        if(distributedRounds < numStatisticRecordingRounds){
	  conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
	  conflict_detection += conflictDetectionPerRound[distributedRounds];
	  compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
	  totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
	  comp_time += conflictDetectionPerRound[distributedRounds];
        } else if(timing){
	  double conflict_detection_serial_round_time = timer() - detection_temp;
	  conflict_detection += conflict_detection_serial_round_time;
	  comp_time += conflict_detection_serial_round_time;
	}
	//do a reduction to determine if we're done
	int globalDone = 0;
	int localDone = recoloringSize_host(0);
	Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
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
        //if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n", distributedRounds);
	int totalBoundarySize = 0;
        int totalVertsPerRound[numStatisticRecordingRounds];
        double finalTotalPerRound[numStatisticRecordingRounds];
        double maxRecoloringPerRound[numStatisticRecordingRounds];
        double finalSerialRecoloringPerRound[numStatisticRecordingRounds];
        double minRecoloringPerRound[numStatisticRecordingRounds];
        double finalCommPerRound[numStatisticRecordingRounds];
        double finalCompPerRound[numStatisticRecordingRounds];
        double finalConflictDetectionPerRound[numStatisticRecordingRounds];
        gno_t finalRecvPerRound[numStatisticRecordingRounds];
        gno_t finalSentPerRound[numStatisticRecordingRounds];
        for(int i = 0; i < numStatisticRecordingRounds; i++) {
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
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,vertsPerRound,totalVertsPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,totalPerRound,finalTotalPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,recoloringPerRound,maxRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,numStatisticRecordingRounds,recoloringPerRound,minRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,serialRecoloringPerRound,finalSerialRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,commPerRound,finalCommPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,compPerRound,finalCompPerRound);
        Teuchos::reduceAll<int,double>(*comm, 
			               Teuchos::REDUCE_MAX,numStatisticRecordingRounds,conflictDetectionPerRound,finalConflictDetectionPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,recvPerRound, finalRecvPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,sentPerRound, finalSentPerRound);
        
        std::cout << "Rank " << comm->getRank() 
                  << ": boundary size: " << localBoundaryVertices << std::endl;
        if(comm->getRank()==0) 
          std::cout << "Total boundary size: " << totalBoundarySize << std::endl;
        for(int i = 0; i < std::min(distributedRounds,numStatisticRecordingRounds); i++){
          std::cout << "Rank " << comm->getRank() 
                    << ": recolor " << vertsPerRound[i]
                    << " vertices in round " << i << std::endl;
          if(comm->getRank()==0) {
            std::cout << "recolored " << totalVertsPerRound[i]
                      << " vertices in round " << i << std::endl;
            std::cout << "total time in round " << i 
                      << ":  " << finalTotalPerRound[i] << std::endl;;
            std::cout << "recoloring time in round " << i
                      << ": " << maxRecoloringPerRound[i] << std::endl;
            std::cout << "serial recoloring time in round " << i
                      << ": " << finalSerialRecoloringPerRound[i] << std::endl;
            std::cout << "min recoloring time in round " << i
                      << ": " << minRecoloringPerRound[i] << std::endl;
            std::cout << "conflict detection time in round " << i
                      << ": " << finalConflictDetectionPerRound[i] << std::endl;
            std::cout << "comm time in round " << i
                      << ": " << finalCommPerRound[i] << std::endl;
            std::cout << "total sent in round " << i
                      << ": " << finalSentPerRound[i] << std::endl;
            std::cout << "total recv in round " << i
                      << ": " << finalRecvPerRound[i] << std::endl;
            std::cout << "comp time in round " << i
                      << ": " << finalCompPerRound[i] << std::endl;
          }
        }
      } else if(timing){
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
          std::cout << "Total Time: " << global_total_time << std::endl;
          std::cout << "Interior Time: " << global_interior_time << std::endl;
          std::cout << "Recoloring Time: " << global_recoloring_time << std::endl;
          std::cout << "Min Recoloring Time: " << global_min_recoloring_time << std::endl;
          std::cout << "Conflict Detection Time: " << global_conflict_detection << std::endl;
          std::cout << "Comm Time: " << global_comm_time << std::endl;
          std::cout << "Comp Time: " << global_comp_time << std::endl;
        }
      }
      if(verbose) std::cout<<comm->getRank()<<": exiting coloring\n";
    }
};

template <typename Adapter>
void AlgDistance1<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
    
  this->env->debug(DETAILED_STATUS, "   building graph model");
  if(verbose) std::cout<<comm->getRank()<<": starting to construct graph model\n";
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
  if(verbose) std::cout<<comm->getRank()<<": done constructing graph model\n";
  this->env->debug(DETAILED_STATUS, "   graph model built");
}


}//end namespace Zoltan2
#endif
