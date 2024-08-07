// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_2GHOSTLAYER_HPP_
#define _ZOLTAN2_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

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
//! \file Zoltan2_2GhostLayer.hpp
//! \brief Abstract base class for any coloring algorithm
//         for methods that use two layers of ghosts.


namespace Zoltan2{

template <typename Adapter>
class AlgTwoGhostLayer : public Algorithm<Adapter> {

  public:

    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t,gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>;
    using device_type = typename femv_t::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using host_exec = typename femv_t::host_view_type::device_type::execution_space;
    using host_mem = typename femv_t::host_view_type::device_type::memory_space;

    double timer(){
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
    }
  private:

    //Entry point for parallel local coloring
    //  nVtx is the number of vertices owned by the current process
    //
    //  adjs_view is the adjacencies, indexed by offset_view
    //
    //  offset_view is the CSR row map, used to index the adjs_view
    //
    //  femv is the FEMultiVector that holds the colors for the vertices
    //       the colors will change in this function.
    //
    //  vertex_list is a list of vertices to recolor
    //
    //  vertex_list_size is the size of the list of vertices to recolor
    //                   vertex_list_size = 0 means recolor all uncolored vertices
    //
    //  recolor decides which KokkosKernels algorithm to use.
    //
    virtual void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*, device_type > adjs_view,
                       Kokkos::View<offset_t*,device_type > offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, device_type> vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false) = 0;

    //Entry point for serial local coloring
    virtual void colorInterior_serial(const size_t nVtx,
                       typename Kokkos::View<lno_t*, device_type >::HostMirror adjs_view,
                       typename Kokkos::View<offset_t*,device_type >::HostMirror offset_view,
                       Teuchos::RCP<femv_t> femv,
		       typename Kokkos::View<lno_t*, device_type>::HostMirror vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false) = 0;
  public:
    //Entry point for parallel conflict detection
    //INPUT ARGS
    //  n_local is the number of vertices owned by the current process
    //
    //  dist_offsets_dev is the device view that holds the CSR offsets
    //                   It holds offsets for all vertices, owned and ghosts.
    //                   It is organized with owned first, and then ghost layers:
    //                   -----------------------------------------------
    //                   | owned     | first       | second            |
    //                   |  vtx      | ghost       | ghost             |
    //                   | offsets   | layer       | layer             |
    //                   -----------------------------------------------
    //                   The allocated length is equal to the sum of owned
    //                   and ghost vertices + 1. Any vertex with LID >= n_local
    //                   is a ghost.
    //
    //  dist_adjs_dev is the device view that holds CSR adjacencies
    //
    //  boundary_verts_view holds the local IDs of vertices in the boundary.
    //                      By definition, boundary vertices are owned vertices
    //                      that have a ghost in their neighborhood.
    //                      (distance-1 coloring = 1-hop neighborhood)
    //                      (distance-2 coloring = 2-hop neighborhood)
    //                      Note: for D1-2GL we do not use this argument.
    //                      It is possible to detect all D1 conflicts by only
    //                      checking the ghost vertices, without constructing
    //                      the boundary.
    //
    //  rand is a device view that holds random numbers generated from GIDs,
    //       they are consistently generated across processes
    //
    //  gid is a device view that holds the GIDs of each vertex on this process.
    //      It is indexable by local ID. It stores both owned and ghost
    //      vertices, and LIDs are ordered so that the structure looks like:
    //      -----------------------------------------
    //      |          |  first    | second         | owned LIDs are smallest
    //      |  owned   |  ghost    | ghost          | ghost LIDs increase based
    //      |   vtx    |  layer    | layer          | on layer (1 < 2)
    //      -----------------------------------------
    //      The allocated size is dist_offsets_dev.extent(0)-1.
    //
    //
    //  ghost_degrees is a device view that holds the degrees of ghost vertices only.
    //                A ghost with local ID n_local will be the first entry in this view.
    //                So, ghost_degrees[i] is the degree of the vtx with
    //                GID = gid[n_local+i].
    //
    //  recolor_degrees is a boolean that determines whether or not we factor in vertex
    //                  degrees on recoloring
    //
    //OUTPUT ARGS
    //  femv_colors is the device color view.
    //              After this function call, conflicts between vertices will
    //              be resolved via setting a vertex's color to 0. The vertex
    //              to be uncolored is determined by rules that are consistent
    //              across multiple processes.
    //
    //  verts_to_recolor_view is a device view that holds the list
    //                          of vertices to recolor. Any vertex
    //                          uncolored by this function will appear in this
    //                          view after the function returns.
    //
    //  verts_to_recolor_size_atomic is an atomic device view that holds the
    //                               size of verts_to_recolor_view.
    //                               Effectively it represents how many
    //                               vertices need to be recolored after
    //                               conflict detection.
    //                               It differs in size from the allocated size of
    //                               verts_to_recolor_view.
    //                               Atomicity is required for building
    //                               verts_to_recolor_view in parallel, which
    //                               is the reason this is a view instead of
    //                               just an integer type.
    //
    //  verts_to_send_view is a device view that holds the list of vertices
    //                       that will need to be sent to remotes after recoloring
    //                       Note: Only locally owned vertices should ever be in
    //                       this view. A process cannot color ghosts correctly.
    //
    //  verts_to_send_size_atomic is an atomic device view that holds the size of
    //                            verts_to_send_view. It differs in size from the
    //                            allocated size of verts_to_send_view.
    //                            Atomicity is required for building
    //                            verts_to_send_view in parallel,
    //                            which is the reason this is a view instead of
    //                            just an integer type.
    //
    //  recoloringSize is a device view that stores the total number of
    //                 vertices that were uncolored by the conflict detection.
    //
    virtual void detectConflicts(const size_t n_local,
		                 Kokkos::View<offset_t*, device_type > dist_offsets_dev,
		                 Kokkos::View<lno_t*, device_type > dist_adjs_dev,
				 Kokkos::View<int*,device_type > femv_colors,
				 Kokkos::View<lno_t*, device_type > boundary_verts_view,
                                 Kokkos::View<lno_t*,
				              device_type > verts_to_recolor_view,
				 Kokkos::View<int*,
				              device_type,
					      Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_size_atomic,
				 Kokkos::View<lno_t*,
				              device_type > verts_to_send_view,
				 Kokkos::View<size_t*,
				              device_type,
					      Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic,
			         Kokkos::View<size_t*, device_type> recoloringSize,
				 Kokkos::View<int*,
				              device_type> rand,
			         Kokkos::View<gno_t*,
				              device_type> gid,
                                 Kokkos::View<gno_t*,
				              device_type> ghost_degrees,
			         bool recolor_degrees) = 0;

    //Entry point for serial conflict detection
    virtual void detectConflicts_serial(const size_t n_local,
		                 typename Kokkos::View<offset_t*, device_type >::HostMirror dist_offsets_host,
		                 typename Kokkos::View<lno_t*, device_type >::HostMirror dist_adjs_host,
				 typename Kokkos::View<int*,device_type >::HostMirror femv_colors,
				 typename Kokkos::View<lno_t*, device_type >::HostMirror boundary_verts_view,
                                 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_recolor_view,
				 typename Kokkos::View<int*,device_type>::HostMirror verts_to_recolor_size_atomic,
				 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_send_view,
				 typename Kokkos::View<size_t*,device_type>::HostMirror verts_to_send_size_atomic,
				 typename Kokkos::View<size_t*, device_type>::HostMirror recoloringSize,
			         typename Kokkos::View<int*,  device_type>::HostMirror rand,
			         typename Kokkos::View<gno_t*,device_type>::HostMirror gid,
                                 typename Kokkos::View<gno_t*,device_type>::HostMirror ghost_degrees,
				 bool recolor_degrees) = 0;
    //Entry point for the construction of the boundary
    //INPUT ARGS
    //  n_local is the number of vertices owned by the current process
    //
    //  dist_offsets_dev is the device view that holds the CSR offsets
    //
    //  dist_adjs_dev is the device view that holds CSR adjacencies
    //
    //  dist_offsets_host is the hostmirror that holds the CSR offsets
    //
    //  dist_adjs_host is the hostmirror that holds the CSR adjacencies
    //
    //OUTPUT ARGS
    //  boundary_verts is an unallocated device view that will hold the
    //                 list of boundary vertices.
    //
    //  verts_to_send_view will hold the list of vertices to send
    //
    //  verts_to_send_size_atomic will hold the number of vertices to
    //                            send currently held in verts_to_send_view.
    //                            verts_to_send_size_atomic differs
    //                            from the allocated size of verts_to_send_view
    //                            Atomicity is required for building
    //                            verts_to_send_view in parallel, which is
    //                            the reason this is a view instead of an
    //                            integer type.
    //
    virtual void constructBoundary(const size_t n_local,
		    		   Kokkos::View<offset_t*, device_type> dist_offsets_dev,
		                   Kokkos::View<lno_t*, device_type> dist_adjs_dev,
				   typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host,
				   typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host,
                                   Kokkos::View<lno_t*, device_type>& boundary_verts,
				   Kokkos::View<lno_t*,
				                device_type > verts_to_send_view,
				   Kokkos::View<size_t*,
				                device_type,
						Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic) = 0;

  protected:
    RCP<const base_adapter_t> adapter;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    bool verbose;
    bool timing;

  private:
    //This function constructs a CSR with complete adjacency information for
    //first-layer ghost vertices. This CSR can contain edges from:
    //      first-layer ghosts to owned;
    //      first-layer ghosts to first-layer ghosts;
    //      first-layer ghosts to second-layer ghosts;
    //1) Each process sends a list of ghost GIDs back to their owning process
    //   to request full adjacency information for that vertex.
    //
    //2) Initially each process sends back degree information to the requesting
    //   process.
    //
    //3) Then, each process can construct send buffers and receive schedules for
    //   the adjacency information for each requested GID it received,
    //   in addition to building the 2GL offsets and relabeling ghost LIDs
    //   to make the final CSR construction easier.
    //
    //   3a)  We can construct the 2GL offsets here because each process has
    //        already received the degree information for each first-layer
    //        ghost vertex.
    //
    //   3b)  The original ghost LIDs may not correspond to the order in which
    //        we will receive the adjacency information. Instead of reordering
    //        the received adjacencies, we relabel the GIDs of first-layer
    //        ghosts with new LIDs that correspond to the order in which we
    //        receive the adjacency information. Only first-layer ghost LIDs
    //        are changed.
    //
    //4) Due to limitations on the size of an MPI message, we split sending and
    //   receiving into rounds to accommodate arbitrarily large graphs.
    //
    //5) As send rounds progress, we construct the 2GL adjacency structure.
    //   Because we relabel ghost LIDs, we do not need to reorder the
    //   data after we receive it.
    //
    //
    //
    //OUTPUT ARGS
    //  ownedPlusGhosts Initially holds the list of GIDs for owned and ghost
    //                  vertices. The only hard constraint on the ordering is
    //                  that owned vertices come before ghost vertices.
    //                  This function can re-assign the LIDs of ghost vertices
    //                  in order to make the final construction of the 2GL
    //                  CSR easier. After this function call, some ghost
    //                  LIDs may be changed, but they will still be greater
    //                  than owned LIDs. No owned LIDs will be changed.
    //
    //  adjs_2GL holds the second ghost layer CSR adjacencies. The 2GL CSR
    //           contains complete adjacency information for the first-layer
    //           ghosts. These adjacencies can hold both owned vertices and
    //           second-layer ghosts, as well as first-layer ghosts.
    //
    //  offsets_2GL holds CSR offsets for vertices in the first ghost layer only.
    //              That is, a vertex with GID = gid[n_local] will be the first
    //              vertex with information in this offset structure.
    //
    //
    //INPUT ARGS
    //  owners holds the owning proc for a given vertex, indexed by local ID.
    //
    //  adjs is the CSR adjacencies of the local graph with a single ghost layer.
    //
    //  offsets is the CSR offsets of the local graph with a single ghost layer.
    //
    //  mapOwned translates from Owned GID to LID. We only need this translation
    //           for owned vertices.
    //
    //TODO: This function uses many vectors of size comm->getSize();
    //      consider changes to increase memory scalability.
    void constructSecondGhostLayer(std::vector<gno_t>& ownedPlusGhosts,
                                   const std::vector<int>& owners,
                                   ArrayView<const gno_t> adjs,
                                   ArrayView<const offset_t> offsets,
                                   RCP<const map_t> mapOwned,
                                   std::vector< gno_t>& adjs_2GL,
                                   std::vector< offset_t>& offsets_2GL) {

      //first, we send ghost GIDs back to owners to receive the
      //degrees of each first-layer ghost.
      std::vector<int> sendcounts(comm->getSize(),0);
      std::vector<size_t> sdispls(comm->getSize()+1,0);
      //loop through owners, count how many vertices we'll send to each processor
      //We send each ghost GID back to its owning process.
      if(verbose) std::cout<<comm->getRank()<<": building sendcounts\n";
      for(size_t i = 0; i < owners.size(); i++){
        if(owners[i] != comm->getRank()&& owners[i] !=-1) sendcounts[owners[i]]++;
      }
      //construct sdispls (for building sendbuf), and sum the total sendcount
      if(verbose) std::cout<<comm->getRank()<<": building sdispls\n";
      size_t sendcount = 0;
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
      //compute send and receive schedules to and from each process
      std::vector<std::vector<uint64_t>> per_proc_round_adj_sums(max_rounds+1,std::vector<uint64_t>(comm->getSize(),0));
      std::vector<std::vector<uint64_t>> per_proc_round_vtx_sums(max_rounds+1,std::vector<uint64_t>(comm->getSize(),0));

      if(verbose) std::cout<<comm->getRank()<<": filling per_proc sums\n";
      //fill out the send schedules
      for(int proc_to_send = 0; proc_to_send < comm->getSize(); proc_to_send++){
        int curr_round = 0;
        for(size_t j = sdispls[proc_to_send]; j < sdispls[proc_to_send+1]; j++){
          if((per_proc_round_adj_sums[curr_round][proc_to_send] + recvDegrees[j])*sizeof(gno_t) > INT_MAX){
            curr_round++;
          }
          per_proc_round_adj_sums[curr_round][proc_to_send] += recvDegrees[j];
          per_proc_round_vtx_sums[curr_round][proc_to_send]++;
        }
      }

      if(verbose) std::cout<<comm->getRank()<<": building recv GID schedule\n";

      //A 3D vector to hold the GIDs so we can see how this process will receive
      //the vertices in each round, from each process. This way we can reorder things
      //locally so that the data we receive will automatically construct the CSR correctly
      //without any other computation. We do the reordering before we start receiving.
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
        for(size_t j = sdispls[i]; j < sdispls[i+1]; j++){
          if(curr_idx > per_proc_round_vtx_sums[curr_round][i]){
            curr_round++;
            curr_idx = 0;
          }
          recv_GID_per_proc_per_round[curr_round][i][curr_idx++] = j;
        }
      }

      if(verbose) std::cout<<comm->getRank()<<": reordering gids and degrees in the order they'll be received\n";
      //reorder the GIDs and degrees locally:
      //  - The way we previously stored GIDs and degrees had no relation
      //    to the order that we'll be receiving the adjacency data.
      //    For the CSR to be complete (and usable), we reorder the LIDs of
      //    ghosts so that they are in the order we will receive the adjacency
      //    data.
      //
      //  - For graphs that are not large enough to trigger sending in multiple
      //    rounds of communication, the LIDs of the ghosts will be ordered
      //    by owning process. If multiple rounds are used, the LIDs of
      //    ghosts will be ordered by rounds in addition to owning process.
      //
      //  - final_gid_vec and final_degree_vec hold the reorganized gids and
      //    degrees, the final re-mapping happens further down
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

      //check to see if the reorganization has happened correctly
      if(verbose){
        //do a quick check to see if we ended up reorganizing anything
        bool reorganized = false;
        for(size_t i = 0; i < sendcount; i++){
          if(final_gid_vec[i] != sendbuf[i]) reorganized = true;
        }

        //if we have more than a single round of communication, we need to reorganize.
        //this alerts of unnecessary reorganization, and a failure to perform necessary
        //reorganization.
        if(!reorganized && (max_rounds > 1)) std::cout<<comm->getRank()<<": did not reorgainze GIDs, but probably should have\n";
        if(reorganized && (max_rounds == 1)) std::cout<<comm->getRank()<<": reorganized GIDs, but probably should not have\n";
      }
      //remap the GIDs so we receive the adjacencies in the same order as the current processes LIDs
      //  originally, the GID->LID mapping has no relevance to how we'll be receiving adj data from
      //  remote processes. Here we make it so that the GID->LID mapping does correspond to the
      //  order we have received degree info and will receive adjacencies. ( LID n_local
      //  corresponds to the first GID whose adjacency info we will receive, LID n_local+1 the
      //  second, etc.)
      //  That way, we don't need to reorder the adjacencies after we receive them.
      //
      //  This next loop is the final mapping step.

      for (size_t i = 0; i < sendcount; i++){
        ownedPlusGhosts[i+offsets.size()-1] = final_gid_vec[i];
      }

      //status report
      if(verbose) {
        std::cout<<comm->getRank()<<": done remapping\n";
        std::cout<<comm->getRank()<<": building ghost offsets\n";
      }
      //start building the second ghost layer's offsets
      std::vector<offset_t> ghost_offsets(sendcount+1,0);
      for(size_t i = 1; i < sendcount+1; i++){
        ghost_offsets[i] = ghost_offsets[i-1] + final_degree_vec[i-1];
      }


      if(verbose) std::cout<<comm->getRank()<<": going through the sending rounds\n";
      //set up counters to keep track of where we are in the sending order
      std::vector<uint64_t> curr_idx_per_proc(comm->getSize(),0);
      for(int i = 0; i < comm->getSize(); i++) curr_idx_per_proc[i] = rdispls[i];
      for(int round = 0; round < max_rounds; round++){
	//send buffers
        std::vector<gno_t> send_adj;
        std::vector<int> send_adj_counts(comm->getSize(),0);
        if(verbose) std::cout<<comm->getRank()<<": round "<<round<<", constructing send_adj\n";
	//construct the adjacencies to send for this round
        for(int curr_proc = 0; curr_proc < comm->getSize(); curr_proc++){
          uint64_t curr_adj_sum = 0;
	  //keep going through adjacencies to send to this process until we're done
          while( curr_idx_per_proc[curr_proc] < (size_t)rdispls[curr_proc+1]){
            lno_t lid = mapOwned->getLocalElement(recvbuf[curr_idx_per_proc[curr_proc]++]);

	    //if the next adjacency would push us over the MPI message size max,
	    //stop for this round
            if((curr_adj_sum + (offsets[lid+1]-offsets[lid]))*sizeof(gno_t) >= INT_MAX){
              break;
            }

	    //add the adjacencies to the send buffer
            curr_adj_sum += (offsets[lid+1] - offsets[lid]);
            for(offset_t j = offsets[lid]; j < offsets[lid+1]; j++){
              send_adj.push_back(adjs[j]);
            }
          }
	  //update the send counts for this round
          send_adj_counts[curr_proc] = curr_adj_sum;
        }
        if(verbose) std::cout<<comm->getRank()<<": round "<<round<<", sending...\n";
	//do the sending...
        Teuchos::ArrayView<gno_t> send_adjs_view = Teuchos::arrayViewFromVector(send_adj);
        Teuchos::ArrayView<int> adjsendcounts_view = Teuchos::arrayViewFromVector(send_adj_counts);
        Teuchos::ArrayRCP<gno_t> ghost_adjs;
        std::vector<int> adjrecvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> adjsrecvcounts_view = Teuchos::arrayViewFromVector(adjrecvcounts);
        Zoltan2::AlltoAllv<gno_t>(*comm, *env, send_adjs_view, adjsendcounts_view, ghost_adjs, adjsrecvcounts_view);
	//Because of the previous reordering, these adjacencies are
	//in the correct order as they arrive on the process.
        for(offset_t i = 0; i< (offset_t)ghost_adjs.size(); i++){
          adjs_2GL.push_back(ghost_adjs[i]);
        }
      }
      if(verbose) std::cout<<comm->getRank()<<": constructing offsets\n";
      //put the offsets we computed into the output argument.
      for(size_t i = 0; i < sendcount+1; i++){
        offsets_2GL.push_back(ghost_offsets[i]);
      }
      if(verbose) std::cout<<comm->getRank()<<": done building 2nd ghost layer\n";
    }

    //Communicates owned vertex colors to remote ghost copies.
    //
    //returns: the amount of time the communication call took.
    //
    //OUTPUT ARGS:
    //  colors: the owned vertices' colors are not changed,
    //          ghost vertex colors are updated from their owners.
    //
    //  total_sent: reports the total size of the send buffer
    //
    //  total_recv: reports the total size of the recv buffer
    //
    //INPUT ARGS:
    //
    //  mapOwnedPlusGhosts: maps global IDs to local IDs and vice-versa.
    //
    //  nVtx: the number of owned vertices on this process
    //
    //  verts_to_send: hostmirror of the list of vertices to send.
    //                 This list will only contain local vertices
    //                 that are ghosted on a remote process.
    //
    //  verts_to_send_size: hostmirror of the size of verts_to_send
    //
    //  procs_to_send: map that takes a local ID and gives a vector of
    //                 process IDs which have a ghost copy of that vertex.
    //
    double doOwnedToGhosts(RCP<const map_t> mapOwnedPlusGhosts,
                           size_t nVtx,
			   typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_send,
			   typename Kokkos::View<size_t*,device_type>::HostMirror verts_to_send_size,
                           Teuchos::RCP<femv_t> femv,
			   const std::unordered_map<lno_t, std::vector<int>>& procs_to_send,
                           gno_t& total_sent, gno_t& total_recvd){

      auto femvColors = femv->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto colors = subview(femvColors, Kokkos::ALL, 0);
      //create vectors to hold send information
      int nprocs = comm->getSize();
      std::vector<int> sendcnts(comm->getSize(), 0);
      std::vector<gno_t> sdispls(comm->getSize()+1, 0);

      //calculate how much data we're sending to each process
      for(size_t i = 0; i < verts_to_send_size(0); i++){
        for(size_t j = 0; j < procs_to_send.at(verts_to_send(i)).size(); j++){
	  sendcnts[procs_to_send.at(verts_to_send(i))[j]] += 2;
	}
      }
      //calculate sendsize and sdispls
      sdispls[0] = 0;
      gno_t sendsize = 0;
      std::vector<int> sentcount(nprocs, 0);

      for(int i = 1; i < comm->getSize()+1; i++){
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
	sendsize += sendcnts[i-1];
      }
      total_sent = sendsize;
      std::vector<gno_t> sendbuf(sendsize,0);
      //construct sendbuf, send each owned vertex's GID
      //and its color to the processes that have a
      //ghost copy of that vertex. If a vertex is not ghosted,
      //it does not get sent anywhere.
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

      //if we're reporting times, remove the computation imbalance from the comm timer
      if(timing) comm->barrier();
      double comm_total = 0.0;
      double comm_temp = timer();

      Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcnts_view, recvbuf, recvcnts_view);
      comm_total += timer() - comm_temp;

      //compute total recvsize for updating local ghost colors
      gno_t recvsize = 0;
      for(int i = 0; i < recvcnts_view.size(); i++){
        recvsize += recvcnts_view[i];
      }
      total_recvd = recvsize;
      //update the local ghost copies with the color we just received.
      for(int i = 0; i < recvsize; i+=2){
        size_t lid = mapOwnedPlusGhosts->getLocalElement(recvbuf[i]);
	colors(lid) = recvbuf[i+1];
      }

      return comm_total;
    }

  public:
    //constructor
    AlgTwoGhostLayer(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_){
      verbose = pl->get<bool>("verbose",false);
      timing = pl->get<bool>("timing", false);
    }
    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution){
      //convert from global graph to local graph
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;

      modelFlag_t flags;
      flags.set(REMOVE_SELF_EDGES);
      const auto model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      // the weights are not used at this point.

      //get the adjacencies in a view that holds GIDs
      ArrayView<const gno_t> adjs;
      //get the CSR offsets
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      model->getEdgeList(adjs, offsets, ewgts);

      //This map is used to map GID to LID initially
      std::unordered_map<gno_t,lno_t> globalToLocal;
      //this vector holds GIDs for owned and ghost vertices
      //The final order of the gids is:
      // -----------------------------------
      // |           | first    | second   |
      // | owned     | layer    | layer    |
      // |           | ghosts   | ghosts   |
      // -----------------------------------
      //
      // constructSecondGhostLayer reorders the first layer ghost
      // GIDs in the particular order that the communication requires.

      std::vector<gno_t> ownedPlusGhosts;

      //this vector holds the owning process for
      //owned vertices and the first ghost layer.
      //owners[i] gives the processor owning ownedPlusGhosts[i],
      //for owned vertices and first layer ghosts.
      //
      //This should only be used BEFORE the call to constructSecondGhostLayer
      std::vector<int> owners;

      //fill these structures using owned vertices
      for(int i = 0; i < vtxIDs.size(); i++){
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
        owners.push_back(comm->getRank());
      }

      //count the initial number of ghosts,
      //map them to a local ID.
      //The globalToLocal map has first layer ghosts
      //from here on.
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

      //construct a Tpetra map for the FEMultiVector
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                           <Tpetra::global_size_t>::invalid();
      RCP<const map_t> mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));

      //GID and owner lookup for ghost vertices
      std::vector<gno_t> ghosts;
      std::vector<int> ghostowners;
      for(size_t i = nVtx; i < nVtx+nGhosts; i++){
        ghosts.push_back(ownedPlusGhosts[i]);
        ghostowners.push_back(-1);
      }

      //get the owners of the ghost vertices
      ArrayView<int> owningProcs = Teuchos::arrayViewFromVector(ghostowners);
      ArrayView<const gno_t> gids = Teuchos::arrayViewFromVector(ghosts);
      mapOwned->getRemoteIndexList(gids, owningProcs);

      //add ghost owners to the total owner vector
      for(size_t i = 0; i < ghostowners.size(); i++){
        owners.push_back(ghostowners[i]);
      }

      //construct the second ghost layer.
      //NOTE: this may reorder the LIDs of the ghosts.
      //
      //these vectors hold the CSR representation of the
      // first ghost layer. This CSR contains:
      //   - edges from first layer ghosts to:
      //       - owned vertices
      //       - first layer ghosts
      //       - second layer ghosts
      //
      //   - first_layer_ghost_adjs uses GIDs that need
      //     translated to LIDs.
      //
      //   - first_layer_ghost_offsets are indexed by LID,
      //     first_layer_ghost_offsets[i] corresponds to
      //     the vertex with GID = ownedPlusGhosts[i+nVtx].
      //     first_layer_ghost_offsets.size() = #first layer ghosts + 1
      //
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

      //use updated globalToLocal map to translate
      //adjacency GIDs to LIDs
      for(int i = 0 ; i < adjs.size(); i++){
        local_adjs[i] = globalToLocal[adjs[i]];
      }

      //at this point, we have ownedPlusGhosts with 1layer ghosts' GIDs.
      //Need to map the second ghost layer GIDs to new local IDs,
      //and add them to the map, as well as translating 2layer ghost
      //adjacencies to use those new LIDs.
      size_t n2Ghosts = 0;

      //local_ghost_adjs is the LID translation of first_layer_ghost_adjs.
      std::vector<lno_t> local_ghost_adjs;
      for(size_t i = 0; i< first_layer_ghost_adjs.size(); i++ ){
        if(globalToLocal.count(first_layer_ghost_adjs[i]) == 0){
          ownedPlusGhosts.push_back(first_layer_ghost_adjs[i]);
          globalToLocal[first_layer_ghost_adjs[i]] = vtxIDs.size() + nGhosts + n2Ghosts;
          n2Ghosts++;
        }
        local_ghost_adjs.push_back(globalToLocal[first_layer_ghost_adjs[i]]);
      }

      //construct data structures necessary for FEMultiVector
      if(verbose) std::cout<<comm->getRank()<<": constructing Tpetra map with copies\n";
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

      //Create random numbers seeded on global IDs to resolve conflicts
      //consistently between processes.
      std::vector<int> rand(ownedPlusGhosts.size());
      for(size_t i = 0; i < rand.size(); i++){
        std::srand(ownedPlusGhosts[i]);
        rand[i] = std::rand();
      }

      // get up-to-date owners for all ghosts.
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
	//do the communication
        Teuchos::ArrayView<int> sendcounts_view = Teuchos::arrayViewFromVector(sendcounts);
        Teuchos::ArrayView<gno_t> sendbuf_view = Teuchos::arrayViewFromVector(sendbuf);
        Teuchos::ArrayRCP<gno_t>  recvbuf;
        std::vector<int> recvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> recvcounts_view = Teuchos::arrayViewFromVector(recvcounts);
        Zoltan2::AlltoAllv<gno_t>(*comm, *env, sendbuf_view, sendcounts_view, recvbuf, recvcounts_view);
        std::vector<int> is_bndry_send(recvbuf.size(),0);

	//send back how many received vertices are in the boundary
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

	//send back the boundary flags
        Teuchos::ArrayView<int> is_bndry_send_view = Teuchos::arrayViewFromVector(is_bndry_send);
        Teuchos::ArrayRCP<int> is_bndry_recv;
        std::vector<int> bndry_recvcounts(comm->getSize(),0);
        Teuchos::ArrayView<int> bndry_recvcounts_view = Teuchos::arrayViewFromVector(bndry_recvcounts);
        Zoltan2::AlltoAllv<int> (*comm, *env, is_bndry_send_view, recvcounts_view, is_bndry_recv, bndry_recvcounts_view);

	//add together the flags to compute how many boundary vertices are in the second ghost layer
        int boundaryverts = 0;
        for(int i = 0; i < is_bndry_recv.size(); i++){
          boundaryverts += is_bndry_recv[i];
        }
	//this cout is guarded by a verbose check.
        std::cout<<comm->getRank()<<": "<<boundaryverts<<" boundary verts out of "<<n2Ghosts<<" verts in 2GL\n";
      }


      //local CSR representation for the owned vertices:
      Teuchos::ArrayView<const lno_t> local_adjs_view = Teuchos::arrayViewFromVector(local_adjs);
      //NOTE: the offset ArrayView was constructed earlier, and is up-to-date.
      //
      //first ghost layer CSR representation:
      Teuchos::ArrayView<const offset_t> ghost_offsets = Teuchos::arrayViewFromVector(first_layer_ghost_offsets);
      Teuchos::ArrayView<const lno_t> ghost_adjacencies = Teuchos::arrayViewFromVector(local_ghost_adjs);
      Teuchos::ArrayView<const int> rand_view = Teuchos::arrayViewFromVector(rand);
      Teuchos::ArrayView<const gno_t> gid_view = Teuchos::arrayViewFromVector(ownedPlusGhosts);
      //these ArrayViews contain LIDs that are ghosted remotely,
      //and the Process IDs that have those ghost copies.
      //An LID may appear in the exportLIDs more than once.
      Teuchos::ArrayView<const lno_t> exportLIDs = importer->getExportLIDs();
      Teuchos::ArrayView<const int> exportPIDs = importer->getExportPIDs();

      //construct a quick lookup datastructure to map from LID to a
      //list of PIDs we need to send data to.
      std::unordered_map<lno_t, std::vector<int>> procs_to_send;
      for(int i = 0; i < exportLIDs.size(); i++){
        procs_to_send[exportLIDs[i]].push_back(exportPIDs[i]);
      }

      //call the coloring algorithm
      twoGhostLayer(nVtx, nVtx+nGhosts, local_adjs_view, offsets, ghost_adjacencies, ghost_offsets,
                    femv, gid_view, rand_view, owners2, mapWithCopies, procs_to_send);

      //copy colors to the output array
      ArrayRCP<int> colors = solution->getColorsRCP();
      auto femvdata = femv->getData(0);
      for(size_t i=0; i<nVtx; i++){
         colors[i] = femvdata[i];
      }

    }
//  private:

    //This is the coloring logic for all 2GL-based methods.
    //
    //INPUT ARGS:
    //
    //  n_local: the number of local owned vertices
    //
    //  n_total: the number of local owned vertices plus first layer ghosts
    //
    //  adjs: the CSR adjacencies for the graph, only including owned vertices
    //        and first layer ghosts
    //
    //  offsets: the CSR offsets for the graph, has owned offsets into adjacencies
    //
    //  ghost_adjs: the adjacency information for first-layer ghost vertices.
    //              Includes all neighbors (owned, first-layer ghost,
    //              second-layer ghost) for each first-layer ghost.
    //
    //
    //  ghost_offsets: the offsets into ghost_adjs, first layer ghost LIDs
    //                 minus n_local are used to index this
    //                 datastructure. I.E. ghost_offsets(0)
    //                 corresponds to the ghost with LID n_local
    //
    //  gids: a vector that holds GIDs for all vertices on this process
    //        (includes owned, and all ghosts, indexable by LID)
    //        The GIDs are ordered like this:
    //        ----------------------------------
    //        |         | first     | second   |
    //        | owned   | layer     | layer    |
    //        |         | ghosts    | ghosts   |
    //        ----------------------------------
    //                  ^           ^
    //                  n_local     n_total
    //
    //  rand: a vector that holds random numbers generated on GID for tie
    //        breaking. Indexed by LID.
    //
    //  ghost_owners: contains the process ID for the owner of each remote vertex.
    //          Indexable by LID. owners[i] = the owning process for vertex
    //          with GID = gids[i+n_local]
    //
    //  mapOwnedPlusGhosts: map that can translate between GID and LID
    //
    //  procs_to_send: for each vertex that is copied on a remote process,
    //                 this map contains the list of processes that have
    //                 a copy of a given vertex. Input LID, get the list
    //                 of PIDs that have a ghost copy of that LID.
    //
    //OUTPUT ARGS:
    //
    //  femv: an FEMultiVector that holds a dual view of the colors.
    //        After this call, femv contains updated colors.
    //
    void twoGhostLayer(const size_t n_local, const size_t n_total,
                       const Teuchos::ArrayView<const lno_t>& adjs,
                       const Teuchos::ArrayView<const offset_t>& offsets,
                       const Teuchos::ArrayView<const lno_t>& ghost_adjs,
                       const Teuchos::ArrayView<const offset_t>& ghost_offsets,
                       const Teuchos::RCP<femv_t>& femv,
                       const Teuchos::ArrayView<const gno_t>& gids,
                       const Teuchos::ArrayView<const int>& rand,
                       const Teuchos::ArrayView<const int>& ghost_owners,
                       RCP<const map_t> mapOwnedPlusGhosts,
		       const std::unordered_map<lno_t, std::vector<int>>& procs_to_send){
      //initialize timing variables
      double total_time = 0.0;
      double interior_time = 0.0;
      double comm_time = 0.0;
      double comp_time = 0.0;
      double recoloring_time=0.0;
      double conflict_detection = 0.0;

      //Number of rounds we are saving statistics for
      //100 is a decent default. Reporting requires --verbose argument.
      const int numStatisticRecordingRounds = 100;

      //includes all ghosts, including the second layer.
      const size_t n_ghosts = rand.size() - n_local;


      //Get the degrees of all ghost vertices
      //This is necessary for recoloring based on degrees,
      //we need ghost vertex degrees for consistency.
      std::vector<int> deg_send_cnts(comm->getSize(),0);
      std::vector<gno_t> deg_sdispls(comm->getSize()+1,0);
      for(int i = 0; i < ghost_owners.size(); i++){
        deg_send_cnts[ghost_owners[i]]++;
      }
      deg_sdispls[0] = 0;
      gno_t deg_sendsize = 0;
      std::vector<int> deg_sentcount(comm->getSize(),0);
      for(int i = 1; i < comm->getSize()+1; i++){
        deg_sdispls[i] = deg_sdispls[i-1] + deg_send_cnts[i-1];
	deg_sendsize += deg_send_cnts[i-1];
      }
      std::vector<gno_t> deg_sendbuf(deg_sendsize,0);
      for(int i = 0; i < ghost_owners.size(); i++){
        size_t idx = deg_sdispls[ghost_owners[i]] + deg_sentcount[ghost_owners[i]];
	deg_sentcount[ghost_owners[i]]++;
	deg_sendbuf[idx] = mapOwnedPlusGhosts->getGlobalElement(i+n_local);
      }
      Teuchos::ArrayView<int> deg_send_cnts_view = Teuchos::arrayViewFromVector(deg_send_cnts);
      Teuchos::ArrayView<gno_t> deg_sendbuf_view = Teuchos::arrayViewFromVector(deg_sendbuf);
      Teuchos::ArrayRCP<gno_t> deg_recvbuf;
      std::vector<int> deg_recvcnts(comm->getSize(),0);
      Teuchos::ArrayView<int> deg_recvcnts_view = Teuchos::arrayViewFromVector(deg_recvcnts);
      Zoltan2::AlltoAllv<gno_t>(*comm, *env, deg_sendbuf_view, deg_send_cnts_view, deg_recvbuf, deg_recvcnts_view);

      //replace GID with the degree the owning process has for this vertex.
      //(this will include all edges present for this vertex globally)
      //This is safe to do, since we sent ghosts to their owners.
      for(int i = 0; i < deg_recvbuf.size(); i++){
        lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_recvbuf[i]);
	deg_recvbuf[i] = offsets[lid+1] - offsets[lid];
      }
      //send modified buffer back
      ArrayRCP<gno_t> ghost_degrees;
      Zoltan2::AlltoAllv<gno_t>(*comm, *env, deg_recvbuf(), deg_recvcnts_view, ghost_degrees, deg_send_cnts_view);

      //create the ghost degree views, for use on host and device.
      Kokkos::View<gno_t*, device_type> ghost_degrees_dev("ghost degree view",ghost_degrees.size());
      typename Kokkos::View<gno_t*, device_type>::HostMirror ghost_degrees_host = Kokkos::create_mirror(ghost_degrees_dev);
      for(int i = 0; i < ghost_degrees.size(); i++){
        lno_t lid = mapOwnedPlusGhosts->getLocalElement(deg_sendbuf[i]);
	ghost_degrees_host(lid-n_local) = ghost_degrees[i];
      }
      Kokkos::deep_copy(ghost_degrees_dev, ghost_degrees_host);

      //track the size of sends and receives through coloring rounds.
      gno_t recvPerRound[numStatisticRecordingRounds];
      gno_t sentPerRound[numStatisticRecordingRounds];

      //find global max degree to determine which
      //coloring algorithm will be the most efficient
      //(For D1-2GL this is important, D2 and PD2 should always use NBBIT
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

      //the initial coloring is able to use a standard CSR representation, so we construct that here.
      //This is a direct translation of the offsets and adjs arguments into Kokkos::Views.
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

      //in order to correctly recolor, all ghost vertices on this process need an entry in the CSR offsets.
      //Otherwise, the color of the ghosts will be ignored, and the coloring will not converge.
      Kokkos::View<offset_t*, device_type> dist_degrees_dev("Owned+Ghost degree view",rand.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_degrees_host = Kokkos::create_mirror(dist_degrees_dev);

      //calculate the local degrees for the owned, first layer ghosts, and second layer ghosts
      //to be used to construct the CSR representation of the local graph.
      //owned
      for(Teuchos_Ordinal i = 0; i < offsets.size()-1; i++) dist_degrees_host(i) = offsets[i+1] - offsets[i];
      //first layer ghosts
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++) dist_degrees_host(i+n_local) = ghost_offsets[i+1] - ghost_offsets[i];
      //second layer ghosts
      for(Teuchos_Ordinal i = 0; i < ghost_adjs.size(); i++){
	//second layer ghosts have LID >= n_total
	if((size_t)ghost_adjs[i] >= n_total ){
          dist_degrees_host(ghost_adjs[i])++;
	}
      }

      //The structure of the final CSR constructed here looks like:
      //
      //                       Index by LID--v
      //                       --------------------------------------------
      //                       |               | first       |second      |
      // dist_offsets_dev/host |owned verts    | layer       |layer       |
      //                       |               | ghosts      |ghosts      |
      //                       --------------------------------------------
      //                                             |indexes
      //                                             v
      //                       --------------------------------------------
      //                       |               |first        |second      |
      // dist_adjs_dev/host    |owned adjs     |layer        |layer       |
      //                       |               |ghost adjs   |ghost adjs  |
      //                       --------------------------------------------
      //
      // We end up with a CSR that has adjacency information for all the vertices on
      // this process. We include only edges on the process, so ghosts have only partial
      // adjacency information.
      //
      // We symmetrize the local graph representation as well, for
      // KokkosKernels to behave as we need it to. This means that edges to
      // second layer ghosts must be symmetrized, so we end up with offsets
      // that correspond to second layer ghosts.
      Kokkos::View<offset_t*, device_type> dist_offsets_dev("Owned+Ghost Offset view", rand.size()+1);
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host = Kokkos::create_mirror(dist_offsets_dev);

      //we can construct the entire offsets using the degrees we computed above
      dist_offsets_host(0) = 0;
      offset_t total_adjs = 0;
      for(Teuchos_Ordinal i = 1; i < rand.size()+1; i++){
        dist_offsets_host(i) = dist_degrees_host(i-1) + dist_offsets_host(i-1);
        total_adjs += dist_degrees_host(i-1);
      }
      Kokkos::View<lno_t*, device_type> dist_adjs_dev("Owned+Ghost adjacency view", total_adjs);
      typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host = Kokkos::create_mirror(dist_adjs_dev);

      //zero out the degrees to use it as a counter.
      //The offsets now allow us to compute degrees,
      //so we aren't losing anything.
      for(Teuchos_Ordinal i = 0; i < rand.size(); i++){
        dist_degrees_host(i) = 0;
      }
      //We can just copy the adjacencies for the owned verts and first layer ghosts
      for(Teuchos_Ordinal i = 0; i < adjs.size(); i++) dist_adjs_host(i) = adjs[i];
      for(Teuchos_Ordinal i = adjs.size(); i < adjs.size() + ghost_adjs.size(); i++) dist_adjs_host(i) = ghost_adjs[i-adjs.size()];

      //We have to symmetrize edges that involve a second layer ghost.
      //Add edges from the second layer ghosts to their neighbors.
      for(Teuchos_Ordinal i = 0; i < ghost_offsets.size()-1; i++){
	//loop through each first layer ghost
        for(offset_t j = ghost_offsets[i]; j < ghost_offsets[i+1]; j++){
          //if an adjacency is a second layer ghost
          if((size_t)ghost_adjs[j] >= n_total){
            //compute the offset to place the symmetric edge, and place it.
            dist_adjs_host(dist_offsets_host(ghost_adjs[j]) + dist_degrees_host(ghost_adjs[j])) = i + n_local;
	    //increment the counter to keep track of how many adjacencies
	    //have been placed.
            dist_degrees_host(ghost_adjs[j])++;
          }
        }
      }
      //copy the host views back to the device views
      Kokkos::deep_copy(dist_degrees_dev,dist_degrees_host);
      Kokkos::deep_copy(dist_offsets_dev,dist_offsets_host);
      Kokkos::deep_copy(dist_adjs_dev, dist_adjs_host);

      //this view represents how many conflicts were found
      Kokkos::View<size_t*, device_type> recoloringSize("Recoloring Queue Size",1);
      typename Kokkos::View<size_t*, device_type>::HostMirror recoloringSize_host = Kokkos::create_mirror(recoloringSize);
      recoloringSize_host(0) = 0;
      Kokkos::deep_copy(recoloringSize, recoloringSize_host);

      //create a view for the tie-breaking random numbers.
      if(verbose) std::cout<<comm->getRank()<<": constructing rand and GIDs views\n";
      Kokkos::View<int*, device_type> rand_dev("Random View", rand.size());
      typename Kokkos::View<int*, device_type>::HostMirror rand_host = Kokkos::create_mirror(rand_dev);
      for(Teuchos_Ordinal i = 0; i < rand.size(); i ++) rand_host(i) = rand[i];
      Kokkos::deep_copy(rand_dev,rand_host);

      //create a view for the global IDs of each vertex.
      Kokkos::View<gno_t*, device_type> gid_dev("GIDs", gids.size());
      typename Kokkos::View<gno_t*, device_type>::HostMirror gid_host = Kokkos::create_mirror(gid_dev);
      for(Teuchos_Ordinal i = 0; i < gids.size(); i++) gid_host(i) = gids[i];
      Kokkos::deep_copy(gid_dev,gid_host);

      //These variables will be initialized by constructBoundary()
      //
      //boundary_verts_dev holds all owned vertices that could possibly conflict
      //with a remote vertex. Stores LIDs.
      Kokkos::View<lno_t*, device_type> boundary_verts_dev;
      //this is the number of vertices in the boundary.
      if(verbose) std::cout<<comm->getRank()<<": constructing communication and recoloring lists\n";

      //We keep track of the vertices that need to get recolored
      //this list can contain ghost vertices, so it needs to be initialized
      //to the number of all vertices on the local process.
      //rand has an entry for each vertex, so its size corresponds to what is needed.
      Kokkos::View<lno_t*, device_type> verts_to_recolor_view("verts to recolor", rand.size());
      typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_recolor_host = create_mirror(verts_to_recolor_view);

      //This view keeps track of the size of the list of vertices to recolor.
      Kokkos::View<int*, device_type> verts_to_recolor_size("verts to recolor size",1);
      Kokkos::View<int*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_size_atomic = verts_to_recolor_size;
      typename Kokkos::View<int*, device_type>::HostMirror verts_to_recolor_size_host = create_mirror(verts_to_recolor_size);

      //initialize the host view
      verts_to_recolor_size_host(0) = 0;
      //copy to device
      Kokkos::deep_copy(verts_to_recolor_size, verts_to_recolor_size_host);

      //verts to send can only include locally owned vertices,
      //so we can safely allocate the view to size n_local.
      //
      //This view is a list of local vertices that are going to be
      //recolored, and need to be sent to their remote copies.
      Kokkos::View<lno_t*, device_type> verts_to_send_view("verts to send", n_local);
      typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send_host = create_mirror(verts_to_send_view);

      //this view keeps track of the size of verts_to_send.
      Kokkos::View<size_t*, device_type> verts_to_send_size("verts to send size",1);
      Kokkos::View<size_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic = verts_to_send_size;
      typename Kokkos::View<size_t*, device_type>::HostMirror verts_to_send_size_host = create_mirror(verts_to_send_size);

      verts_to_send_size_host(0) = 0;
      Kokkos::deep_copy(verts_to_send_size, verts_to_send_size_host);

      if(verbose) std::cout<<comm->getRank()<<": Constructing the boundary\n";

      //this function allocates and initializes the boundary_verts_dev view,
      //sets the boundary_size accordingly, and adds vertices to the
      //verts_to_send_atomic view, updating the size view as well.
      constructBoundary(n_local, dist_offsets_dev, dist_adjs_dev, dist_offsets_host, dist_adjs_host, boundary_verts_dev,
	                verts_to_send_view, verts_to_send_size_atomic);

      //this boolean chooses which KokkosKernels algorithm to use.
      //It was experimentally chosen for distance-1 coloring.
      //Should be disregarded for distance-2 colorings.
      bool use_vbbit = (global_max_degree < 6000);
      //Done initializing, start coloring!

      //use a barrier if we are reporting timing info
      if(timing) comm->barrier();
      interior_time = timer();
      total_time = timer();
      //give the entire local graph to KokkosKernels to color
      this->colorInterior(n_local, adjs_dev, offsets_dev, femv,adjs_dev,0,use_vbbit);
      interior_time = timer() - interior_time;
      comp_time = interior_time;

      //ghost_colors holds the colors of only ghost vertices.
      //ghost_colors(0) holds the color of a vertex with LID n_local.
      //To index this with LIDs, subtract n_local. If an LID is < n_local,
      //it will not have a color stored in this view.
      Kokkos::View<int*,device_type> ghost_colors("ghost color backups", n_ghosts);

      //communicate the initial coloring.
      if(verbose) std::cout<<comm->getRank()<<": communicating\n";

      //update the host views for the verts to send,
      //doOwnedToGhosts needs host memory views, but doesn't
      //change the values in them, so no need to copy afterwards
      Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
      Kokkos::deep_copy(verts_to_send_size_host, verts_to_send_size);
      gno_t recv,sent;
      comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,n_local,verts_to_send_host,verts_to_send_size_host,femv,procs_to_send,sent,recv);
      sentPerRound[0] = sent;
      recvPerRound[0] = recv;

      //store ghost colors so we can restore them after recoloring.
      //the local process can't color ghosts correctly, so we
      //reset the colors to avoid consistency issues.
      //get the color view from the FEMultiVector
      auto femvColors = femv->getLocalViewDevice(Tpetra::Access::ReadWrite);
      auto femv_colors = subview(femvColors, Kokkos::ALL, 0);
      Kokkos::parallel_for("get femv colors",
        Kokkos::RangePolicy<execution_space, int>(0,n_ghosts),
        KOKKOS_LAMBDA(const int& i){
          ghost_colors(i) = femv_colors(i+n_local);
        });
      Kokkos::fence();

      double temp = timer();
      //detect conflicts only for ghost vertices
      bool recolor_degrees = this->pl->template get<bool>("recolor_degrees",false);
      if(verbose) std::cout<<comm->getRank()<<": detecting conflicts\n";

      //these sizes will be updated by detectConflicts, zero them out beforehand
      verts_to_send_size_host(0) = 0;
      verts_to_recolor_size_host(0) = 0;
      recoloringSize_host(0) = 0;
      Kokkos::deep_copy(verts_to_send_size, verts_to_send_size_host);
      Kokkos::deep_copy(verts_to_recolor_size, verts_to_recolor_size_host);
      Kokkos::deep_copy(recoloringSize, recoloringSize_host);

      detectConflicts(n_local, dist_offsets_dev, dist_adjs_dev, femv_colors, boundary_verts_dev,
	       verts_to_recolor_view, verts_to_recolor_size_atomic, verts_to_send_view, verts_to_send_size_atomic,
	       recoloringSize, rand_dev, gid_dev, ghost_degrees_dev, recolor_degrees);

      //these sizes were updated by detectConflicts,
      //copy the device views back into host memory
      Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
      Kokkos::deep_copy(verts_to_send_size_host, verts_to_send_size);
      Kokkos::deep_copy(recoloringSize_host, recoloringSize);
      Kokkos::deep_copy(verts_to_recolor_size_host, verts_to_recolor_size);

      if(comm->getSize() > 1){
        conflict_detection = timer() - temp;
	comp_time += conflict_detection;
      }
      //all conflicts detected!
      if(verbose) std::cout<<comm->getRank()<<": starting to recolor\n";
      //variables for statistics
      double totalPerRound[numStatisticRecordingRounds];
      double commPerRound[numStatisticRecordingRounds];
      double compPerRound[numStatisticRecordingRounds];
      double recoloringPerRound[numStatisticRecordingRounds];
      double conflictDetectionPerRound[numStatisticRecordingRounds];
      uint64_t vertsPerRound[numStatisticRecordingRounds];
      uint64_t incorrectGhostsPerRound[numStatisticRecordingRounds];
      int distributedRounds = 1;
      totalPerRound[0] = interior_time + comm_time + conflict_detection;
      recoloringPerRound[0] = 0;
      commPerRound[0] = comm_time;
      compPerRound[0] = interior_time + conflict_detection;
      conflictDetectionPerRound[0] = conflict_detection;
      recoloringPerRound[0] = 0;
      vertsPerRound[0] = 0;
      incorrectGhostsPerRound[0]=0;
      typename Kokkos::View<int*, device_type>::HostMirror ghost_colors_host;
      typename Kokkos::View<lno_t*, device_type>::HostMirror boundary_verts_host;
      size_t serial_threshold = this->pl->template get<int>("serial_threshold",0);
      //see if recoloring is necessary.
      size_t totalConflicts = 0;
      size_t localConflicts = recoloringSize_host(0);
      Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_SUM, 1, &localConflicts, &totalConflicts);
      bool done = !totalConflicts;
      if(comm->getSize()==1) done = true;

      //recolor until no conflicts are left
      while(!done){
	//if the number of vertices left to recolor is less than
	//the serial threshold set by the user, finish the coloring
	//only using host views in a serial execution space.
	if(recoloringSize_host(0) < serial_threshold) break;
        if(distributedRounds < numStatisticRecordingRounds) {
          vertsPerRound[distributedRounds] = verts_to_recolor_size_host(0);
        }

        if(timing) comm->barrier();
        double recolor_temp = timer();
        //recolor using KokkosKernels' coloring function
        if(verts_to_recolor_size_host(0) > 0){
	  this->colorInterior(femv_colors.size(), dist_adjs_dev, dist_offsets_dev,femv,verts_to_recolor_view,verts_to_recolor_size_host(0),use_vbbit);
	}

	if(distributedRounds < numStatisticRecordingRounds){
          recoloringPerRound[distributedRounds] = timer() - recolor_temp;
          recoloring_time += recoloringPerRound[distributedRounds];
          comp_time += recoloringPerRound[distributedRounds];
          compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
          totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	} else if(timing){
	  double recoloring_round_time = timer() - recolor_temp;
	  recoloring_time += recoloring_round_time;
	  comp_time += recoloring_round_time;
	}

	//reset the ghost colors to what they were before recoloring
	//the recoloring does not have enough information to color
	//ghosts correctly, so we set the colors to what they were before
	//to avoid consistency issues.
        Kokkos::parallel_for("set femv colors",
          Kokkos::RangePolicy<execution_space, int>(0,n_ghosts),
          KOKKOS_LAMBDA(const int& i){
            femv_colors(i+n_local) = ghost_colors(i);
          });
        Kokkos::fence();

	//send views are up-to-date, they were copied after conflict detection.
        //communicate the new colors

        // Reset device views
        femvColors = decltype(femvColors)();
        femv_colors = decltype(femv_colors)();
        double curr_comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,n_local,verts_to_send_host,verts_to_send_size_host,femv,procs_to_send,sent,recv);
	comm_time += curr_comm_time;

        if(distributedRounds < numStatisticRecordingRounds){
          commPerRound[distributedRounds] = curr_comm_time;
	  recvPerRound[distributedRounds] = recv;
          sentPerRound[distributedRounds] = sent;
          totalPerRound[distributedRounds] += commPerRound[distributedRounds];
	}

	//store ghost colors before we uncolor and recolor them.
	//this process doesn't have enough info to correctly color
	//ghosts, so we set them back to what they were before to
	//remove consistency issues.
        femvColors = femv->getLocalViewDevice(Tpetra::Access::ReadWrite);
        femv_colors = subview(femvColors, Kokkos::ALL, 0);
        Kokkos::parallel_for("get femv colors 2",
          Kokkos::RangePolicy<execution_space, int>(0,n_ghosts),
          KOKKOS_LAMBDA(const int& i){
            ghost_colors(i) = femv_colors(i+n_local);
          });
        Kokkos::fence();


        //these views will change in detectConflicts, they need
	//to be zeroed out beforehand
        verts_to_send_size_host(0) = 0;
	verts_to_recolor_size_host(0) = 0;
	recoloringSize_host(0) = 0;
	Kokkos::deep_copy(verts_to_send_size, verts_to_send_size_host);
	Kokkos::deep_copy(verts_to_recolor_size, verts_to_recolor_size_host);
	Kokkos::deep_copy(recoloringSize, recoloringSize_host);

        //check for further conflicts
        double detection_temp = timer();

        detectConflicts(n_local, dist_offsets_dev, dist_adjs_dev,femv_colors, boundary_verts_dev,
	         verts_to_recolor_view, verts_to_recolor_size_atomic, verts_to_send_view, verts_to_send_size_atomic,
		 recoloringSize, rand_dev, gid_dev, ghost_degrees_dev, recolor_degrees);

        //copy the updated device views back into host memory where necessary
        Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
        Kokkos::deep_copy(verts_to_send_size_host, verts_to_send_size);
	//we only use the list of verts to recolor on device, no need to copy to host.
        Kokkos::deep_copy(verts_to_recolor_size_host, verts_to_recolor_size);
	Kokkos::deep_copy(recoloringSize_host, recoloringSize);

	if(distributedRounds < numStatisticRecordingRounds){
	  conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
          conflict_detection += conflictDetectionPerRound[distributedRounds];
          compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
          totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
          comp_time += conflictDetectionPerRound[distributedRounds];
        } else if(timing){
	  double conflict_detection_round_time = timer() - detection_temp;
	  conflict_detection += conflict_detection_round_time;
	  comp_time += conflict_detection_round_time;
	}

        distributedRounds++;
        size_t localDone = recoloringSize_host(0);
        size_t globalDone = 0;
        Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
        done = !globalDone;

      }//end device coloring


      //If we reach this point and still have vertices to color,
      //the rest of the coloring is happening in serial on the host.
      //
      //First, we need to copy device-only views to host mirrors
      if(recoloringSize_host(0) > 0 || !done){
	ghost_colors_host = Kokkos::create_mirror_view_and_copy(host_mem(),ghost_colors,"ghost_colors host mirror");
	boundary_verts_host = Kokkos::create_mirror_view_and_copy(host_mem(),boundary_verts_dev,"boundary_verts host mirror");
      }

      //Now we do a similar coloring loop to before,
      //but using only host views in a serial execution space.
      // Reset device views
      femvColors = decltype(femvColors)();
      femv_colors = decltype(femv_colors)();
      while(recoloringSize_host(0) > 0 || !done){
	auto femvColors_host = femv->getLocalViewHost(Tpetra::Access::ReadWrite);
	auto colors_host = subview(femvColors_host, Kokkos::ALL, 0);
	if(distributedRounds < numStatisticRecordingRounds){
	  vertsPerRound[distributedRounds] = recoloringSize_host(0);
	}
	if(verbose) std::cout<<comm->getRank()<<": starting to recolor, serial\n";
        if(timing) comm->barrier();

	double recolor_temp = timer();
	if(verts_to_recolor_size_host(0) > 0){
	  this->colorInterior_serial(colors_host.size(), dist_adjs_host, dist_offsets_host, femv,
			             verts_to_recolor_host, verts_to_recolor_size_host(0), true);
	}
	if(distributedRounds < numStatisticRecordingRounds){
	  recoloringPerRound[distributedRounds] = timer() - recolor_temp;
	  recoloring_time += recoloringPerRound[distributedRounds];
	  comp_time += recoloringPerRound[distributedRounds];
	  compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	  totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	} else if(timing){
	  double recoloring_serial_round_time = timer() - recolor_temp;
	  recoloring_time += recoloring_serial_round_time;
	  comp_time += recoloring_serial_round_time;
	}

	//reset the ghost colors to their previous values to avoid
	//consistency issues.
	for(size_t i = 0; i < n_ghosts; i++){
	  colors_host(i+n_local) = ghost_colors_host(i);
	}

	double curr_comm_time = doOwnedToGhosts(mapOwnedPlusGhosts, n_local,verts_to_send_host, verts_to_send_size_host, femv, procs_to_send, sent,recv);
	comm_time += curr_comm_time;

	if(distributedRounds < numStatisticRecordingRounds){
	  commPerRound[distributedRounds] = curr_comm_time;
	  recvPerRound[distributedRounds] = recv;
	  sentPerRound[distributedRounds] = sent;
	  totalPerRound[distributedRounds] += commPerRound[distributedRounds];
	}

	//store updated ghost colors we received from their owners
	//before conflict detection and recoloring changes them locally.
	for(size_t i = 0; i < n_ghosts; i++){
	  ghost_colors_host(i) = colors_host(i+n_local);
	}

	if(timing) comm->barrier();
	double detection_temp = timer();

	//zero these out, they'll be updated by detectConflicts_serial
        verts_to_recolor_size_host(0) = 0;
	verts_to_send_size_host(0) = 0;
	recoloringSize_host(0) = 0;

        detectConflicts_serial(n_local,dist_offsets_host, dist_adjs_host, colors_host, boundary_verts_host,
	                       verts_to_recolor_host, verts_to_recolor_size_host, verts_to_send_host, verts_to_send_size_host,
		               recoloringSize_host, rand_host, gid_host, ghost_degrees_host, recolor_degrees);

	//no need to update the host views, we're only using host views now.

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

	size_t globalDone = 0;
	size_t localDone = recoloringSize_host(0);
	Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_SUM, 1, &localDone, &globalDone);
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

        //if(comm->getRank() == 0) printf("did %d rounds of distributed coloring\n", distributedRounds);
        uint64_t totalVertsPerRound[numStatisticRecordingRounds];
        uint64_t totalBoundarySize = 0;
        uint64_t totalIncorrectGhostsPerRound[numStatisticRecordingRounds];
        double finalTotalPerRound[numStatisticRecordingRounds];
        double maxRecoloringPerRound[numStatisticRecordingRounds];
        double minRecoloringPerRound[numStatisticRecordingRounds];
        double finalCommPerRound[numStatisticRecordingRounds];
        double finalCompPerRound[numStatisticRecordingRounds];
        double finalConflictDetectionPerRound[numStatisticRecordingRounds];
        gno_t finalRecvPerRound[numStatisticRecordingRounds];
        gno_t finalSentPerRound[numStatisticRecordingRounds];
        for(int i = 0; i < numStatisticRecordingRounds; i++){
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
        Teuchos::reduceAll<int,uint64_t>(*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,vertsPerRound,totalVertsPerRound);
        Teuchos::reduceAll<int,uint64_t>(*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,incorrectGhostsPerRound,totalIncorrectGhostsPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,totalPerRound, finalTotalPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,recoloringPerRound,maxRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,numStatisticRecordingRounds,recoloringPerRound,minRecoloringPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,commPerRound,finalCommPerRound);
        Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,numStatisticRecordingRounds,compPerRound,finalCompPerRound);
        Teuchos::reduceAll<int,double>(*comm,
			               Teuchos::REDUCE_MAX,numStatisticRecordingRounds,conflictDetectionPerRound,finalConflictDetectionPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,recvPerRound,finalRecvPerRound);
        Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,numStatisticRecordingRounds,sentPerRound,finalSentPerRound);
        std::cout << "Rank " << comm->getRank()
                  << ": boundary size: " << localBoundaryVertices << std::endl;
        if(comm->getRank() == 0)
          std::cout << "Total boundary size: " << totalBoundarySize << std::endl;
        for(int i = 0; i < std::min((int)distributedRounds,numStatisticRecordingRounds); i++){
          std::cout << "Rank " << comm->getRank()
                    << ": recolor " << vertsPerRound[i]
                    << " vertices in round " << i << std::endl;
          std::cout << "Rank " << comm->getRank()
                    << " sentbuf had " << sentPerRound[i]
                    << " entries in round " << i << std::endl;
          if(comm->getRank()==0){
            std::cout << "recolored " << totalVertsPerRound[i]
                      << " vertices in round " << i << std::endl;
            std::cout << totalIncorrectGhostsPerRound[i]
                      << " inconsistent ghosts in round " << i << std::endl;
            std::cout << "total time in round " << i
                      << ": " << finalTotalPerRound[i] << std::endl;
            std::cout << "recoloring time in round " << i
                      << ": " << maxRecoloringPerRound[i] << std::endl;
            std::cout << "min recoloring time in round " << i
                      << ": " << minRecoloringPerRound[i] << std::endl;
            std::cout << "conflict detection time in round " << i
                      << ": " << finalConflictDetectionPerRound[i] << std::endl;
            std::cout << "comm time in round " << i
                      << ": " << finalCommPerRound[i] << std::endl;
            std::cout << "recvbuf size in round " << i
                      << ": " << finalRecvPerRound[i] << std::endl;
            std::cout << "sendbuf size in round " << i
                      << ": " << finalSentPerRound[i] << std::endl;
            std::cout << "comp time in round " << i
                      << ": " << finalCompPerRound[i] << std::endl;
          }
        }
      } else if (timing){
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
          std::cout << "Total Time: " << global_total_time << std::endl;
          std::cout << "Interior Time: " << global_interior_time << std::endl;
          std::cout << "Recoloring Time: " << global_recoloring_time << std::endl;
          std::cout << "Min Recoloring Time: " << global_min_recoloring_time << std::endl;
          std::cout << "Conflict Detection Time: " << global_conflict_detection << std::endl;
          std::cout << "Comm Time: " << global_comm_time << std::endl;
          std::cout << "Comp Time: " << global_comp_time << std::endl;
        }
      }
    }
}; //end class

}//end namespace Zoltan2

#endif
