//@HEADER
#ifndef _ZOLTAN2_ICESHEET_HPP_
#define _ZOLTAN2_ICESHEET_HPP_
#include<unordered_set>
#include<utility>
#include<vector>
#include<Zoltan2_GraphModel.hpp>
#include<Zoltan2_Algorithm.hpp>
#include<Zoltan2_PartitioningSolution.hpp>
#include<Zoltan2_Util.hpp>
#include<Zoltan2_TPLTraits.hpp>
#include<Zoltan2_IceProp.hpp>
#include<Zoltan2_AlltoAll.hpp>

/*! \file Zoltan2_IceSheet.hpp
 * \brief The user interface for ice-sheet connectivity functionality.
 *
 * This file contains the entry point for the ice sheet connectivity
 * label propagation algorithm. 
*/

namespace Zoltan2{
  /*! Utility function for detecting degenerate features of ice sheets
   *  
   *  \param (input) problemComm A Teuchos::RCP<const Teuchos::Comm<int>> 
   *                 representing the communicator for the problem.
   *
   *  \param (input) adapter A Teuchos::RCP<const GraphAdapter<>> 
   *                 representing the structure of the ice sheet
   *                 mesh's bottom layer.
   *
   *  \param (input) basalFriction An ArrayView of boolean flags for each 
   *                 local vertex with true representing a grounded vertex 
   *                 and false representing a floating vertex.
   *                 These indicies must be consistent with the local vertex 
   *                 identifiers for the input graph.
   *
   *  \param (input) boundary_edges An ArrayView of global vertex identifiers 
   *                 representing all of the edges that the current process 
   *                 knows about which are on the boundary between the ice and 
   *                 the water. The IDs must be present in the overlap map of 
   *                 the current process if the ID is not owned by this process.
   *                 Edges should not be specified more than once, locally. 
   *                 Edges may be specified more than once globally.
   *
   *  \param (output) vertex_status is an output argument, and must have a size 
   *                  corresponding to the number of local vertices. After 
   *                  calling DetectDegenerateVertices, this ArrayView will be 
   *                  filled with an integer indicating the corresponding 
   *                  vertex's status. 
   *                        Zoltan2::IceGrounded (-2) indicates a vertex has 
   *                        a sufficient connection to ground, and should be
   *                        kept.
   *
   *                        Zoltan2::IceFloating (-1) indicates a vertex is 
   *                        floating and should be removed. 
   *
   *                        Zoltan2::IceHinged   (0) indicates a vertex is 
   *                        part of a floating hinge, and should be removed.
   *
   *  \param (output) hinge_vertices is an output argument, and must have a 
   *                  size corresponding to the number of local vertices.
   *                  For a vertex with a vertex_status of 
   *                  Zoltan2::IceHinged(0), the entry in this array 
   *                  corresponds to the hinge vertex. If multiple hinges 
   *                  are chained together and removed, this entry will 
   *                  indicate the first hinge vertex that has a vertex_status
   *                  of Zoltan2::IceGrounded (-2). Note: all of the entries 
   *                  contained in this ArrayView are valid global vertex 
   *                  identifiers. It is necessary to check the vertex_status 
   *                  before using this output.
   *
   */
  template<typename Adapter>
  void DetectDegenerateVertices( 
       const RCP<const Comm<int>> &problemComm,
       const Adapter &adapter,
       const Teuchos::ArrayView<const bool> &basalFriction,
       const Teuchos::ArrayView<const typename Adapter::gno_t> &boundary_edges,
       const Teuchos::ArrayView<IcePropVtxStatus> &vertex_status,
       const Teuchos::ArrayView<typename Adapter::gno_t> &hinge_vertices){

    typedef typename Adapter::base_adapter_t base_adapter_t;
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;
    typedef typename Adapter::part_t part_t;
    typedef typename Adapter::user_t user_t;
    typedef typename Adapter::userCoord_t userCoord_t;
    typedef Tpetra::Map<lno_t, gno_t> map_t;
    
    RCP<const Environment> env = rcp(new Environment(problemComm));
    RCP<const base_adapter_t> b_adapter = rcpFromRef(adapter);
    modelFlag_t flags;
    flags.reset();
    
    flags.set(REMOVE_SELF_EDGES);
  
    env->debug(DETAILED_STATUS, "	building graph model");
    //const RCP<const GraphAdapter<user_t>> g_adapter = adapter;
    RCP<const GraphModel<base_adapter_t>> model = rcp(new 
                                      GraphModel<base_adapter_t>(b_adapter,env,
			                                  problemComm, flags));
    env->debug(DETAILED_STATUS, "	graph model built");
    
    //Get vertex GIDs, in a locally indexed array
    ArrayView<const gno_t> vtxIDs;
    ArrayView<StridedData<lno_t, scalar_t> > vwgts;
    size_t nVtx = model->getVertexList(vtxIDs, vwgts);
    //at this point, weights are not used.
 
    //Get the edge information
    ArrayView<const gno_t> adjs; 
    ArrayView<const offset_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > ewgts;
    size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);
    //edge weights are not used either.
  
    //copy data over to regular arrays
    gno_t* out_edges = NULL;
    typename map_t::local_ordinal_type* out_offsets = NULL;
    gno_t* global_ids = NULL;
    TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
    TPL_Traits<typename map_t::local_ordinal_type, 
               const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);
    TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&global_ids, vtxIDs);
  
    //get the rank of the current process
    int me = problemComm->getRank();
  
    //create the Tpetra map for the global indices
    Tpetra::global_size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const map_t > map = rcp(new map_t(dummy, vtxIDs, 0, problemComm)); 
  
    //sentinel value to detect unsuccessful searches
    lno_t fail = Teuchos::OrdinalTraits<lno_t>::invalid();

    //map to make sure ghosted vertices aren't duplicated in the Tpetra Maps
    std::unordered_map<gno_t,lno_t> ghost_map;
    //The new local IDs for each ghosted vertex
    std::vector<lno_t> ghostGIDs;
    int nGhosts = 0;
    //count how many ghosts are in the edge list, to create the mapWithCopies.
    for(size_t i = 0; i < nEdge; i++){
      if(map->getLocalElement(adjs[i]) == fail){
        if(ghost_map.count(adjs[i]) == 0){
          ghost_map[adjs[i]] = nVtx + nGhosts;
          ghostGIDs.push_back(adjs[i]);
          nGhosts++;
        }
      }
    }
  
    //build nVtx + nGhosts size array of grounding flags
    bool* grounding = new bool[nVtx];
    for(size_t i = 0; i < nVtx; i++){
        grounding[i] = basalFriction[i];
    }

  
    //use the count of ghosts + owned to make mapWithCopies, need to create a Teuchos array first.
    Teuchos::Array<gno_t> gids(nVtx+nGhosts);
    for(size_t i = 0; i < nVtx; i++){
      gids[i] = vtxIDs[i];
    }
    for(size_t i = nVtx; i < nVtx+nGhosts; i++){
      gids[i] = ghostGIDs[i-nVtx];
    }
  
    //create the Tpetra map with copies
    Tpetra::global_size_t dummy_2 = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const map_t > mapWithCopies = rcp(new map_t(dummy_2, gids, 0, 
                                                    problemComm));
 
    //communicate boundary edges back to owning processors
    gno_t num_boundary_edges = boundary_edges.size()/2;
    Teuchos::Array<int> owners(boundary_edges.size(),0);
  
    //get the owning process for each vertex in the boundary list
    mapWithCopies->getRemoteIndexList(boundary_edges,owners());
    //store the number of processes for future use
    int nprocs = problemComm->getSize();
  
    //send boundary edges to any remote processes that own an endpoint
    Teuchos::Array<int> sendcnts(nprocs,0);
    Teuchos::Array<int> recvcnts(nprocs,0);
    for(int i = 0; i < num_boundary_edges*2; i+= 2){
      if(owners[i] != -1){
        sendcnts[owners[i]] +=2;
      }
      if(owners[i+1] != -1 && owners[i] != owners[i+1]){
        sendcnts[owners[i+1]] +=2;
      }
    }
    //send the number of vertex IDs to be sent.
  
    //create sdispls and sentcount to build sendbuf
    Teuchos::Array<int> sdispls(nprocs+1, 0);
    for(int i = 1; i < nprocs+1; i++){
      sdispls[i] = sdispls[i-1] + sendcnts[i-1];
    }

    std::vector<int> sentcount;
    for(int i = 0; i < nprocs; i++){
      sentcount.push_back(0);
    }
    int sendsize = 0;
    int recvsize = 0;
    for(int i = 0; i < nprocs; i++){
      sendsize += sendcnts[i];
      recvsize += recvcnts[i];
    }
  
    Teuchos::Array<gno_t> sendbuf(sendsize,0);
    for(int i = 0; i > num_boundary_edges*2; i+=2){
      if(owners[i] != me){
        int proc_to_send = owners[i];
        int sendbufidx = sdispls[proc_to_send] + sentcount[proc_to_send];
        sentcount[proc_to_send] += 2;
        sendbuf[sendbufidx++] = boundary_edges[i];
        sendbuf[sendbufidx++] = boundary_edges[i+1];
      }
      if(owners[i+1] != me && owners[i] != owners[i+1]){
        int proc_to_send = owners[i];
        int sendbufidx = sdispls[proc_to_send] + sentcount[proc_to_send];
        sentcount[proc_to_send] += 2;
        sendbuf[sendbufidx++] = boundary_edges[i];
        sendbuf[sendbufidx++] = boundary_edges[i+1];
      }
    }
 
    //Do the final communication back to each remote process 
    Teuchos::ArrayRCP<gno_t> recvbuf;
    Zoltan2::AlltoAllv<gno_t>(*problemComm,*env,sendbuf(),sendcnts(),
                              recvbuf,recvcnts());

    //hash for the boundary-edge set
    struct pair_hash {
      inline std::size_t operator()(const std::pair<gno_t,gno_t>& v) const {
        return v.first * 10 + v.second;
      }
    };

    //create the set to see if we've counted this boundary edge before
    std::unordered_set<std::pair<gno_t,gno_t>,pair_hash> edge_set;

    //int* local_boundary_counts = new int[nVtx];
    //for(size_t i = 0; i < nVtx; i++){
    //  local_boundary_counts[i] = 0;
    //}
    Teuchos::Array<int> local_boundary_counts(nVtx+1,0);
    //insert the local boundary edges into the set, small global endpoint first.
    for(int i = 0; i < num_boundary_edges*2; i+=2){
      if(owners[i] == me || owners[i+1] == me){
        if(boundary_edges[i] < boundary_edges[i+1]){
          edge_set.insert(std::make_pair(boundary_edges[i],
                                         boundary_edges[i+1]));
        } else {
          edge_set.insert(std::make_pair(boundary_edges[i+1],
                                         boundary_edges[i]));
        }
      }
    }
    //factor in the locally-owned boundary edges
    for(int i = 0; i < num_boundary_edges*2; i++){
      if(map->getLocalElement(boundary_edges[i]) != fail){
        local_boundary_counts[map->getLocalElement(boundary_edges[i])]++;
      }
    }
    
    //use the received boundary edges, but check to make sure 
    //they haven't been used yet.
    for(int i = 0; i < recvsize; i+=2){
      //make sure the edge has not been used (if count returns 
      //1 for one of these, the edge has been used)
      if(edge_set.count(std::make_pair(recvbuf[i], recvbuf[i+1])) == 0 && 
         edge_set.count(std::make_pair(recvbuf[i+1],recvbuf[i])) == 0){
        //check which endpoints are owned locally, and set the boundary 
        //counts appropriately
        if(map->getLocalElement(recvbuf[i]) != fail){
          local_boundary_counts[map->getLocalElement(recvbuf[i])]++;
        }
        if(map->getLocalElement(recvbuf[i+1]) != fail){
          local_boundary_counts[map->getLocalElement(recvbuf[i+1])]++;
        }
        //add the received edge to the set, to ensure no received duplicates 
        //are counted
        if(recvbuf[i] < recvbuf[i+1]){
          edge_set.insert(std::make_pair(recvbuf[i],recvbuf[i+1]));
        } else {
  	  edge_set.insert(std::make_pair(recvbuf[i+1],recvbuf[i]));
        }
      }
    }
    //we should now have complete local knowledge of the boundary.

    //convert adjacency array to use local identifiers instead of global.
  
    Teuchos::Array<typename map_t::local_ordinal_type> out_edges_lid(nEdge,0);
    for(size_t i = 0; i < nEdge; i++){
      out_edges_lid[i] = mapWithCopies->getLocalElement(out_edges[i]);
    }
    icePropGraph<typename map_t::local_ordinal_type, 
                 offset_t> local_graph = {nVtx, nEdge, out_edges_lid,offsets};
    Zoltan2::iceSheetPropagation<typename map_t::local_ordinal_type, 
                                 typename map_t::global_ordinal_type,
                                 offset_t,
                                 map_t> prop(problemComm, map, 
                                                mapWithCopies, &local_graph, 
                                                local_boundary_counts, 
                                                grounding, nVtx, nGhosts);
    prop.propagate(vertex_status, hinge_vertices);
    
    delete [] grounding;
    
  }


}//end namespace Zoltan2

#endif// _ZOLTAN2_ICESHEET_HPP
