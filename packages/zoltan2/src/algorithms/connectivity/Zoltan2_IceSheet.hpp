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
#include<Zoltan2_VtxLabel.hpp>
#include<Zoltan2_AlltoAll.hpp>

namespace Zoltan2{
  template <typename Adapter>
  class IceProp : public Algorithm<Adapter> {
    public:
      typedef typename Adapter::base_adapter_t base_adapter_t;
      typedef typename Adapter::lno_t lno_t;
      typedef typename Adapter::gno_t gno_t;
      typedef typename Adapter::offset_t offset_t;
      typedef typename Adapter::scalar_t scalar_t;
      typedef typename Adapter::part_t part_t;
      typedef typename Adapter::user_t user_t;
      typedef typename Adapter::userCoord_t userCoord_t;
      typedef Tpetra::Map<lno_t, gno_t> map_t;
      //typedef map_t::global_ordinal_type mgno_t;
      
      //Arguments: problemComm is the communicator we will use for the propagation
      //           adapter is the graph adapter that represents the ice mesh's bottom layer
      //           basalFriction is the array of basalFriction values. Nonzero indicates grounded
      //           boundary_edges is the array of edges on the current processor that are on the
      //                          boundary, using global identifiers
      //           num_boundary_edges is the number of boundary edges there are
      IceProp(const RCP<const Comm<int> > &problemComm__,
	      const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__, 
	      bool* basalFriction, gno_t* boundary_edges__, int num_boundary_edges__):
         adapter(adapter__), grounding_flags(basalFriction), boundary_edges(boundary_edges__),
         num_boundary_edges(num_boundary_edges__), problemComm(problemComm__)
      {
	int me = problemComm->getRank();

	std::cout<<me<<": Creating new Environment\n";
	env = rcp(new Environment(problemComm));
	std::cout<<me<<": Created new Environment\n";
        modelFlag_t flags;
	flags.reset();

	std::cout<<me<<": Building the Model\n";
	///build the graph model from the GraphAdapter.	 
	buildModel(flags);
	std::cout<<me<<": Built the Model\n";
      }
      //This should probably return an RCP, but for now it's an int*
      //This will return a number of flags consistent with the 
      //number of locally owned vertices.
      int* getDegenerateFeatureFlags();
    private:
      const RCP<const base_adapter_t> adapter;
      void buildModel(modelFlag_t &flags);

      bool* grounding_flags;
      gno_t*  boundary_edges;
      int   num_boundary_edges;
      const RCP<const Comm<int> > problemComm;
      RCP<Environment> env;
      RCP<const GraphModel<base_adapter_t> > model;
  };
  
template <typename Adapter>
void IceProp<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
  //flags.set(GENERATE_CONSECUTIVE_IDS);
  
  this->env->debug(DETAILED_STATUS, "	building graph model");
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
			                           this->problemComm, flags));
  this->env->debug(DETAILED_STATUS, "	graph model built");
}
  

template <typename Adapter>
int* IceProp<Adapter>::getDegenerateFeatureFlags() {
  //using the GraphModel, construct 
  //mapOwned and mapWithCopies, as well as a 
  //csr graph, could just use the graph.h from 
  //the test directory. Should be a simple conversion.
  //Run the propagation, and return the flags.
  
  //get the number of local vertices and edges.
  //const size_t modelVerts = model->getLocalNumVertices();
  //const size_t modelEdges = model->getLocalNumEdges();

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

  gno_t* out_edges = NULL;
  typename map_t::local_ordinal_type* out_offsets = NULL;
  gno_t* global_ids = NULL;
  TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
  TPL_Traits<typename map_t::local_ordinal_type, const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);
  TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&global_ids, vtxIDs);
  
  
  
  int me = problemComm->getRank();

  std::cout<<me<<": Local vertices:\n";
  for(int i = 0; i < nVtx; i++){
    std::cout<<"\t"<<vtxIDs[i]<<"\n";
  }
  std::cout<<me<<": Local edges:\n";
  for(int i = 0; i < nVtx; i++){
    int outDegree = offsets[i+1] - offsets[i];
    for(int j = offsets[i]; j < offsets[i+1]; j++){
      std::cout<<me<<": "<<vtxIDs[i]<<" -- "<<adjs[j]<<"\n";
    }
  }
  
  //create the Tpetra map for the global indices
  Tpetra::global_size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  RCP<const map_t > map = rcp(new map_t(dummy, vtxIDs, 0, problemComm)); 
  
  //sentinel value to detect unsuccessful searches
  lno_t fail = Teuchos::OrdinalTraits<lno_t>::invalid();

  
  std::unordered_map<gno_t,lno_t> ghost_map;
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
  bool* grounding = new bool[nVtx+nGhosts];
  for(int i = 0; i < nVtx+nGhosts; i++){
    if(i < nVtx){
      grounding[i] = grounding_flags[i];
      if(grounding[i]){
        std::cout<<me<<": global vert "<<vtxIDs[i]<<" is initially grounded\n";
      }
    }
    else grounding[i] = false;
  }

  //std::cout<<me<<": number of ghosts = "<<nGhosts<<"\n";
  
  //use the count of ghosts + owned to make mapWithCopies, need to create a Teuchos array first.
  Teuchos::Array<gno_t> gids(nVtx+nGhosts);
  for(size_t i = 0; i < nVtx; i++){
    gids[i] = vtxIDs[i];
  }
  for(size_t i = nVtx; i < nVtx+nGhosts; i++){
    gids[i] = ghostGIDs[i-nVtx];
  }
  /*int ghostCount = 0;
  for(size_t i = 0; i < nEdge; i++){
    if(map->getLocalElement(adjs[i]) == fail){
      gids[nVtx+ghostCount] = adjs[i];
      ghostCount++;
    }
  }*/
  Tpetra::global_size_t dummy_2 = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  RCP<const map_t > mapWithCopies = rcp(new map_t(dummy_2, gids, 0, problemComm));
 
  //communicate boundary edges back to owning processors
  //
  //Teuchos::ArrayView<int> owners;
  //Teuchos::ArrayView<gno_t> boundary(num_boundary_edges*2);
  std::vector<int> owners_vec;
  std::vector<gno_t> boundary_vec;
  for(int i = 0; i < num_boundary_edges*2; i++){
    boundary_vec.push_back(boundary_edges[i]);
    owners_vec.push_back(0);
  }
  Teuchos::ArrayView<gno_t> boundary = Teuchos::arrayViewFromVector(boundary_vec);
  Teuchos::ArrayView<int> owners = Teuchos::arrayViewFromVector(owners_vec);
  mapWithCopies->getRemoteIndexList(boundary,owners);
  
  int nprocs = problemComm->getSize();
  std::cout<<me<<": nprocs = "<<nprocs<<"\n";
  int maxOwner = 0;
  int minOwner = 0;
  for(int i = 0; i < owners.size(); i++){
    if(owners[i] == -1) std::cout<<"Vertex "<<boundary[i]<<" does not have an owner\n";
    //if(owners[i] <minOwner) minOwner = owners[i];
    //if(owners[i] >maxOwner) maxOwner = owners[i];
  }
  std::cout<<me<<": maxOwner = "<<maxOwner<<", minOwner = "<<minOwner<<"\n";
  
  //send boundary edges to any remote processes that own an endpoint
  //COMMENT THIS BLOCK BACK IN
  std::vector<int> sendcnts;
  for(int i = 0; i < nprocs; i++) sendcnts.push_back(0);
  for(int i = 0; i < num_boundary_edges*2; i+= 2){
    if(owners[i] != -1){
      sendcnts[owners[i]] +=2;
    }
    if(owners[i+1] != -1 && owners[i] != owners[i+1]){
      sendcnts[owners[i+1]] +=2;
    }
  }
  std::vector<int> recvcnts;
  for(int i = 0; i < nprocs; i++) recvcnts.push_back(0);
  //send the number of vertex IDs to be sent.
  Teuchos::ArrayView<const int> sendCount = Teuchos::arrayViewFromVector(sendcnts);
  Teuchos::ArrayView<int> recvCount = Teuchos::arrayViewFromVector(recvcnts);
  //Zoltan2::AlltoAllCount(*problemComm,*env,sendCount, recvCount);
  
  //create sdispls and sentcount to build sendbuf
  //COMMENT THIS BACK IN
  std::vector<int> sdispls;
  sdispls.push_back(0);
  for(int i = 1; i < nprocs; i++){
    sdispls.push_back(sdispls[i-1] + sendcnts[i-1]);
  }

  std::vector<int> sentcount;
  for(int i = 0; i < nprocs; i++){
    sentcount.push_back(0);
  }
  int sendsize = 0;
  int recvsize = 0;
  for(int i = 0; i < nprocs; i++){
    sendsize += sendCount[i];
    recvsize += recvCount[i];
  }
  
  std::vector<gno_t> sendbuf(sendsize,0);
  for(int i = 0; i > num_boundary_edges*2; i+=2){
    if(owners[i] != me){
      int proc_to_send = owners[i];
      int sendbufidx = sdispls[proc_to_send] + sentcount[proc_to_send];
      sentcount[proc_to_send] += 2;
      sendbuf[sendbufidx++] = boundary[i];
      sendbuf[sendbufidx++] = boundary[i+1];
    }
    if(owners[i+1] != me && owners[i] != owners[i+1]){
      int proc_to_send = owners[i];
      int sendbufidx = sdispls[proc_to_send] + sentcount[proc_to_send];
      sentcount[proc_to_send] += 2;
      sendbuf[sendbufidx++] = boundary[i];
      sendbuf[sendbufidx++] = boundary[i+1];
    }
  }
  
  Teuchos::ArrayRCP<gno_t> recvbuf;
  Teuchos::ArrayView<const gno_t> sendBuf = Teuchos::arrayViewFromVector(sendbuf);
  Zoltan2::AlltoAllv(*problemComm,*env,sendBuf,sendCount,recvbuf,recvCount);

  //create the set to see if we've counted this boundary edge before
  struct pair_hash {
    inline std::size_t operator()(const std::pair<gno_t,gno_t>& v) const {
      return v.first * 10 + v.second;
    }
  };

  std::unordered_set<std::pair<gno_t,gno_t>,pair_hash> edge_set;
  int* local_boundary_counts = new int[nVtx];
  for(size_t i = 0; i < nVtx; i++){
    local_boundary_counts[i] = 0;
  }
  
  //insert the local boundary edges into the set, small global endpoint first.
  for(int i = 0; i < num_boundary_edges*2; i+=2){
    if(owners[i] == me || owners[i+1] == me){
      if(boundary_edges[i] < boundary_edges[i+1]){
        edge_set.insert(std::make_pair(boundary_edges[i],boundary_edges[i+1]));
      } else {
        edge_set.insert(std::make_pair(boundary_edges[i+1],boundary_edges[i]));
      }
    }
  }
  //factor in the locally-owned boundary edges
  for(int i = 0; i < num_boundary_edges*2; i++){
    if(map->getLocalElement(boundary_edges[i]) != fail){
      local_boundary_counts[map->getLocalElement(boundary_edges[i])]++;
    }
  }
  
  //use the received boundary edges, but check to make sure they haven't been used yet.
  for(int i = 0; i < recvsize; i+=2){
    //make sure the edge has not been used (if count returns 1 for one of these, the edge has been used)
    if(edge_set.count(std::make_pair(recvbuf[i], recvbuf[i+1])) == 0 && edge_set.count(std::make_pair(recvbuf[i+1],recvbuf[i])) == 0){
      //check which endpoints are owned locally, and set the boundary counts appropriately
      if(map->getLocalElement(recvbuf[i]) != fail){
        local_boundary_counts[map->getLocalElement(recvbuf[i])]++;
      }
      if(map->getLocalElement(recvbuf[i+1]) != fail){
        local_boundary_counts[map->getLocalElement(recvbuf[i+1])]++;
      }
      //add the received edge to the set, to ensure no received duplicates are counted
      if(recvbuf[i] < recvbuf[i+1]){
        edge_set.insert(std::make_pair(recvbuf[i],recvbuf[i+1]));
      } else {
	edge_set.insert(std::make_pair(recvbuf[i+1],recvbuf[i]));
      }
    }
  }
  std::cout<<me<<": done constructing boundary edges\n";
  //we should now have complete local knowledge of the boundary.

  //convert adjacency array to use local identifiers instead of global.
  
  typename map_t::local_ordinal_type* out_edges_lid = new typename map_t::local_ordinal_type[nEdge];
  for(size_t i = 0; i < nEdge; i++){
    out_edges_lid[i] = mapWithCopies->getLocalElement(out_edges[i]);
  }
  std::cout<<me<<": done creating out edges, creating csr graph\n";
  iceProp::graph<typename map_t::local_ordinal_type>* g = new iceProp::graph<typename map_t::local_ordinal_type>({nVtx, nEdge, out_edges_lid,out_offsets, 0,0.0});
  std::cout<<me<<": constructing propagation object\n";
  Zoltan2::iceSheetPropagation<map_t> prop(problemComm, map, mapWithCopies, g, local_boundary_counts, grounding, nVtx, nGhosts);
  std::cout<<me<<": starting propagation\n";  
  int* removed = prop.propagate();
  
  std::cout<<me<<": done propagating\n";
  //delete [] local_boundary_counts;
  //delete g;
  //delete [] grounding;

  return removed;
}

}//end namespace Zoltan2

#endif// _ZOLTAN2_ICESHEET_HPP
