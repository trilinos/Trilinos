//@HEADER
#ifndef _ZOLTAN2_ICESHEET_HPP_
#define _ZOLTAN2_ICESHEET_HPP_

#include<Zoltan2_GraphModel.hpp>
#include<Zoltan2_Algorithm.hpp>
#include<Zoltan2_PartitioningSolution.hpp>
#include<Zoltan2_Util.hpp>
#include<Zoltan2_TPLTraits.hpp>
#include<Zoltan2_VtxLabel.hpp>

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
      typedef Tpetra::Map<> map_t;
      typedef map_t::global_ordinal_type mgno_t;
      
      //Arguments: problemComm is the communicator we will use for the propagation
      //           adapter is the graph adapter that represents the ice mesh's bottom layer
      //           basalFriction is the array of basalFriction values. Nonzero indicates grounded
      //           boundary_edges is the array of edges on the current processor that are on the
      //                          boundary, using global identifiers
      //           num_boundary_edges is the number of boundary edges there are
      IceProp(const RCP<const Comm<int> > &problemComm__,
	      const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__, 
	      int* basalFriction, int* boundary_edges__, int num_boundary_edges__):
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
      void buildModel(modelFlag_t &flags);

      int* grounding_flags;
      int*  boundary_edges;
      int   num_boundary_edges;
      const RCP<const Comm<int> > problemComm;
      const RCP<const Environment> env;
      const RCP<const base_adapter_t> adapter;
      RCP<const GraphModel<base_adapter_t> > model;
  };
  
template <typename Adapter>
void IceProp<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
  flags.set(GENERATE_CONSECUTIVE_IDS);
  
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
  const size_t modelVerts = model->getLocalNumVertices();
  const size_t modelEdges = model->getLocalNumEdges();
  int num_verts = (int)modelVerts;
  long num_edges = (long)modelEdges;

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

  int* out_edges = NULL;
  unsigned* out_offsets = NULL;
  unsigned long* global_ids = NULL;
  TPL_Traits<int, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
  TPL_Traits<unsigned, const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);
  TPL_Traits<unsigned long, const gno_t>::ASSIGN_ARRAY(&global_ids, vtxIDs);
  
  int me = problemComm->getRank();

  /*std::cout<<me<<": Local vertices:\n";
  for(int i = 0; i < nVtx; i++){
    std::cout<<"\t"<<vtxIDs[i]<<"\n";
  }
  std::cout<<me<<": Local edges:\n";
  for(int i = 0; i < nVtx; i++){
    int outDegree = offsets[i+1] - offsets[i];
    for(int j = offsets[i]; j < offsets[i+1]; j++){
      std::cout<<me<<": "<<vtxIDs[i]<<" -- "<<adjs[j]<<"\n";
    }
  }*/
  
  //create the Tpetra map for the global indices
  Tpetra::global_size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  RCP<const Tpetra::Map<> > map = rcp(new Tpetra::Map<>(dummy, vtxIDs, 0, problemComm)); 
  
  //sentinel value to detect unsuccessful searches
  lno_t fail = Teuchos::OrdinalTraits<lno_t>::invalid();

  int nGhosts = 0;
  //count how many ghosts are in the edge list, to create the mapWithCopies.
  for(int i = 0; i < nEdge; i++){
    if(map->getLocalElement(adjs[i]) == fail)
      nGhosts++;
  }
  
  //std::cout<<me<<": number of ghosts = "<<nGhosts<<"\n";
  
  //use the count of ghosts + owned to make mapWithCopies, need to create a Teuchos array first.
  Teuchos::Array<gno_t> gids(nVtx+nGhosts);
  for(int i = 0; i < nVtx; i++){
    gids[i] = vtxIDs[i];
  }
  int ghostCount = 0;
  for(int i = 0; i < nEdge; i++){
    if(map->getLocalElement(adjs[i]) == fail){
      gids[nVtx+ghostCount] = adjs[i];
      ghostCount++;
    }
  }
  RCP<const Tpetra::Map<> > mapWithCopies = rcp(new Tpetra::Map<>(dummy, gids, 0, problemComm));
  
  int* local_boundary_counts = new int[nVtx];
  for(int i = 0; i < nVtx; i++){
    local_boundary_counts[i] = 0;
  }
  for(int i = 0; i < num_boundary_edges*2; i++){
    if(map->getLocalElement(boundary_edges[i]) != fail){
      local_boundary_counts[map->getLocalElement(boundary_edges[i])]++;
    }
  }
  
  //convert adjacency array to use local identifiers instead of global.
  for(int i = 0; i < nEdge; i++){
    out_edges[i] = mapWithCopies->getLocalElement(out_edges[i]);
  }

  graph* g = new graph({nVtx, nEdge, out_edges,out_offsets, 0,0.0});

  iceProp::iceSheetPropagation prop(problemComm, map, mapWithCopies, g, local_boundary_counts, grounding_flags, nVtx, nGhosts);
  
  int* removed = prop.propagate();
  
  delete local_boundary_counts;
  delete g;

  return removed;
}

}//end namespace Zoltan2

#endif// _ZOLTAN2_ICESHEET_HPP
