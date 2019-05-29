#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin, Manne, and Boman

namespace Zoltan2{

template <typename Adapter>

class AlgHybridGMB : public Algorithm<Adapter>
{
  private:
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;
    typedef typename Adapter::base_adapter_t base_adapter_t;
    typedef Tpetra::Map<lno_t, gno_t> map_t;

    void buildModel(modelFlag_t &flags);

    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;

  public:
    AlgHybridGMB(
      const RCP<const base_adapter_t> &adapter_, 
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_) {
      
      modelFlag_t flags;
      flags.reset();
      buildModel(flags);
      
    }

    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution ) {
      
      //this will color the global graph in a manner similar to Zoltan
      
      //get vertex GIDs in a locally indexed array (stolen from Ice-Sheet interface)
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      //we do not use weights at this point

      //get edge information from the model
      ArrayView<const gno_t> adjs;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);
      //again, weights are not used

      //create maps for FEMultiVector
      gno_t* out_edges = NULL;
      offset_t* out_offsets = NULL;
      gno_t* global_ids = NULL;
      TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
      TPL_Traits<offset_t, const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);
      TPL_Traits<gno_t, const gno_t>::ASSIGN_ARRAY(&global_ids, vtxIDs);

      //create Tpetra map for the global indices
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      RCP<const map_t> map = rcp(new map_t(dummy, vtxIDs, 0, comm));
      
      //sentinel value to detect unsuccessful searches
      lno_t fail = Teuchos::OrdinalTraits<lno_t>::invalid();

      int nGhosts = 0;
      //count how many ghosts are in the edge list to create the mapWithCopies
      for(size_t i = 0; i< nEdge; i++){
        if(map->getLocalElement(adjs[i]) == fail)
          nGhosts++;
      }

      //use the count of ghosts + owned to make mapWithCopies
      Teuchos::Array<gno_t> gids(nVtx+nGhosts);
      for(size_t i = 0; i < nVtx; i++){
        gids[i] = vtxIDs[i];
      }
      int ghostCount = 0;
      for(size_t i = 0; i < nEdge; i++){
        if(map->getLocalElement(adjs[i]) == fail){
          gids[nVtx+ghostCount] = adjs[i];
          ghostCount++;
        }
      }
      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy, gids, 0, comm));
      
      /*printf("--Rank %d has global(local) verts: ",comm->getRank());
      for(size_t i = 0; i < nVtx+nGhosts;i++){
        printf("%u(%u) ", gids[i], mapWithCopies->getLocalElement(gids[i]));
      }
      printf("\n");*/
      
      //convert edges to local ID, as I'm pretty sure Kokkos-Kernel's coloring expects this.
      printf("--Rank %d local adjacencies\n",comm->getRank());
      std::vector<lno_t> local_adjs_vec;
      for(size_t i = 0; i < nEdge; i++){
        //printf("\t%u\n",mapWithCopies->getLocalElement(out_edges[i]));
        local_adjs_vec.push_back( mapWithCopies->getLocalElement(out_edges[i]));
      }
      //for(size_t i = 0; i < nEdge; i++) printf("\t%u\n",local_adjs_vec.at(i));
      ArrayView<const lno_t> local_adjs = Teuchos::arrayViewFromVector(local_adjs_vec);
      //TODO: create FEMultiVector of some type, need to figure out what is appropriate.
      //relevant lines from VtxLabel:
      //   typedef Tpetra::Import<lno_t, gno_t> import_t;
      //   typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;
      //                                 ^-----Need to figure out what type this should be
      //   Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, mapWithCopies))
      //   femv_t femv = rcp(new femv_t(mapOwned, importer, 1, true));
      //                                                    ^---could change, potentially? (#vectors in multivector)


      //Get color array to fill
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
        colors[i] = 0;
      } 
            
      
      // call actual coloring function
      // THESE ARGUMENTS WILL NEED TO CHANGE,
      // THESE ARE A COPY OF THE EXISTING ALGORITHM CLASS.
      hybridGMB(nVtx, local_adjs, offsets,colors);
    }
    
    void hybridGMB(const size_t nVtx, ArrayView<const lno_t> adjs, 
                   ArrayView<const offset_t> offsets, ArrayRCP<int> colors){
      //print statements for debugging purposes:
       
      //color the interior vertices (maybe)

      //color boundary vertices using FEMultiVector (definitely)

      //color interior vertices if not colored yet. (must alter Kokkos-Kernels to allow partial colorings)
      //As long as the initial coloring is untouched, we can pass the whole graph to the kokkos-kernels coloring.
      
    }
};

template <typename Adapter>
void AlgHybridGMB<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
  
  this->env->debug(DETAILED_STATUS, "   building graph model");
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
  this->env->debug(DETAILED_STATUS, "   graph model built");
}

}//end namespace Zoltan2
#endif
