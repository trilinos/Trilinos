#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include <KokkosKernels_Handle.hpp>
#include "KokkosKernels_IOUtils.hpp"
//#include "KokkosKernels_TestParameters.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"
//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin, Manne, and Boman

namespace Zoltan2{

template <typename Adapter>

class AlgHybridGMB : public Algorithm<Adapter>
{
  public:
  
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;
    typedef typename Adapter::base_adapter_t base_adapter_t;
    typedef Tpetra::Map<lno_t, gno_t> map_t;

  private:

    void buildModel(modelFlag_t &flags);
    
    //template <typename ExecutionSpace, typename TempMemorySpace, typename MemorySpace>
    //void colorInterior(const size_t, Teuchos::ArrayView<const lno_t>, 
    //                   Teuchos::ArrayView<const offset_t>, Teuchos::ArrayRCP<int>);
    //This function is templated twice to avoid explicitly coding all possible
    //Execution space situations. We can just template the function call and have the compiler do that work.
    template <typename ExecutionSpace, typename TempMemorySpace, typename MemorySpace>
    void colorInterior(const size_t nVtx, Teuchos::ArrayView<const lno_t> adjs,
                                       Teuchos::ArrayView<const offset_t> offsets, Teuchos::ArrayRCP<int> colors){ 
      typedef Kokkos::Device<ExecutionSpace, MemorySpace> device;
      Kokkos::View<offset_t*, device> offset_view("degree_offsets", offsets.size());
      Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
      //copy from the arrayviews into the Kokkos::Views
      Kokkos::parallel_for("Init_Offsets",offsets.size(), KOKKOS_LAMBDA (const int& i) {
        offset_view(i) = offsets[i];
      });
      Kokkos::parallel_for("Init_Adjacencies", adjs.size(), KOKKOS_LAMBDA (const int& i) {
        adjs_view(i) = adjs[i];
      });
      
      //default values are taken from KokkosKernels_TestParameters.hpp
      
      int algorithm = pl->get<int>("Hybrid_algorithm",0);
      int repeat = pl->get<int>("Hybrid_repeat",1); //probably not using this until we test performance
      int chunk_size = pl->get<int>("Hybrid_chunk_size",-1); 
      int shmemsize = pl->get<int>("Hybrid_shmemsize", 16128);
      int team_size = pl->get<int>("Hybrid_team_size", -1);
      int use_dynamic_scheduling = pl->get<int>("Hybrid_use_dynamic_scheduling",0);
      int verbose = pl->get<int>("Hybrid_verbose",0);

      int vector_size = pl->get<int>("Hybrid_vector_size",-1);

      typedef KokkosKernels::Experimental::KokkosKernelsHandle
          <size_t, lno_t, lno_t, ExecutionSpace, TempMemorySpace, MemorySpace> KernelHandle;
      KernelHandle kh;

      kh.set_team_work_size(chunk_size);
      kh.set_shmem_size(shmemsize);
      kh.set_suggested_team_size(team_size);
      kh.set_suggested_vector_size(vector_size);
  
      kh.set_dynamic_scheduling(use_dynamic_scheduling);
      kh.set_verbose(verbose);
    
      switch(algorithm) {
        case 1:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
          break;
        case 2:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_SERIAL);
          break;
        case 3:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_VB);
          break;
        case 4:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
          break;
        case 5:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBCS);
          break;
        case 6:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
          break;
        case 7:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBD);
          break;
        case 8:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBDBIT);
          break;
        default:
          kh.create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
      }
     
      KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx,nVtx,offset_view,adjs_view);
      std::cout << std::endl << 
          "Time: " << kh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
          "Num colors: " << kh.get_graph_coloring_handle()->get_num_colors() << " "
          "Num Phases: " << kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
      std::cout << "\t"; KokkosKernels::Impl::print_1Dview(kh.get_graph_coloring_handle()->get_vertex_colors());
    }
    
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
      
      comm->barrier();
    }
    
    void hybridGMB(const size_t nVtx, Teuchos::ArrayView<const lno_t> adjs, 
                   Teuchos::ArrayView<const offset_t> offsets, Teuchos::ArrayRCP<int> colors){
      
      //need to create a subgraph that has only interior vertices.
      //this graph should get passed to the Kokkos calls.
      //Need to have a map from the interior ID to the local ID, so that the colors can be
      //put in the colors array correctly.
      
       
      bool use_cuda = pl->get<bool>("Hybrid_use_cuda",false);
      bool use_openmp = pl->get<bool>("Hybrid_use_openmp",false);
      bool use_serial = pl->get<bool>("Hybrid_use_serial",false);
      //do the kokkos stuff in here (assume that the user checked for the option they selected)
      //Might want to enable multi-memory stuff?
      if(use_cuda + use_openmp + use_serial == 0){
        //use the default spaces to run the KokkosKernels coloring
        this->colorInterior<Kokkos::DefaultExecutionSpace,Kokkos::DefaultExecutionSpace::memory_space,
                      Kokkos::DefaultExecutionSpace::memory_space> (nVtx, adjs, offsets, colors);
      } else if (use_cuda) {
        //use the cuda spaces
        #ifdef KOKKOS_ENABLE_CUDA
        this->colorInterior<Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>
                     (nVtx, adjs, offsets, colors);
        #endif
      } else if (use_openmp) {
        //use openmp spaces
        #ifdef KOKKOS_ENABLE_OPENMP
        this->colorInterior<Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>
                     (nVtx, adjs, offsets, colors);
        #endif
      } else if (use_serial) {
        //use serial spaces
        this->colorInterior<Kokkos::Serial, Kokkos::Serial::memory_space, Kokkos::Serial::memory_space>
                     (nVtx, adjs, offsets, colors);
      }
       
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
