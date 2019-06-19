#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include <KokkosKernels_Handle.hpp>
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"

#include "zz_rand.h"
//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin, Manne, 
//!        and Boman

namespace Zoltan2{

//new class, colorLabel?
class VtxColor {
  public:
    unsigned long long rand;
    unsigned long long color;
    std::vector<VtxColor*> nbors;
    //constructors
    VtxColor(unsigned long long r, unsigned long long c):rand(r),color(c){}
    VtxColor(){
      rand = 0;
      color = 0;
    }
    //constructors for FEMultiVector's health
    VtxColor(const unsigned long long l){
      rand = l;
      color = 0;
    }
    VtxColor(volatile const VtxColor& other){
      rand = other.rand;
      color = other.color;
      //nbors = other.nbors;
    }
    VtxColor(const VtxColor& other){
      rand = other.rand;
      color = other.color;
      //nbors = other.nbors;
    }
    volatile VtxColor operator=(const VtxColor& other) volatile {
      rand = other.rand;
      color = other.color;
      //nbors = other.nbors;
      return *this;
    }
    VtxColor& operator=(const VtxColor& other){
      rand = other.rand;
      color = other.color;
      //nbors = other.nbors;
      return *this;
    }
    VtxColor& operator+=(const VtxColor& copy){
      //conflict detection is here.
      color = copy.color;
      return *this;
    }
    friend bool operator==(const VtxColor& lhs, const VtxColor& rhs){
      return (lhs.rand == rhs.rand) && (lhs.color == rhs.color) &&
             (lhs.nbors == rhs.nbors);
    }
    friend std::ostream& operator<<(std::ostream& os, const VtxColor& a){
      os<<"color: "<<a.color<<", rand: "<<a.rand;
      return os;
    }
   
};

}//end namespace Zoltan2

namespace Kokkos{
  namespace Details {
    template<>
    class ArithTraits<Zoltan2::VtxColor>{
    public:
      typedef Zoltan2::VtxColor val_type;
      typedef unsigned long long mag_type;
      
      static const bool is_specialized = true;
      static const bool is_signed = false;
      static const bool is_integer = true;
      static const bool is_exact = true;
      static const bool is_complex = false;

      static KOKKOS_FORCEINLINE_FUNCTION bool isInf(const val_type&){
        return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isNan(const val_type &){
        return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x){
        return x.rand;
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type zero() {return 0;}
      static KOKKOS_FORCEINLINE_FUNCTION mag_type one() {return 1;}
      static KOKKOS_FORCEINLINE_FUNCTION mag_type min() {return LONG_MIN;}
      static KOKKOS_FORCEINLINE_FUNCTION mag_type max() {return LONG_MAX;}
      static KOKKOS_FORCEINLINE_FUNCTION mag_type nan() {
        return (unsigned long long)-1;
      }
      
      //Backwards compatibility with Teuchos::ScalarTraits
      typedef mag_type magnitudeType;
      static const bool isComplex = false;
      static const bool isOrdinal = true;
      static const bool isComparable = true;
      static const bool hasMachineParameters = false;
      static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude(
        const val_type &x){
        return abs(x);
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf(const val_type&){
        return false;
      }
      static std::string name() {return "Zoltan2::VtxColor";}
    };
  }//end namespace Details
}//end namespace Kokkos

namespace Teuchos {
template<typename Ordinal>
struct SerializationTraits<Ordinal, Zoltan2::VtxColor> :
       public Teuchos::DirectSerializationTraits<Ordinal, Zoltan2::VtxColor>
{};
}//end namespace Teuchos

namespace Zoltan2 {

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
    typedef int femv_scalar_t;
    typedef Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t> femv_t;
    
  private:

    void buildModel(modelFlag_t &flags);
    std::vector<gno_t> reorderGraph(
                                 ArrayView<const gno_t> ownedGIDs,
                                 ArrayView<const gno_t> globalAdjs, 
                                 ArrayView<const offset_t> offsets, 
                                 std::vector<lno_t>& reorderAdjs, 
                                 std::vector<offset_t>&reorderOffsets,
                                 std::vector<lno_t>& reorderToLocal,
                                 lno_t& nInterior) {
      std::unordered_map<gno_t, lno_t> globalToLocal; 
      
      printf("--Rank %d: setting up map\n", comm->getRank()); 
      //map globally owned ids to local ids. O(n)
      for( size_t i = 0; i < ownedGIDs.size(); i++){
        globalToLocal[ownedGIDs[i]] = i;
      }
      
      //map ghosts to local ids. O(m)
      lno_t ghostCount = 0;
      std::vector<gno_t> finalLocalToGlobal(ownedGIDs.size());
      for(size_t i = 0; i < globalAdjs.size(); i++){
        if(globalToLocal.count(globalAdjs[i])==0){
          //this is a ghost GID, put it in the map, assign it a local ID
          globalToLocal[globalAdjs[i]] = ownedGIDs.size() + ghostCount;
          ghostCount++;
          //also put it in the correct spot in the final vector.
          finalLocalToGlobal.push_back(globalAdjs[i]);
          printf("--Rank %d: put %u at index %u in reorderGIDs\n",
                 comm->getRank(), globalAdjs[i], finalLocalToGlobal.size());
        }
      }
      
      printf("--Rank %d: deterimining interior and boundary\n",comm->getRank());
      //determine interior vs. boundary vertices. O(m)
      bool* interiorFlags = new bool[ownedGIDs.size()];
      nInterior = 0;
      for(size_t i = 0; i < ownedGIDs.size(); i++){
        bool interior = true;
        for(size_t j = offsets[i]; j < offsets[i+1]; j++){
          if(globalToLocal[globalAdjs[j]] >= ownedGIDs.size()){
            //not interior
            interior = false;
          }
        }
        //can rewrite this for readability.
        nInterior += (interiorFlags[i] = interior);
      }

      printf("--Rank %d: changing to reordered LIDs\n", comm->getRank());
      //change the map to use reordered LIDs, and set up other necessary
      //structures for reorderOffsets. O(n)
      lno_t interiorCount = 0;
      lno_t boundaryCount = nInterior;
      offset_t* reorderDegrees = new offset_t[ownedGIDs.size()];
      for(size_t i = 0; i < ownedGIDs.size(); i++){
        if(interiorFlags[i]){
          //interiorCount is the reordered LID
          finalLocalToGlobal[interiorCount] = ownedGIDs[i];
          globalToLocal[ownedGIDs[i]] = interiorCount;
          reorderDegrees[interiorCount++] = offsets[i+1] - offsets[i]; 
          printf("--Rank %d: put %u at index %d of reorderGIDs\n",
                 comm->getRank(), ownedGIDs[i], interiorCount-1);
        } else {
          //boundaryCount is the reordered LID
          finalLocalToGlobal[boundaryCount] = ownedGIDs[i];
          globalToLocal[ownedGIDs[i]] = boundaryCount;
          reorderDegrees[boundaryCount++] = offsets[i+1] - offsets[i];
          printf("--Rank %d: put %u at index %d of reorderGIDs\n",
                 comm->getRank(), ownedGIDs[i], boundaryCount-1);
        }
      }
      
      printf("--Rank %d: computing reordered offsets\n",comm->getRank());
      //compute reorderOffsets. O(n)
      reorderOffsets[0] = 0;
      for(size_t i = 0; i < ownedGIDs.size(); i++){
        reorderOffsets[i+1] = reorderOffsets[i] + reorderDegrees[i];
      }
 
      printf("--Rank %d: computing reordered adjacencies\n",comm->getRank());
      if(reorderAdjs.size() != globalAdjs.size()){
        printf("--Rank %d: size mismatch between original and reordered adjs\n"
              , comm->getRank());
      }
      //compute reorderAdjs. O(m)
      for(size_t i = 0; i < ownedGIDs.size(); i++){
        for(size_t j = 0; j < offsets[i+1]-offsets[i]; j++){
          reorderAdjs.at(reorderOffsets[globalToLocal[ownedGIDs[i]]]+j) = 
            globalToLocal[globalAdjs[offsets[i]+j]];
        }
      }
      
      reorderToLocal.resize(ownedGIDs.size());
      for(size_t i = 0; i < ownedGIDs.size(); i++){
        //printf("--Rank %d: reordered ID %d is local ID %d\n",comm->getRank(),globalToLocal[ownedGIDs[i]],i);
        reorderToLocal[globalToLocal[ownedGIDs[i]]] = i;
      }
      
      printf("--Rank %d: done reordering graph\n", comm->getRank());
      return finalLocalToGlobal;
    } 
    
    
    
    //This function is templated twice to avoid explicitly coding all possible
    //Execution space situations. We can just template the function call and 
    //have the compiler do that work.
    template <class ExecutionSpace, typename TempMemorySpace, 
              typename MemorySpace>
    void colorInterior(const size_t nVtx, 
                       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > adjs_view,
                       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > offset_view, 
                       Teuchos::ArrayRCP<int> colors,
                       Teuchos::RCP<femv_t> femv){
      //typedef Kokkos::Device<ExecutionSpace, MemorySpace> device;
      //Kokkos::View<offset_t*, device> offset_view("degree_offsets", 
      //                                            offsets.size());
      //Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
      //copy from the arrayviews into the Kokkos::Views
      //Kokkos::parallel_for("Init_Offsets",offsets.size(),
      //                     KOKKOS_LAMBDA (const int& i) {
      //  offset_view(i) = offsets[i];
      //});
      //Kokkos::parallel_for("Init_Adjacencies", adjs.size(), 
      //                     KOKKOS_LAMBDA (const int& i) {
      //  adjs_view(i) = adjs[i];
      //});
      
      //default values are taken from KokkosKernels_TestParameters.hpp
      
      int algorithm = pl->get<int>("Hybrid_algorithm",0);
      int repeat = pl->get<int>("Hybrid_repeat",1); //probably not using this 
                                                    //until we test performance
      int chunk_size = pl->get<int>("Hybrid_chunk_size",-1); 
      int shmemsize = pl->get<int>("Hybrid_shmemsize", 16128);
      int team_size = pl->get<int>("Hybrid_team_size", -1);
      int use_dynamic_scheduling = pl->get<int>("Hybrid_use_dynamic_scheduling",
                                                0);
      int verbose = pl->get<int>("Hybrid_verbose",0);

      int vector_size = pl->get<int>("Hybrid_vector_size",-1);

      typedef KokkosKernels::Experimental::KokkosKernelsHandle
          <size_t, lno_t, lno_t, ExecutionSpace, TempMemorySpace, 
           MemorySpace> KernelHandle;
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
     
      //set the initial coloring of the kh.get_graph_coloring_handle() to be
      //the data view from the femv.
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<MemorySpace>();
      //printf("--Rank %d: femv view reference count before subview: %d\n",comm->getRank(),femvColors.use_count());
      Kokkos::View<int*, Tpetra::Map<>::device_type >  sv = subview(femvColors, Kokkos::ALL, 0);
      //printf("--Rank %d: femv view reference count after subview: %d\n",comm->getRank(),femvColors.use_count());
      //Kokkos::View<int*,Tpetra::Map<>::device_type > device_color("device_Color",sv.extent(0));
      //Kokkos::deep_copy(device_color, sv); 
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      //printf("--Rank %d: femv view reference count after set vertex colors: %d\n",comm->getRank(),femvColors.use_count());
      
      KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx, nVtx, 
                                                      offset_view, adjs_view);
      Kokkos::fence();
      std::cout << std::endl << 
          "Time: " << 
          kh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
          "Num colors: " << 
          kh.get_graph_coloring_handle()->get_num_colors() << " "
          "Num Phases: " << 
          kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
      std::cout << "\t"; 
      KokkosKernels::Impl::print_1Dview(
                           kh.get_graph_coloring_handle()->get_vertex_colors());
      
      numColors = kh.get_graph_coloring_handle()->get_num_colors();
      //Kokkos::deep_copy(sv, device_color); 
      //all this code should be unnecessary if we can pass in the Kokkos::View to the
      //coloring correctly.
      //auto host_view = Kokkos::create_mirror_view(
      //                   kh.get_graph_coloring_handle()->get_vertex_colors());
      //Kokkos::deep_copy(host_view, 
      //                   kh.get_graph_coloring_handle()->get_vertex_colors());
      //auto nr = host_view.extent(0);
      //for(auto i = 0; i < nr; i++){ 
        //colors[i] = host_view(i);
        //femv->replaceLocalValue(i,0, host_view(i)); 
      //}
    }
    
    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    int numColors;
    
  public:
    AlgHybridGMB(
      const RCP<const base_adapter_t> &adapter_, 
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_) {

      numColors = 4;
      modelFlag_t flags;
      flags.reset();
      buildModel(flags);

    }


    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution ) {
      
      //this will color the global graph in a manner similar to Zoltan
      
      //get vertex GIDs in a locally indexed array (stolen from Ice-Sheet 
      //interface)
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

      //reorder the graph so that boundary vertices are in the
      //end of the offset array (saves time later) 
      lno_t nInterior;
      std::vector<lno_t> reorderAdjs_vec(adjs.size());
      std::vector<offset_t> reorderOffsets_vec(offsets.size());
      std::vector<lno_t> reorderToLocal;
      std::vector<gno_t> reorderGIDs = reorderGraph(vtxIDs,
                                                      adjs,
                                                      offsets,
                                                      reorderAdjs_vec,
                                                      reorderOffsets_vec,
                                                      reorderToLocal,
                                                      nInterior); 
      ArrayView<const lno_t> reorderAdjs = Teuchos::arrayViewFromVector(
                                                         reorderAdjs_vec);
      ArrayView<const offset_t> reorderOffsets = Teuchos::arrayViewFromVector(
                                                         reorderOffsets_vec);
     
      std::vector<gno_t> ownedReorderGIDs;
      for(int i = 0; i < nVtx; i++){
        ownedReorderGIDs.push_back(reorderGIDs[i]);
      }
       
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                             <Tpetra::global_size_t>::invalid();
      RCP<const map_t> mapOwned = rcp(new map_t(dummy, 
                                      Teuchos::arrayViewFromVector(ownedReorderGIDs),
                                           0, comm));
      RCP<const map_t> mapWithCopies = rcp(new map_t(dummy, 
                                      Teuchos::arrayViewFromVector(reorderGIDs),
                                           0, comm));
      
       
      //TODO: create FEMultiVector of some type, need to figure out what is
      //appropriate.
      //relevant lines from VtxLabel:
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      //   typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;
      //                                 ^-----Need to figure out what type 
      //                                       this should be
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, 
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, 
                                                    importer, 1, true));
      //                                 could change,--------^ 
      //                                 potentially? (#vectors in multivector)

      //Get color array to fill
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
        colors[i] = 0;
      } 
      
      //taken directly from the Zoltan coloring implementation 
      std::vector<int> rand(reorderGIDs.size());
      for(int i = 0; i < reorderGIDs.size(); i++){
        Zoltan_Srand((unsigned int) reorderGIDs[i], NULL);
        rand[i] = (int) (((double) Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
      }

      // call actual coloring function
      // THESE ARGUMENTS WILL NEED TO CHANGE,
      // THESE ARE A COPY OF THE EXISTING ALGORITHM CLASS.
      hybridGMB(nVtx, nInterior, reorderAdjs, reorderOffsets,colors,femv,reorderGIDs,rand);
      
      for(int i = 0; i < colors.size(); i++){
        colors[reorderToLocal[i]] = femv->getData(0)[i];
      }
     
      //test validation on simple
      /*if(comm->getRank() == 0){
        colors[reorderToLocal[12]] = 3;
      }
      if(comm->getRank() == 1){
        colors[reorderToLocal[11]] = 3;
      }*/
       
      comm->barrier();
    }
    
    void hybridGMB(const size_t nVtx,lno_t nInterior, Teuchos::ArrayView<const lno_t> adjs, 
                   Teuchos::ArrayView<const offset_t> offsets, 
                   Teuchos::ArrayRCP<int> colors, Teuchos::RCP<femv_t> femv,
                   std::vector<gno_t> reorderGIDs,
                   std::vector<int> rand){
      
      //make host views to deep_copy into the device views
      
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> host_offsets("Host Offset view", offsets.size());
      for(int i = 0; i < offsets.size(); i++){
        host_offsets(i) = offsets[i];
      }
      Kokkos::View<lno_t*, Tpetra::Map<>::device_type> host_adjs("Host Adjacencies view", adjs.size());
      for(int i = 0; i < adjs.size(); i++){
        host_adjs(i) = adjs[i];
      }
      
      this->colorInterior<Tpetra::Map<>::execution_space,
                          Tpetra::Map<>::memory_space,
                          Tpetra::Map<>::memory_space>
                 (nInterior, host_adjs, host_offsets, colors, femv);
      
      /*bool use_cuda = pl->get<bool>("Hybrid_use_cuda",false);
      bool use_openmp = pl->get<bool>("Hybrid_use_openmp",false);
      bool use_serial = pl->get<bool>("Hybrid_use_serial",false);
      //do the kokkos stuff in here (assume that the user checked for the option they selected)
      //Might want to enable multi-memory stuff?
      if(use_cuda + use_openmp + use_serial == 0){
        typedef Kokkos::Device<Kokkos::DefaultExecutionSpace, 
                               Kokkos::DefaultExecutionSpace::memory_space> device;
        Kokkos::View<offset_t*, device> offset_view("degree_offsets", 
                                                    offsets.size());
        Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
        //copy from the arrayviews into the Kokkos::Views
        Kokkos::deep_copy(offset_view,host_offsets);
        Kokkos::deep_copy(adjs_view, host_adjs);

        //use the default spaces to run the KokkosKernels coloring
        this->colorInterior<Kokkos::DefaultExecutionSpace,
                            Kokkos::DefaultExecutionSpace::memory_space,
                            Kokkos::DefaultExecutionSpace::memory_space> 
                     (nInterior, adjs_view, offset_view, colors,femv);
        //try to use femv's execution space stuff?
        this->colorInterior<femv_t::execution_space,
                            femv_t::device_type::memory_space,
                            femv_t::device_type::memory_space>
                   (nInterior, adjs_view, offset_view, colors, femv);
      } else if (use_cuda) {
        //use the cuda spaces
        #ifdef KOKKOS_ENABLE_CUDA
        typedef Kokkos::Device<Kokkos::Cuda, 
                               Kokkos::Cuda::memory_space> device;
        Kokkos::View<offset_t*, device> offset_view("degree_offsets", 
                                                    offsets.size());
        Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
        //copy from the arrayviews into the Kokkos::Views
        Kokkos::deep_copy(offset_view, host_offsets);
        Kokkos::deep_copy(adjs_view, host_adjs);

        this->colorInterior<Kokkos::Cuda, Kokkos::Cuda::memory_space, 
                            Kokkos::Cuda::memory_space>
                     (nInterior, adjs_view, offset_view, colors,femv);
        #endif
      } else if (use_openmp) {
        //use openmp spaces
        #ifdef KOKKOS_ENABLE_OPENMP
        typedef Kokkos::Device<Kokkos::OpenMP, 
                               Kokkos::OpenMP::memory_space> device;
        Kokkos::View<offset_t*, device> offset_view("degree_offsets", 
                                                    offsets.size());
        Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
        //copy from the arrayviews into the Kokkos::Views
        Kokkos::deep_copy(offset_view, host_offsets);
        Kokkos::deep_copy(adjs_view, host_adjs);

        this->colorInterior<Kokkos::OpenMP, Kokkos::OpenMP::memory_space, 
                            Kokkos::OpenMP::memory_space>
                     (nInterior, adjs_view, offset_view, colors,femv);
        #endif
      } else if (use_serial) {
        //use serial spaces

        typedef Kokkos::Device<Kokkos::Serial, 
                               Kokkos::Serial::memory_space> device;
        Kokkos::View<offset_t*, device> offset_view("degree_offsets", 
                                                    offsets.size());
        Kokkos::View<lno_t*, device> adjs_view("adjacency_list",adjs.size());
        //copy from the arrayviews into the Kokkos::Views
        Kokkos::deep_copy(offset_view, host_offsets);
        Kokkos::deep_copy(adjs_view, host_adjs);

        this->colorInterior<Kokkos::Serial, Kokkos::Serial::memory_space, 
                            Kokkos::Serial::memory_space> 
                     (nInterior, adjs_view, offset_view, colors,femv);
      }*/
       
      int batch_size = pl->get<int>("Hybrid_batch_size",100);
      //color boundary vertices using FEMultiVector (definitely)
      //create a queue of vertices to color
      std::queue<lno_t> recoloringQueue;
      std::queue<lno_t> conflictQueue;
      Teuchos::ArrayRCP<bool> forbiddenColors;
      forbiddenColors.resize(numColors);
      for(int i = 0; i < numColors; i++) forbiddenColors[i] = false;
      printf("--Rank %d: batch size: %d\n",comm->getRank(),batch_size);
      //for each batch of boundary vertices
      //femv->switchActiveMultiVector();
      bool done = false;
      int i = nInterior;
      //  while the queue is not empty
          while(recoloringQueue.size() > 0 || !done){
        //  add next batch (an argument to this class?) to the queue
            for(size_t j = i; j < i+batch_size; j++){
              if(j < nVtx) {
                printf("--Rank %d: pushing %d on the queue\n",comm->getRank(),j);
                recoloringQueue.push(j);
              }
            }
            if(i < nVtx) i+= batch_size;
      //    color everything in the recoloring queue, put everything on conflict queue
            while(recoloringQueue.size() > 0) {
               printf("--Rank %d: coloring vertex %u\n",comm->getRank(),recoloringQueue.front());
               lno_t currVtx = recoloringQueue.front();
               recoloringQueue.pop();
               conflictQueue.push(currVtx);
               for(offset_t nborIdx = offsets[currVtx]; nborIdx < offsets[currVtx+1]; nborIdx++){
                 
                 int nborColor = femv->getData(0)[adjs[nborIdx]];
                 if(nborColor > 0){
                   forbiddenColors[nborColor-1] = true;
                 } 
               }
               //pick the color for currVtx based on the forbidden array
               bool colored = false;
               for(int i = 0; i < numColors; i++){
                 if(!forbiddenColors[i]) {
                   femv->replaceLocalValue(currVtx,0,i+1);
                   //colors[currVtx] = i+1;
                   colored = true;
                   break;
                 }
               }
               for(int i = 0; i < numColors; i++){
                 forbiddenColors[i] = false;
               }
               if(!colored){
                 //colors[currVtx] = numColors+1;
                 femv->replaceLocalValue(currVtx,0,numColors+1);
                 numColors++;
                 forbiddenColors.resize(numColors);
               }
               printf("--Rank %d: colored vtx %u color %d\n",comm->getRank(),currVtx,colors[currVtx]);
            }
            
      //    communicate
            femv->switchActiveMultiVector();
            femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
            femv->switchActiveMultiVector();
      //    check for conflicts with the vertices that were just colored
            while(conflictQueue.size() > 0){
      //      if conflicts were detected put them on the recoloring queue
              lno_t currVtx = conflictQueue.front();
              conflictQueue.pop();
              bool conflict = false;
              int currColor = femv->getData(0)[currVtx];
              for(offset_t nborIdx = offsets[currVtx]; nborIdx < offsets[currVtx+1]; nborIdx++){
                int nborColor = femv->getData(0)[adjs[nborIdx]];
                if(nborColor == currColor && rand[currVtx] > rand[adjs[nborIdx]]) {
                  printf("--Rank %d: vtx %u conflicts with vtx %u\n",comm->getRank(),currVtx,adjs[nborIdx]);
                  conflict = true;
                }
              }
              if(conflict) {
                recoloringQueue.push(currVtx);
                printf("--Rank %d: putting vertex %u on the recoloring queue\n",comm->getRank(),currVtx);
              }
            }
            //do a reduction to determine if we're done
            int globalDone = 0;
            int localDone = recoloringQueue.size() + (nVtx > i);
            //comm->reduceAll(Teuchos::REDUCE_SUM,sizeof(int),&localDone, &globalDone);
            Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_MAX,1, &localDone, &globalDone);
            //We're only allowed to stop once everyone has no work to do.
            //collectives will hang if one process exits. 
            done = !globalDone;
          }
        
       
      //need to handle half-done batches?
      
      //color interior vertices if not colored yet. 
      //(must alter Kokkos-Kernels to allow partial colorings)
      //As long as the initial coloring is untouched, we can pass the whole 
      //graph to the kokkos-kernels coloring.
      
      //print the coloring
      for(int i = 0; i < femv->getData(0).size(); i++){
        printf("--Rank %d: local vtx %u is color %d\n",comm->getRank(),i,femv->getData(0)[i]);
      }
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
