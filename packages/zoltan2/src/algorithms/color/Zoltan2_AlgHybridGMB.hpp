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
      
      //printf("--Rank %d: setting up map\n", comm->getRank()); 
      //map globally owned ids to local ids. O(n)
      for( uint64_t i = 0; i < ownedGIDs.size(); i++){
        globalToLocal[ownedGIDs[i]] = i;
      }
      
      //map ghosts to local ids. O(m)
      lno_t ghostCount = 0;
      std::vector<gno_t> finalLocalToGlobal(ownedGIDs.size());
      for(uint64_t i = 0; i < globalAdjs.size(); i++){
        if(globalToLocal.count(globalAdjs[i])==0){
          //this is a ghost GID, put it in the map, assign it a local ID
          globalToLocal[globalAdjs[i]] = ownedGIDs.size() + ghostCount;
          ghostCount++;
          //also put it in the correct spot in the final vector.
          finalLocalToGlobal.push_back(globalAdjs[i]);
          //printf("--Rank %d: put %u at index %u in reorderGIDs\n",
          //      comm->getRank(), globalAdjs[i], finalLocalToGlobal.size());
        }
      }
      
      //printf("--Rank %d: deterimining interior and boundary\n",comm->getRank());
      //determine interior vs. boundary vertices. O(m)
      bool* interiorFlags = new bool[ownedGIDs.size()];
      nInterior = 0;
      for(uint64_t i = 0; i < ownedGIDs.size(); i++){
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

      //printf("--Rank %d: changing to reordered LIDs\n", comm->getRank());
      //change the map to use reordered LIDs, and set up other necessary
      //structures for reorderOffsets. O(n)
      lno_t interiorCount = 0;
      lno_t boundaryCount = nInterior;
      offset_t* reorderDegrees = new offset_t[ownedGIDs.size()];
      for(uint64_t i = 0; i < ownedGIDs.size(); i++){
        if(interiorFlags[i]){
          //interiorCount is the reordered LID
          finalLocalToGlobal[interiorCount] = ownedGIDs[i];
          globalToLocal[ownedGIDs[i]] = interiorCount;
          reorderDegrees[interiorCount++] = offsets[i+1] - offsets[i]; 
          //printf("--Rank %d: put %u at index %d of reorderGIDs\n",
          //       comm->getRank(), ownedGIDs[i], interiorCount-1);
        } else {
          //boundaryCount is the reordered LID
          finalLocalToGlobal[boundaryCount] = ownedGIDs[i];
          globalToLocal[ownedGIDs[i]] = boundaryCount;
          reorderDegrees[boundaryCount++] = offsets[i+1] - offsets[i];
          //printf("--Rank %d: put %u at index %d of reorderGIDs\n",
          //       comm->getRank(), ownedGIDs[i], boundaryCount-1);
        }
      }
      
      //printf("--Rank %d: computing reordered offsets\n",comm->getRank());
      //compute reorderOffsets. O(n)
      reorderOffsets[0] = 0;
      for(uint64_t i = 0; i < ownedGIDs.size(); i++){
        reorderOffsets[i+1] = reorderOffsets[i] + reorderDegrees[i];
      }
 
      //printf("--Rank %d: computing reordered adjacencies\n",comm->getRank());
      if(reorderAdjs.size() != globalAdjs.size()){
        //printf("--Rank %d: size mismatch between original and reordered adjs\n"
        //      , comm->getRank());
      }
      //compute reorderAdjs. O(m)
      for(uint64_t i = 0; i < ownedGIDs.size(); i++){
        for(uint64_t j = 0; j < offsets[i+1]-offsets[i]; j++){
          reorderAdjs.at(reorderOffsets[globalToLocal[ownedGIDs[i]]]+j) = 
            globalToLocal[globalAdjs[offsets[i]+j]];
        }
      }
      
      reorderToLocal.resize(ownedGIDs.size());
      for(uint64_t i = 0; i < ownedGIDs.size(); i++){
        //printf("--Rank %d: reordered ID %d is local ID %d\n",comm->getRank(),globalToLocal[ownedGIDs[i]],i);
        reorderToLocal[globalToLocal[ownedGIDs[i]]] = i;
      }
      
      //printf("--Rank %d: done reordering graph\n", comm->getRank());
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
      Kokkos::View<int*, Tpetra::Map<>::device_type >  sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      
      KokkosGraph::Experimental::graph_color_symbolic(&kh, nVtx, nVtx, 
                                                      offset_view, adjs_view);
      Kokkos::fence();
      //std::cout << std::endl << 
      //    "Time: " << 
      //    kh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      //    "Num colors: " << 
      //    kh.get_graph_coloring_handle()->get_num_colors() << " "
      //    "Num Phases: " << 
      //    kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
      //std::cout << "\t"; 
      //KokkosKernels::Impl::print_1Dview(
      //                     kh.get_graph_coloring_handle()->get_vertex_colors());
      
      numColors = kh.get_graph_coloring_handle()->get_num_colors();
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
      
      Teuchos::RCP<Teuchos::Time>
             timeReorder(Teuchos::TimeMonitor::getNewTimer("00 REORDER")),
             timeInterior(Teuchos::TimeMonitor::getNewTimer("01 INTERIOR")),
             timeBoundary(Teuchos::TimeMonitor::getNewTimer("02 BOUNDARY"));
      //this will color the global graph in a manner similar to Zoltan
      
      timeReorder->start();
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
      
      RCP<const map_t> mapOwned;
      RCP<const map_t> mapWithCopies;
      
      std::vector<gno_t> finalGIDs;
      std::vector<offset_t> finalOffset_vec;
      std::vector<lno_t> finalAdjs_vec;      

      lno_t nInterior = 0;
      std::vector<lno_t> reorderToLocal;
      for(int i = 0;  i< nVtx; i++) reorderToLocal.push_back(i);

      //printf("Starting to create local graph\n");
      std::string kokkos_only_interior = pl->get<std::string>("Kokkos_only_interior","false");
      if(comm->getSize() == 1 || kokkos_only_interior=="false") {
        //Set up a typical local mapping here. Need to use the first step of reordering, I think
        std::unordered_map<gno_t,lno_t> globalToLocal;
        std::vector<gno_t> ownedPlusGhosts;
        for (int i = 0; i < vtxIDs.size(); i++){
          globalToLocal[vtxIDs[i]] = i;
          ownedPlusGhosts.push_back(vtxIDs[i]);
        }
        int nghosts = 0;
        for (int i = 0; i < adjs.size(); i++){
          if(globalToLocal.count(adjs[i]) == 0){
            //new unique ghost found
            ownedPlusGhosts.push_back(adjs[i]);
            globalToLocal[adjs[i]] = vtxIDs.size() + nghosts;
            nghosts++;
            
          }
        }
        
        finalAdjs_vec.resize(adjs.size()); 
        for(int i = 0; i < finalAdjs_vec.size();i++){
          //printf("finalAdjs index %d is local vtx %d\n",i,globalToLocal[adjs[i]]);
          finalAdjs_vec[i] = globalToLocal[adjs[i]];
        }
        //finalAdjs = Teuchos::arrayViewFromVector(finalAdjs_vec);
        for(int i = 0; i < offsets.size(); i++) finalOffset_vec.push_back(offsets[i]);
        //finalOffsets = offsets;
        finalGIDs = ownedPlusGhosts;
        Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                             <Tpetra::global_size_t>::invalid();
        mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));
        mapWithCopies = rcp(new map_t(dummy, 
                                  Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                  0, comm)); 
                                      
      } else {
        //printf("Doing the reordering\n");
        //reorder the graph so that boundary vertices are in the
        //end of the offset array.
        finalAdjs_vec.resize(adjs.size());
        finalOffset_vec.resize(offsets.size());
        std::vector<gno_t> reorderGIDs = reorderGraph(vtxIDs,
                                                        adjs,
                                                        offsets,
                                                        finalAdjs_vec,
                                                        finalOffset_vec,
                                                        reorderToLocal,
                                                        nInterior); 
        finalGIDs = reorderGIDs;
        //set up necessary data structures to build the FEMultiVector
        std::vector<gno_t> ownedReorderGIDs;
        for(int i = 0; i < nVtx; i++){
          ownedReorderGIDs.push_back(reorderGIDs[i]);
        }
         
        Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                               <Tpetra::global_size_t>::invalid();
        mapOwned = rcp(new map_t(dummy, Teuchos::arrayViewFromVector(ownedReorderGIDs),
                                             0, comm));
        mapWithCopies = rcp(new map_t(dummy, 
                                        Teuchos::arrayViewFromVector(reorderGIDs),
                                        0, comm));
      }
      //create the FEMultiVector for the distributed communication.
      //We also use the views from this datastructure as arguments to
      //KokkosKernels coloring functions.
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, 
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, 
                                                    importer, 1, true));
      //printf("Done creating local graph\n");
      timeReorder->stop();
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
      for(int i = 0; i < finalGIDs.size(); i++){
        Zoltan_Srand((unsigned int) finalGIDs[i], NULL);
        rand[i] = (int) (((double) Zoltan_Rand(NULL)/(double) ZOLTAN_RAND_MAX)*100000000);
      }

      
      ArrayView<const offset_t> finalOffsets = Teuchos::arrayViewFromVector(finalOffset_vec);
      ArrayView<const lno_t> finalAdjs = Teuchos::arrayViewFromVector(finalAdjs_vec);

      // call coloring function
      hybridGMB(nVtx, nInterior, finalAdjs, finalOffsets,colors,femv,finalGIDs,rand,timeInterior,timeBoundary);
      
      //copy colors to the output array.
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
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
      comm->barrier();
    }
     
    void hybridGMB(const size_t nVtx,lno_t nInterior, Teuchos::ArrayView<const lno_t> adjs, 
                   Teuchos::ArrayView<const offset_t> offsets, 
                   Teuchos::ArrayRCP<int> colors, Teuchos::RCP<femv_t> femv,
                   std::vector<gno_t> reorderGIDs,
                   std::vector<int> rand,
                   Teuchos::RCP<Teuchos::Time> timeInterior,
                   Teuchos::RCP<Teuchos::Time> timeBoundary){
     Teuchos::RCP<Teuchos::Time> 
             timeBoundaryComm(Teuchos::TimeMonitor::getNewTimer("03 BOUNDARY-COMM")),
             timeBoundaryComp(Teuchos::TimeMonitor::getNewTimer("04 BOUNDARY-COMP"));
      
      timeInterior->start();
      //make views out of arrayViews 
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> host_offsets("Host Offset view", offsets.size());
      for(int i = 0; i < offsets.size(); i++){
        host_offsets(i) = offsets[i];
      }
      Kokkos::View<lno_t*, Tpetra::Map<>::device_type> host_adjs("Host Adjacencies view", adjs.size());
      for(int i = 0; i < adjs.size(); i++){
        host_adjs(i) = adjs[i];
      }

      std::string kokkos_only_interior = pl->get<std::string>("Kokkos_only_interior","false");
      size_t kokkosVerts = nVtx;
      if(kokkos_only_interior == "true" && comm->getSize() != 1) {
        kokkosVerts = nInterior;
      }
      //call the KokkosKernels coloring function with the Tpetra default spaces.
      this->colorInterior<Tpetra::Map<>::execution_space,
                          Tpetra::Map<>::memory_space,
                          Tpetra::Map<>::memory_space>
                 (kokkosVerts, host_adjs, host_offsets, colors, femv);
      
      timeInterior->stop();
      //set the batch size to a reasonable default
      timeBoundary->start();
      timeBoundaryComp->start();
      int batch_size = pl->get<int>("Hybrid_batch_size",100);
      
      if(batch_size < 0) batch_size = nVtx;
      //color boundary vertices using FEMultiVector
      //create a queue of vertices to color
      //std::queue<lno_t> recoloringQueue;
      //std::queue<lno_t> conflictQueue;
      typedef Tpetra::Map<>::device_type device_type;
      typedef Tpetra::Map<>::execution_space execution_space;
      typedef Tpetra::Map<>::memory_space memory_space;
      //This is the Kokkos version of two queues. These will attempt to be used in parallel.
      Kokkos::View<lno_t*, device_type> recoloringQueue("recoloringQueue",nVtx);
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
        recoloringQueue(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringQueue_atomic=recoloringQueue;
      Kokkos::View<int[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize;
      Kokkos::View<lno_t*, device_type> conflictQueue("conflictQueue",nVtx);
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
        conflictQueue(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > conflictQueue_atomic = conflictQueue;
      Kokkos::View<int[1], device_type> conflictSize("Conflict Queue Size");
      conflictSize(0) = 0;
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > conflictSize_atomic = conflictSize;
      
      Kokkos::View<int*,device_type> host_rand("randVec",rand.size());
      for(int i = 0; i < rand.size(); i++){
        host_rand(i) = rand[i];
      }
      
      Kokkos::View<int[1], device_type> host_maxDegree("max degree");
      
      //get a subview of the colors:
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
      //bootstrap distributed coloring, add conflicting vertices to the recoloring queue.
      if(kokkos_only_interior == "false"||comm->getSize() > 1){
        timeBoundaryComp->stop();
        timeBoundaryComm->start();
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();
        timeBoundaryComm->stop();
        timeBoundaryComp->start();
        //printf("Adding conflicts to the recoloring queue\n");
        /*for(int i = 0; i < nVtx; i++){
          for(int j = offsets[i]; j < offsets[i+1];j++){
            //printf("getting color for vertex %d\n",i);
            int currColor = femv->getData(0)[i];
            //printf("getting color for vertex %d, at index %d\n",adjs[j],j);
            int nborColor = femv->getData(0)[adjs[j]];
            if(currColor == nborColor && rand[i] > rand[adjs[j]]){
              recoloringQueue.push(i);
              break;
            } 
          }
        }*/
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA (const int& i){
          for(offset_t j = host_offsets(i); j < host_offsets(i+1); j++){
            int currColor = femv_colors(i);
            int nborColor = femv_colors(host_adjs(j));
            if(currColor == nborColor && host_rand(i) > host_rand(host_adjs(j))){
              recoloringQueue_atomic(recoloringSize_atomic(0)++) = i;
              //conflictQueue_atomic(conflictSize_atomic(0)++) = i;
              break;
            }
          }
        });
        //printf("Finished adding conflicts to the recoloring queue\n");
      }
      
      int maxDegree = 0;
      for(int i = 0; i< nVtx; i++){
        if( offsets[i+1]-offsets[i] > maxDegree){
          maxDegree = offsets[i+1] - offsets[i];
        }
      }
      int globalMaxDegree = 0;
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &maxDegree, &globalMaxDegree);
      host_maxDegree(0) = globalMaxDegree;
      Teuchos::ArrayRCP<bool> forbiddenColors;
      forbiddenColors.resize(numColors);
      Kokkos::View<bool*, device_type> forbidden("Forbidden colors",globalMaxDegree);
      Kokkos::parallel_for(globalMaxDegree, KOKKOS_LAMBDA(const int& i){
        forbidden(i) = false;
      });
      
      for(int i = 0; i < numColors; i++) forbiddenColors[i] = false;
      //printf("--Rank %d: batch size: %d\n",comm->getRank(),batch_size);
      int vertsPerRound[100];
      bool done = false; //We're only done when all processors are done
      if(comm->getSize() == 1) done = true;
      int distributedRounds = 0; //this is the same across all processors
      int i = nInterior; //we only need to worry about boundary vertices here.
      //  while the queue is not empty
          while(recoloringSize(0) > 0 || !done){
            //printf("-- Rank %d: Recoloring %d verts\n",comm->getRank(),recoloringQueue.size());
        //  add next batch to the queue
            if (kokkos_only_interior == "true"){
              for(size_t j = i; j < i+batch_size; j++){
                if(j < nVtx) {
                  //printf("--Rank %d: pushing %d on the queue\n",comm->getRank(),j);
                  //recoloringQueue.push(j); 
                  recoloringQueue(recoloringSize(0)++) = j;
                }
              }
              if(i < nVtx) i+= batch_size;
            }else i = nVtx;
      //    color everything in the recoloring queue, put everything on conflict queue
            if(distributedRounds < 100) vertsPerRound[distributedRounds] = recoloringSize(0);
            //eventually turn this to a Kokkos::parallel_for (could make more conflicts, but eh?
            /*while(recoloringSize(0) > 0) {
               //printf("--Rank %d: coloring vertex %u\n",comm->getRank(),recoloringQueue.front());
               lno_t currVtx = recoloringQueue(recoloringSize(0)-1);
               //recoloringQueue.pop();
               recoloringSize(0)--;
               conflictQueue(conflictSize(0)++) = currVtx;
               for(offset_t nborIdx = offsets[currVtx]; nborIdx < offsets[currVtx+1]; nborIdx++){
                 
                 int nborColor = femv_colors(adjs[nborIdx]);
                 if(nborColor > 0){
                   if(nborColor > numColors){
                     forbiddenColors.resize(nborColor);
                     numColors = nborColor;
                   }
                   forbidden(nborColor-1) = true;
                 } 
               }
               //pick the color for currVtx based on the forbidden array
               bool colored = false;
               <for(int i = 0; i < maxDegree; i++){
                 if(!forbidden(i)) {
                   femv->replaceLocalValue(currVtx,0,i+1);
                   colored = true;
                   break;
                 }
               }
               Kokkos::parallel_for(maxDegree,KOKKOS_LAMBDA(const int& i){
                 forbidden(i) = false;
               });
               if(!colored){
                 femv->replaceLocalValue(currVtx,0,numColors+1);
	 numColors++;
	 forbiddenColors.resize(numColors);
               }
               //printf("--Rank %d: colored vtx %u color %d\n",comm->getRank(),currVtx, femv->getData(0)[currVtx]);
            }*/
            //printf("--Rank %d: recoloring\n\t",comm->getRank());
            //for(int i = 0; i< recoloringSize(0); i++){
              //printf("%d ",recoloringQueue(i));
            //}
            //printf("\n");
            typedef KokkosKernels::Experimental::KokkosKernelsHandle
                <size_t, lno_t, lno_t, execution_space, memory_space, 
                 memory_space> KernelHandle;
            //KernelHandle kh;
            //kh->create_graph_handle(KokkosGraph::COLORING_VB);
            /*typename KokkosGraph::Impl::GraphColor_VB<typename KernelHandle::GraphColoringHandleType,
                                             Kokkos::View<offset_t*,device_type>,
                                             Kokkos::View<lno_t*, device_type> >
                       ::functorGreedyColor gc(nVtx,host_offsets,host_adjs,femv_colors,
                                                     recoloringQueue,recoloringSize(0),0);
            Kokkos::parallel_for("Recoloring",recoloringSize(0),gc);*/
            //need to provide means of expanding forbidden array.

            typedef Kokkos::TeamPolicy<execution_space>::member_type member_type;
            Kokkos::parallel_for(Kokkos::TeamPolicy<>(recoloringSize(0),1).set_scratch_size(0,Kokkos::PerTeam(0),Kokkos::PerThread(2*host_maxDegree(0)*sizeof(bool))),
                    KOKKOS_LAMBDA (const member_type& team_member){
              bool* dev_forbidden = (bool*) team_member.team_shmem().get_shmem(host_maxDegree(0)*sizeof(bool));
              for(int j=0; j < host_maxDegree(0); j++) dev_forbidden[j] = false;
              lno_t currVtx = recoloringQueue(team_member.league_rank());
              conflictQueue_atomic(conflictSize_atomic(0)++) = currVtx;
              for(offset_t nborIdx = host_offsets(currVtx); nborIdx < host_offsets(currVtx+1); nborIdx++){
                int nborColor = femv_colors(host_adjs(nborIdx));
                if(nborColor > 0) dev_forbidden[nborColor-1] = true;
              }
              for( int j = 0; j < 2*host_maxDegree(0); j++){
                if(!dev_forbidden[j]){
                  femv_colors(currVtx)=j+1;
                  break;
                }
              }
            });
            /*Kokkos::parallel_for(recoloringSize(0), KOKKOS_LAMBDA(const int& i){
              bool* dev_forbidden = new bool[2*host_maxDegree(0)];
              for(int j = 0; j < 2*host_maxDegree(0); j++) dev_forbidden[j]=false;
              lno_t currVtx = recoloringQueue(i);
              conflictQueue_atomic(conflictSize_atomic(0)++) = currVtx;
              for(offset_t nborIdx = host_offsets(currVtx); nborIdx < host_offsets(currVtx+1); nborIdx++){
                int nborColor=femv_colors(host_adjs(nborIdx));
                if(nborColor > 0) dev_forbidden[nborColor-1] = true;
              }
              for(int j = 0; j < 2*host_maxDegree(0); j++){
                if(!dev_forbidden[j]){
                  femv_colors(currVtx) = j+1;
                  //printf("\t--Rank %d recolored vertex %u color %d\n",comm->getRank(),currVtx,j+1);
                }
              }
              delete [] dev_forbidden;
            });*/
            recoloringSize(0) = 0;
            timeBoundaryComp->stop();
            timeBoundaryComm->start();
      //    communicate
            //printf("--Rank %d: communicating\n",comm->getRank());
            femv->switchActiveMultiVector();
            femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
            femv->switchActiveMultiVector();
            
            timeBoundaryComm->stop();
            timeBoundaryComp->start();
      //    check for conflicts with the vertices that were just colored
            /*while(conflictQueue.size() > 0){
      //      if conflicts were detected put them on the recoloring queue
              lno_t currVtx = conflictQueue.front();
              conflictQueue.pop();
              bool conflict = false;
              int currColor = femv->getData(0)[currVtx];
              for(offset_t nborIdx = offsets[currVtx]; nborIdx < offsets[currVtx+1]; nborIdx++){
                int nborColor = femv->getData(0)[adjs[nborIdx]];
                if(reorderGIDs[currVtx] == 117961 && reorderGIDs[adjs[nborIdx]] == 150790){
                  printf("--Rank %d: checking 117961 (color %d) against 150790 (color %d)\n",comm->getRank(),
                          femv->getData(0)[currVtx], femv->getData(0)[adjs[nborIdx]]);
                }
                if(nborColor == currColor && rand[currVtx] > rand[adjs[nborIdx]]) {
                  //printf("--Rank %d: vtx %u conflicts with vtx %u\n",comm->getRank(),currVtx,adjs[nborIdx]);
                  conflict = true;
                }
              }
              if(conflict) {
                recoloringQueue.push(currVtx);
                //printf("--Rank %d: putting vertex %u on the recoloring queue\n",comm->getRank(),currVtx);
              }
            }*/
            //printf("--Rank %d: Detecting conflicts\n",comm->getRank());
            Kokkos::parallel_for(conflictSize(0),KOKKOS_LAMBDA(const int& i){
              lno_t currVtx = conflictQueue(i);
              bool conflict = false;
              int currColor = femv_colors(currVtx);
              for(offset_t nborIdx = host_offsets(currVtx); nborIdx < host_offsets(currVtx+1); nborIdx++){
                int nborColor = femv_colors(host_adjs(nborIdx));
                if(nborColor == currColor && host_rand(currVtx) > host_rand(host_adjs(nborIdx))){
                  conflict = true;
                  break;
                }
              }
              if(conflict) recoloringQueue_atomic(recoloringSize_atomic(0)++) = currVtx;
            });
            conflictSize(0) = 0;
            //do a reduction to determine if we're done
            int globalDone = 0;
            int localDone = 0;
            if(kokkos_only_interior=="true") localDone= recoloringSize(0) + (nVtx > i);
            else localDone = recoloringSize(0);
            Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
            //We're only allowed to stop once everyone has no work to do.
            //collectives will hang if one process exits. 
            printf("globalDone = %d\n",globalDone);
            distributedRounds++;
            done = !globalDone;
          }
        
      timeBoundaryComp->stop();
      timeBoundary->stop();
      //color interior vertices if not colored yet. 
      //(must alter Kokkos-Kernels to allow partial colorings)
      //As long as the initial coloring is untouched, we can pass the whole 
      //graph to the kokkos-kernels coloring.
      
      //print the coloring
      //for(int i = 0; i < femv->getData(0).size(); i++){
      //  printf("--Rank %d: local vtx %u is color %d\n",comm->getRank(),i,femv->getData(0)[i]);
      //}
      
      //print how many rounds of speculating/correcting happened (this should be the same for all ranks):
      if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n", distributedRounds);
      int totalVertsPerRound[100];
      for(int i = 0; i < 100; i++) totalVertsPerRound[i]=0;
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,100,vertsPerRound,totalVertsPerRound);
      for(int i = 0; i < std::min(distributedRounds,100); i++){
        if(comm->getRank()==0) printf("recolored %d vertices in round %d\n",totalVertsPerRound[i],i);
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
