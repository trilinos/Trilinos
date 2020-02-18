#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <unistd.h>

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

    //function to invoke KokkosKernels distance-1 coloring    
    template <class ExecutionSpace, typename TempMemorySpace, 
              typename MemorySpace>
    void colorInterior(const size_t nVtx, 
                       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > adjs_view,
                       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > offset_view, 
                       Teuchos::ArrayRCP<int> colors,
                       Teuchos::RCP<femv_t> femv,
                       bool recolor=false){
      
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
    
      if(recolor){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
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
      
      numColors = kh.get_graph_coloring_handle()->get_num_colors();
    }
    
    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    int numColors;
    
  public:
    //constructor for the  hybrid distributed distance-1 algorithm
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
      int rank = comm->getRank(); 
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
      for(size_t i = 0;  i< nVtx; i++) reorderToLocal.push_back(i);
      
      //Set up a typical local mapping here.
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
      for(size_t i = 0; i < finalAdjs_vec.size();i++){
        finalAdjs_vec[i] = globalToLocal[adjs[i]];
      }
      for(int i = 0; i < offsets.size(); i++) finalOffset_vec.push_back(offsets[i]);
      finalGIDs = ownedPlusGhosts;

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
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, 
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, 
                                                    importer, 1, true));
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
      for(size_t i = 0; i < finalGIDs.size(); i++){
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
     
      //RESULT REPORTING INSTRUMENTATION
      //Teuchos::TimeMonitor::summarize();
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
      typedef Tpetra::Map<>::device_type device_type;
      typedef Tpetra::Map<>::execution_space execution_space;
      typedef Tpetra::Map<>::memory_space memory_space;
     
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> dist_degrees("Owned+Ghost degree view",rand.size());
      for(int i = 0; i < adjs.size(); i++){
        dist_degrees(adjs[i])++;
      }
      for(int i = 0; i < offsets.size()-1; i++){
        dist_degrees(i) = offsets[i+1] - offsets[i];
      }
      Kokkos::View<offset_t*, Tpetra::Map<>::device_type> dist_offsets("Owned+Ghost Offset view", rand.size()+1);
      dist_offsets(0) = 0;
      int total_adjs = adjs.size();
      for(size_t i = 1; i < rand.size()+1; i++){
        dist_offsets(i) = dist_degrees(i-1) + dist_offsets(i-1);
        total_adjs+= dist_degrees(i-1);
      }
      Kokkos::View<lno_t*, Tpetra::Map<>::device_type> dist_adjs("Owned+Ghost adjacency view", total_adjs);
      for(size_t i = 0; i < rand.size(); i++){
        dist_degrees(i) = 0;
      }
      for(int i = 0; i < adjs.size(); i++) dist_adjs(i) = adjs[i];
      for(size_t i = 0; i < nVtx; i++){
        for(size_t j = offsets[i]; j < offsets[i+1]; j++){
          if( (size_t)adjs[j] >= nVtx){
            dist_adjs(dist_offsets(adjs[j]) + dist_degrees(adjs[j])) = i;
            dist_degrees(adjs[j])++;
          }
        }
      }
      //This is the Kokkos version of two queues. These will attempt to be used in parallel.
      Kokkos::View<lno_t*, device_type> recoloringQueue("recoloringQueue",nVtx);
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
        recoloringQueue(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringQueue_atomic=recoloringQueue;
      Kokkos::View<int[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize; 
      Kokkos::View<int*,device_type> host_rand("randVec",rand.size());
      for(size_t i = 0; i < rand.size(); i++){
        host_rand(i) = rand[i];
      }
      
      
      //bootstrap distributed coloring, add conflicting vertices to the recoloring queue.
      if(comm->getSize() > 1){
        timeBoundaryComp->stop();
        timeBoundaryComm->start();
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();
        timeBoundaryComm->stop();
        timeBoundaryComp->start();
        //get a subview of the colors:
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        
        //detect conflicts from the initial coloring
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA (const int& i){
          for(offset_t j = host_offsets(i); j < host_offsets(i+1); j++){
            int currColor = femv_colors(i);
            int nborColor = femv_colors(host_adjs(j));
            if(currColor == nborColor ){
              if(host_rand(i) > host_rand(host_adjs(j))){
                femv_colors(i) = 0;
                recoloringSize_atomic(0)++;
                break;
              }
              if(host_rand(host_adjs(j)) > host_rand(i)){
                femv_colors(host_adjs(j)) = 0;
                recoloringSize_atomic(0)++;
              }
            }
          }
        });
        //ensure the parallel_for has completed before continuing.
        Kokkos::fence();
      }
      
      
      int vertsPerRound[100];
      bool done = false; //We're only done when all processors are done
      if(comm->getSize() == 1) done = true;
      int distributedRounds = 0; //this is the same across all processors
      //while the queue is not empty
      while(recoloringSize(0) > 0 || !done){
        //color everything in the recoloring queue, put everything on conflict queue
        if(distributedRounds < 100) vertsPerRound[distributedRounds] = recoloringSize(0);
        //get a subview of the colors:
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        
        //use KokkosKernels to recolor the conflicting vertices.  
        this->colorInterior<Tpetra::Map<>::execution_space,
                            Tpetra::Map<>::memory_space,
                            Tpetra::Map<>::memory_space>
                            (femv_colors.size(),dist_adjs,dist_offsets,colors,femv,true);
            
        recoloringSize(0) = 0;
        timeBoundaryComp->stop();
        timeBoundaryComm->start();
        //communicate
        femv->switchActiveMultiVector();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        femv->switchActiveMultiVector();
            
        timeBoundaryComm->stop();
        timeBoundaryComp->start();
        
        //detect conflicts in parallel. For a detected conflict,
        //reset the vertex-to-be-recolored's color to 0, in order to
        //allow KokkosKernels to recolor correctly.
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA (const int& i){
          for(offset_t j = host_offsets(i); j < host_offsets(i+1); j++){
            int currColor = femv_colors(i);
            int nborColor = femv_colors(host_adjs(j));
            if(currColor == nborColor ){
              if(host_rand(i) > host_rand(host_adjs(j))){
                femv_colors(i) = 0;
                recoloringSize_atomic(0)++;
                break;
              }
              if(host_rand(host_adjs(j)) > host_rand(i)){
                femv_colors(host_adjs(j)) = 0;
                recoloringSize_atomic(0)++;
              }
            }
          }
        });
        //For Cuda, this fence is necessary to ensure the Kokkos::parallel_for is finished
        //before continuing with the coloring. 
        Kokkos::fence();
        
        //do a reduction to determine if we're done
        int globalDone = 0;
        int localDone = recoloringSize(0);
        Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
        //We're only allowed to stop once everyone has no work to do.
        //collectives will hang if one process exits. 
        distributedRounds++;
        done = !globalDone;
      }
        
      timeBoundaryComp->stop();
      timeBoundary->stop();
      
      
      //print how many rounds of speculating/correcting happened (this should be the same for all ranks):
      /* RESULT REPORTING INSTRUMENTATION
      if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n", distributedRounds);
      int totalVertsPerRound[100];
      for(int i = 0; i < 100; i++) totalVertsPerRound[i]=0;
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,100,vertsPerRound,totalVertsPerRound);
      for(int i = 0; i < std::min(distributedRounds,100); i++){
        if(comm->getRank()==0) printf("recolored %d vertices in round %d\n",totalVertsPerRound[i],i);
      }*/
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
