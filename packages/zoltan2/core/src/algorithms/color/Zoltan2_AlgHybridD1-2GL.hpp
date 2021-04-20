#ifndef _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_
#define _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <sys/time.h>

#include "Zoltan2_Algorithm.hpp"
#include "Zoltan2_GraphModel.hpp"
#include "Zoltan2_ColoringSolution.hpp"
#include "Zoltan2_Util.hpp"
#include "Zoltan2_TPLTraits.hpp"
#include "Zoltan2_AlltoAll.hpp"
#include "Zoltan2_AlgHybrid2GL.hpp"

#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"

/////////////////////////////////////////////////
//! \file Zoltan2_Distance1_2GhostLayer.hpp
//! \brief A Communication Avoidant Distance-1 Coloring Algorithm


namespace Zoltan2{

template <typename Adapter>
class AlgDistance1TwoGhostLayer : public AlgTwoGhostLayer<Adapter> {
  public:
    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t,gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>;
    using device_type = Tpetra::Map<>::device_type;
    using execution_space = Tpetra::Map<>::execution_space;
    using memory_space = Tpetra::Map<>::memory_space;

  private:
    
    virtual void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,device_type > adjs_view,
                       Kokkos::View<offset_t*,device_type > offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, device_type > vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){

      //setup types to be used by KokkosKernels
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, execution_space, memory_space, memory_space>;
      using lno_row_view_t = Kokkos::View<offset_t*, device_type>;
      using lno_nnz_view_t = Kokkos::View<lno_t*, device_type>;
      KernelHandle kh;
  
      //pick which KokkosKernels algorithm to use, based on max graph degree
      if(recolor){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
      }
      //vertex_list_size indicates whether we have provided a list of vertices to recolor.
      //Only makes a difference if the algorithm to be used is VBBIT
      if(vertex_list_size != 0){
        kh.get_graph_coloring_handle()->set_vertex_list(vertex_list,vertex_list_size);
      }
  
      kh.set_verbose(this->verbose);
  
      //set the initial coloring to be the colors from femv.
      auto femvColors = femv->template getLocalView<memory_space>();
      auto sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      kh.get_graph_coloring_handle()->set_tictoc(this->verbose);
  
      KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t> (&kh, nVtx,nVtx, offset_view, adjs_view);

      this->numColors = kh.get_graph_coloring_handle()->get_num_colors();

      if(this->verbose){
        std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
      }
    }
   
    virtual void colorInterior_serial(const size_t nVtx,
                       typename Kokkos::View<lno_t*, device_type >::HostMirror adjs_view,
                       typename Kokkos::View<offset_t*,device_type >::HostMirror offset_view,
                       Teuchos::RCP<femv_t> femv,
                       typename Kokkos::View<lno_t*, device_type>::HostMirror vertex_list,
                       size_t vertex_list_size = 0,
                       bool recolor=false) {
        //fill in colorInterior_serial
	using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
            <offset_t, lno_t, lno_t, Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, Kokkos::HostSpace>;
	using lno_row_view_t = Kokkos::View<offset_t*, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>>;
	using lno_nnz_view_t = Kokkos::View<lno_t*, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>>;
	KernelHandle kh;

	if(recolor){
	  kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
	} else {
	  kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
	}

	if(vertex_list_size != 0){
	  kh.get_graph_coloring_handle()->set_vertex_list(vertex_list, vertex_list_size);
	}

	kh.set_verbose(this->verbose);

	auto femvColors = femv->getLocalViewHost();
	auto sv = subview(femvColors, Kokkos::ALL, 0);
	kh.get_graph_coloring_handle()->set_vertex_colors(sv);
	kh.get_graph_coloring_handle()->set_tictoc(this->verbose);

	KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t>(&kh, nVtx, nVtx, offset_view, adjs_view);

	this->numColors = kh.get_graph_coloring_handle()->get_num_colors();

	if(this->verbose){
	  std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
	}
    }
  public:
    virtual void detectConflicts(const size_t n_local,
                                 Kokkos::View<offset_t*, device_type > dist_offsets_dev,
                                 Kokkos::View<lno_t*, device_type > dist_adjs_dev,
                                 Kokkos::View<int*,device_type > femv_colors,
                                 Kokkos::View<lno_t*, device_type > boundary_verts_view,
                                 gno_t boundary_size,
                                 Kokkos::View<lno_t*,
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_recolor_atomic,
                                 Kokkos::View<int[1],
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_size_atomic,
                                 Kokkos::View<lno_t*,
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_atomic,
                                 Kokkos::View<size_t[1],
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic,
				 Kokkos::View<gno_t[1], device_type> recoloringSize,
                                 Kokkos::View<int*,
                                              device_type> rand,
                                 Kokkos::View<gno_t*,
                                              device_type> gid,
                                 Kokkos::View<gno_t*,
                                              device_type> ghost_degrees,
				 bool recolor_degrees){
      Kokkos::RangePolicy<execution_space> policy(n_local,rand.size());
      Kokkos::parallel_reduce("conflict detection",policy, KOKKOS_LAMBDA (const int& i, gno_t& recoloring_size){
        lno_t localIdx = i;
        int currColor = femv_colors(localIdx);
        int currDegree = ghost_degrees(i-n_local);
        for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
          int nborColor = femv_colors(dist_adjs_dev(j));
          int nborDegree = 0;
          if((size_t)dist_adjs_dev(j) < n_local) nborDegree = dist_offsets_dev(dist_adjs_dev(j)+1) - dist_offsets_dev(dist_adjs_dev(j));
          else nborDegree = ghost_degrees(dist_adjs_dev(j) - n_local);
          if(currColor == nborColor ){
            if(currDegree < nborDegree && recolor_degrees){
              femv_colors(localIdx) = 0;
              recoloring_size++;
              break;
            }else if(nborDegree < currDegree && recolor_degrees){
              femv_colors(dist_adjs_dev(j)) = 0;
              recoloring_size++;
            }else if(rand(localIdx) > rand(dist_adjs_dev(j))){
              recoloring_size++;
              femv_colors(localIdx) = 0;
              break;
            }else if(rand(dist_adjs_dev(j)) > rand(localIdx)){
                recoloring_size++;
                femv_colors(dist_adjs_dev(j)) = 0;
            } else {
              if (gid(localIdx) >= gid(dist_adjs_dev(j))){
                femv_colors(localIdx) = 0;
                recoloring_size++;
                break;
              } else {
                  femv_colors(dist_adjs_dev(j)) = 0;
                  recoloring_size++;
              }
            }
          }
        }
      },recoloringSize(0));
      Kokkos::fence();
      Kokkos::parallel_for(femv_colors.size(), KOKKOS_LAMBDA (const size_t& i){
        if(femv_colors(i) == 0){
          if(i < n_local){
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
          }
          verts_to_recolor_atomic(verts_to_recolor_size_atomic(0)++) = i;
        }
      });
      Kokkos::fence();
    }

    virtual void detectConflicts_serial(const size_t n_local,
                                 typename Kokkos::View<offset_t*, device_type >::HostMirror dist_offsets_host,
                                 typename Kokkos::View<lno_t*, device_type >::HostMirror dist_adjs_host,
                                 typename Kokkos::View<int*,device_type >::HostMirror femv_colors,
                                 typename Kokkos::View<lno_t*, device_type >::HostMirror boundary_verts_view,
                                 gno_t boundary_size,
                                 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_recolor,
                                 typename Kokkos::View<int[1],device_type>::HostMirror verts_to_recolor_size,
                                 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_send,
                                 typename Kokkos::View<size_t[1],device_type>::HostMirror verts_to_send_size,
				 typename Kokkos::View<gno_t[1], device_type>::HostMirror recoloringSize,
                                 typename Kokkos::View<int*,  device_type>::HostMirror rand,
                                 typename Kokkos::View<gno_t*,device_type>::HostMirror gid,
                                 typename Kokkos::View<gno_t*,device_type>::HostMirror ghost_degrees,
				 bool recolor_degrees) {
        //fill in detectConflicts_serial
	for(size_t i = n_local; i < rand.size(); i++){
          lno_t localIdx = i;
          int currColor = femv_colors(localIdx);
          int currDegree = ghost_degrees(i-n_local);
          for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
            int nborColor = femv_colors(dist_adjs_host(j));
            int nborDegree = 0;
            if((size_t)dist_adjs_host(j) < n_local){
              nborDegree = dist_offsets_host(dist_adjs_host(j)+1) - dist_offsets_host(dist_adjs_host(j));
            } else {
              nborDegree = ghost_degrees(dist_adjs_host(j) - n_local);
            }
            if(currColor == nborColor) {
              if(currDegree < nborDegree && recolor_degrees){
                femv_colors(localIdx) = 0;
                recoloringSize(0)++;
                break;
              }else if(nborDegree < currDegree && recolor_degrees){
                femv_colors(dist_adjs_host(j)) = 0;
                recoloringSize(0)++;
              }else if(rand(localIdx) > rand(dist_adjs_host(j))){
                femv_colors(localIdx) = 0;
                recoloringSize(0)++;
                break;
              } else if (rand(dist_adjs_host(j)) > rand(localIdx)) {
                femv_colors(dist_adjs_host(j)) = 0;
                recoloringSize(0)++;
              } else {
                if(gid(localIdx) >= gid(dist_adjs_host(j))){
                  femv_colors(localIdx) = 0;
                  recoloringSize(0)++;
                  break;
                } else {
                  femv_colors(dist_adjs_host(j)) = 0;
                  recoloringSize(0)++;
                }
              }
            }
          }
        }
        for(size_t i = 0; i < femv_colors.size(); i++){
          if(femv_colors(i) == 0){
            if(i < n_local){
              verts_to_send(verts_to_send_size(0)++) = i;
            }
            verts_to_recolor(verts_to_recolor_size(0)++) = i;
          }
        }

    }

    virtual void constructBoundary(const size_t n_local,
                                   Kokkos::View<offset_t*, device_type> dist_offsets_dev,
                                   Kokkos::View<lno_t*, device_type> dist_adjs_dev,
				   typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host,
				   typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host,
                                   Kokkos::View<lno_t*, device_type>& boundary_verts,
                                   gno_t& boundary_size,
                                   Kokkos::View<lno_t*,
                                                device_type,
                                                Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_atomic,
                                   Kokkos::View<size_t[1],
                                                device_type,
                                                Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic){

      Kokkos::parallel_for(n_local, KOKKOS_LAMBDA(const int& i){
        for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
          if((size_t)dist_adjs_dev(j) >= n_local){
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
            break;
          }
          bool found = false;
          for(offset_t k = dist_offsets_dev(dist_adjs_dev(j)); k < dist_offsets_dev(dist_adjs_dev(j)+1); k++){
            if((size_t)dist_adjs_dev(k) >= n_local){
              verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
              found = true;
              break;
            }
          }
          if(found) break;
        }
      });
      Kokkos::fence();
      boundary_size = verts_to_send_size_atomic(0);
    }


  public:
    AlgDistance1TwoGhostLayer(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : AlgTwoGhostLayer<Adapter>(adapter_,pl_,env_,comm_){
    }
    

    
}; //end class



}//end namespace Zoltan2

#endif
