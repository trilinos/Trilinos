#ifndef _ZOLTAN2_DISTANCE2_HPP_
#define _ZOLTAN2_DISTANCE2_HPP_

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


#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "Kokkos_Core.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance2Color.hpp"
#include "KokkosGraph_Distance2ColorHandle.hpp"

/////////////////////////////////////////////////
//! \file Zoltan2_Distance1_2GhostLayer.hpp
//! \brief A hybrid Distance-2 Coloring Algorithm


namespace Zoltan2{

template <typename Adapter>
class AlgDistance2 : public AlgTwoGhostLayer<Adapter> {

  public:
    
    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t,gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>; 
    using device_type = Tpetra::Map<>::device_type;//Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>;
    using execution_space = Tpetra::Map<>::execution_space;//Kokkos::Cuda;
    using memory_space = Tpetra::Map<>::memory_space;//Kokkos::Cuda::memory_space;
  private:
    
    virtual void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,device_type> adjs_view,
                       Kokkos::View<offset_t*, device_type> offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, device_type> vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, execution_space, memory_space, memory_space>;
  
      KernelHandle kh;
      kh.set_verbose(this->verbose); 
      if(recolor)kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_VB_BIT);
      else kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_NB_BIT);
      if(vertex_list_size != 0){
        kh.get_distance2_graph_coloring_handle()->set_vertex_list(vertex_list, vertex_list_size);
      }
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, Kokkos::Device<execution_space,memory_space>> sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_distance2_graph_coloring_handle()->set_vertex_colors(sv);

      KokkosGraph::Experimental::graph_color_distance2(&kh, nVtx, offset_view, adjs_view);
      this->numColors = kh.get_distance2_graph_coloring_handle()->get_num_colors();
      if(this->verbose) std::cout<<"\nKokkosKernels Coloring: "<<kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time()<<"\n";
    
    }
    
    virtual void colorInterior_serial(const size_t nVtx,
		       typename Kokkos::View<lno_t*, device_type>::HostMirror adjs_view,
		       typename Kokkos::View<offset_t*, device_type>::HostMirror offset_view,
		       Teuchos::RCP<femv_t> femv,
		       typename Kokkos::View<lno_t*, device_type>::HostMirror vertex_list,
		       size_t vertex_list_size = 0,
		       bool recolor=false){
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, Kokkos::DefaultHostExecutionSpace, memory_space, memory_space>;

      KernelHandle kh;
      kh.set_verbose(this->verbose);
      if(recolor)kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_VB_BIT);
      else kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_NB_BIT);
      if(vertex_list_size != 0){
        kh.get_distance2_graph_coloring_handle()->set_vertex_list(vertex_list, vertex_list_size);
      }
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, Kokkos::Device<Kokkos::DefaultHostExecutionSpace,memory_space>> sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_distance2_graph_coloring_handle()->set_vertex_colors(sv);

      KokkosGraph::Experimental::graph_color_distance2(&kh, nVtx, offset_view, adjs_view);
      this->numColors = kh.get_distance2_graph_coloring_handle()->get_num_colors();
      if(this->verbose) std::cout<<"\nKokkosKernels Coloring: "<<kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time()<<"\n";

    }

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
    //fill in detectConflicts
    Kokkos::parallel_reduce(boundary_size, KOKKOS_LAMBDA(const uint64_t& i,gno_t& recoloring_size){
        const size_t curr_lid = boundary_verts_view(i);
        const int curr_color = femv_colors(curr_lid);
        const size_t vid_d1_adj_begin = dist_offsets_dev(curr_lid);
        const size_t vid_d1_adj_end   = dist_offsets_dev(curr_lid+1);
        const size_t curr_degree = vid_d1_adj_end - vid_d1_adj_begin;
        for(size_t vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++){
          size_t vid_d1 = dist_adjs_dev(vid_d1_adj);
          size_t vid_d1_degree = 0;
          if(vid_d1 < n_local){
            vid_d1_degree = dist_offsets_dev(vid_d1+1) - dist_offsets_dev(vid_d1);
          } else {
            vid_d1_degree = ghost_degrees(vid_d1-n_local);
          }
          if( vid_d1 != curr_lid && femv_colors(vid_d1) == curr_color){
            if(curr_degree < vid_d1_degree && recolor_degrees){
              femv_colors(curr_lid) = 0;
              recoloring_size++;
              break;
            } else if (vid_d1_degree < curr_degree && recolor_degrees){
              femv_colors(vid_d1) = 0;
              recoloring_size++;
            } else if(rand(curr_lid) < rand(vid_d1)){
              femv_colors(curr_lid) = 0;
              recoloring_size++;
              break;
            } else if(rand(vid_d1) < rand(curr_lid)){
              femv_colors(vid_d1) = 0;
              recoloring_size++;
            } else{
              if(gid(curr_lid) >= gid(vid_d1)){
                femv_colors(curr_lid) = 0;
                recoloring_size++;
                break;
              } else {
                femv_colors(vid_d1) = 0;
                recoloring_size++;
              }
            }
          }
          size_t d2_adj_begin = 0;
          size_t d2_adj_end   = 0;
          d2_adj_begin = dist_offsets_dev(vid_d1);
          d2_adj_end   = dist_offsets_dev(vid_d1+1);
          bool found = false;
          for(size_t vid_d2_adj = d2_adj_begin; vid_d2_adj < d2_adj_end; vid_d2_adj++){
            const size_t vid_d2 = dist_adjs_dev(vid_d2_adj);
            size_t vid_d2_degree = 0;
            if(vid_d2 < n_local){
              vid_d2_degree = dist_offsets_dev(vid_d2+1) - dist_offsets_dev(vid_d2);
            } else {
              vid_d2_degree = ghost_degrees(vid_d2-n_local);
            }
            if(curr_lid != vid_d2 && femv_colors(vid_d2) == curr_color){
              if(curr_degree < vid_d2_degree && recolor_degrees){
                found = true;
                femv_colors(curr_lid) = 0;
                recoloring_size++;
                break;
              } else if(vid_d2_degree < curr_degree && recolor_degrees){
                femv_colors(vid_d2) = 0;
                recoloring_size++;
              } else if(rand(curr_lid) < rand(vid_d2)){
                found = true;
                femv_colors(curr_lid) = 0;
                recoloring_size++;
                break;
              } else if(rand(vid_d2) < rand(curr_lid)){
                femv_colors(vid_d2) = 0;
                recoloring_size++;
              } else {
                if(gid(curr_lid) >= gid(vid_d2)){
                  found = true;
                  femv_colors(curr_lid) = 0;
                  recoloring_size++;
                  break;
                } else {
                  femv_colors(vid_d2) = 0;
                  recoloring_size++;
                }
              }
            }
          }
          if(found) break;
        }
      },recoloringSize(0));
      Kokkos::fence();
      verts_to_send_size_atomic(0) = 0;
      verts_to_recolor_size_atomic(0) = 0;
      Kokkos::parallel_for(femv_colors.size(), KOKKOS_LAMBDA(const uint64_t& i){
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
    
    for(int i = 0; i <boundary_size; i++){
      const size_t curr_lid = boundary_verts_view(i);
      const int curr_color = femv_colors(curr_lid);
      const size_t vid_d1_adj_begin = dist_offsets_host(curr_lid);
      const size_t vid_d1_adj_end = dist_offsets_host(curr_lid+1);
      const size_t curr_degree = vid_d1_adj_end - vid_d1_adj_begin;
      for(size_t vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++){
        size_t vid_d1 = dist_adjs_host(vid_d1_adj);
        size_t vid_d1_degree = 0;
        if(vid_d1 < n_local){
          vid_d1_degree = dist_offsets_host(vid_d1+1) - dist_offsets_host(vid_d1);
        } else {
          vid_d1_degree = ghost_degrees(vid_d1-n_local);
        }
        if(vid_d1 != curr_lid && femv_colors(vid_d1) == curr_color){
              if(curr_degree < vid_d1_degree && recolor_degrees){
                femv_colors(curr_lid) = 0;
                recoloringSize(0)++;
                break;
              } else if (vid_d1_degree < curr_degree && recolor_degrees){
                femv_colors(vid_d1) = 0;
                recoloringSize(0)++;
              } else if(rand(curr_lid) < rand(vid_d1)){
                femv_colors(curr_lid) = 0;
                recoloringSize(0)++;
                break;
              } else if(rand(vid_d1) < rand(curr_lid)){
                femv_colors(vid_d1) = 0;
                recoloringSize(0)++;
              } else {
                if(gid(curr_lid) >= gid(vid_d1)){
                  femv_colors(curr_lid) = 0;
                  recoloringSize(0)++;
                  break;
                } else {
                  femv_colors(vid_d1) = 0;
                  recoloringSize(0)++;
                }
              }
            }
            size_t d2_adj_begin = dist_offsets_host(vid_d1);
            size_t d2_adj_end = dist_offsets_host(vid_d1+1);
            bool found=false;
            for(size_t vid_d2_adj = d2_adj_begin; vid_d2_adj < d2_adj_end; vid_d2_adj++){
              const size_t vid_d2 = dist_adjs_host(vid_d2_adj);
              size_t vid_d2_degree = 0;
              if(vid_d2 < n_local){
                vid_d2_degree = dist_offsets_host(vid_d2+1) - dist_offsets_host(vid_d2);
              } else {
                vid_d2_degree = ghost_degrees(vid_d2-n_local);
              }
              if(curr_lid != vid_d2 && femv_colors(vid_d2) == curr_color){
                if(curr_degree < vid_d2_degree && recolor_degrees){
                  found = true;
                  femv_colors(curr_lid) = 0;
                  recoloringSize(0)++;
                  break;
                } else if(vid_d2_degree < curr_degree && recolor_degrees){
                  femv_colors(vid_d2) = 0;
                } else if(rand(curr_lid) < rand(vid_d2)){
                  found = true;
                  femv_colors(curr_lid) = 0;
                  recoloringSize(0)++;
                  break;
                } else if(rand(vid_d2) < rand(curr_lid)){
                  femv_colors(vid_d2) = 0;
                  recoloringSize(0)++;
                } else {
                  if(gid(curr_lid) >= gid(vid_d2)){
                    found = true;
                    femv_colors(curr_lid) = 0;
                    recoloringSize(0)++;
                    break;
                  } else {
                    femv_colors(vid_d2) = 0;
                    recoloringSize(0)++;
                  }
                }
              }
            }
            if(found) break;
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
    
    
    boundary_verts = Kokkos::View<lno_t*, device_type>("boundary verts",n_local);
    typename Kokkos::View<lno_t*, device_type>::HostMirror boundary_verts_host = Kokkos::create_mirror_view(boundary_verts);

    for(size_t i = 0; i < n_local; i++){
      for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
        if((size_t)dist_adjs_host(j) >= n_local){
          boundary_verts_host(boundary_size++) = i;
          break;
        }
        bool found = false;
        for(offset_t k = dist_offsets_host(dist_adjs_host(j)); k < dist_offsets_host(dist_adjs_host(j)+1); k++){
          if((size_t)dist_adjs_host(k) >= n_local){
            boundary_verts_host(boundary_size++) = i;
            found = true;
            break;
          }
        }
        if(found) break;
      }
    }
    
    Kokkos::deep_copy(boundary_verts, boundary_verts_host);

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
  }

  public:
    AlgDistance2(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : AlgTwoGhostLayer<Adapter>(adapter_,pl_,env_,comm_){}

}; //end class


}//end namespace Zoltan2

#endif
