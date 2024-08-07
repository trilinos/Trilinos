// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_
#define _ZOLTAN2_DISTANCE1_2GHOSTLAYER_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

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
    using device_type = typename femv_t::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using host_exec = typename femv_t::host_view_type::device_type::execution_space;
    using host_mem = typename femv_t::host_view_type::device_type::memory_space;

  private:
   
    template <class ExecutionSpace, typename MemorySpace>
    void localColoring(const size_t nVtx,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> adjs_view,
		       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> offset_view,
		       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > vertex_list,
		       size_t vertex_list_size = 0,
		       bool use_vertex_based_coloring=false){
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
	  <offset_t, lno_t, lno_t, ExecutionSpace, MemorySpace, MemorySpace>;
      using lno_row_view_t = Kokkos::View<offset_t*,Kokkos::Device<ExecutionSpace, MemorySpace>>;
      using lno_nnz_view_t = Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;
      KernelHandle kh;

      //pick which KokkosKernels algorithm to use.
      //this boolean's value is based on max degree.
      if(use_vertex_based_coloring){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
      }

      //vertex_list_size indicates whether we have provided a list of vertices to recolor
      //only makes a difference if the algorithm to be used is VBBIT
      if(vertex_list_size != 0){
        kh.get_graph_coloring_handle()->set_vertex_list(vertex_list, vertex_list_size);
      }
      
      //the verbose argument should carry through the local coloring
      kh.set_verbose(this->verbose);

      //set initial colors to be the colors from the femv
      auto femvColors = femv->template getLocalView<Kokkos::Device<ExecutionSpace,MemorySpace> >(Tpetra::Access::ReadWrite);
      auto sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      
      //if running verbosely, also report local color timing breakdown
      kh.get_graph_coloring_handle()->set_tictoc(this->verbose);
      
      //call coloring
      KokkosGraph::Experimental::graph_color_symbolic
	      <KernelHandle, lno_row_view_t, lno_nnz_view_t> (&kh, 
			                                      nVtx, 
							      nVtx, 
							      offset_view, 
							      adjs_view);

      
      //output total time and #iterations
      if(this->verbose){
        std::cout<<"\nKokkosKernels Coloring: "
		 <<kh.get_graph_coloring_handle()->get_overall_coloring_time()
		 <<" iterations: "
		 <<kh.get_graph_coloring_handle()->get_num_phases()
		 <<"\n\n";
      }
    }

    virtual void colorInterior(const size_t nVtx,
                       Kokkos::View<lno_t*,device_type > adjs_view,
                       Kokkos::View<offset_t*,device_type > offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, device_type > vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){

      this->localColoring<execution_space, memory_space>(nVtx,
		                                         adjs_view,
							 offset_view,
							 femv,
							 vertex_list,
							 vertex_list_size,
							 recolor);
    }
   
    virtual void colorInterior_serial(const size_t nVtx,
                       typename Kokkos::View<lno_t*, device_type >::HostMirror adjs_view,
                       typename Kokkos::View<offset_t*,device_type >::HostMirror offset_view,
                       Teuchos::RCP<femv_t> femv,
                       typename Kokkos::View<lno_t*, device_type>::HostMirror vertex_list,
                       size_t vertex_list_size = 0,
                       bool recolor=false) {
	this->localColoring<host_exec, host_mem>(nVtx,
		                                 adjs_view,
				  	         offset_view,
						 femv,
						 vertex_list,
						 vertex_list_size,
						 recolor);	    
    }

  public:
    template <class ExecutionSpace, typename MemorySpace>
    void detectD1Conflicts(const size_t n_local,
	                   Kokkos::View<offset_t*, 
                                        Kokkos::Device<ExecutionSpace, MemorySpace>> dist_offsets,
                           Kokkos::View<lno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace>> dist_adjs,
			   Kokkos::View<int*,
			                Kokkos::Device<ExecutionSpace, MemorySpace>> femv_colors, 
			   Kokkos::View<lno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace>> boundary_verts_view,
			   Kokkos::View<lno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace> > verts_to_recolor_view,
			   Kokkos::View<int*,
			                Kokkos::Device<ExecutionSpace, MemorySpace>,
					Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_recolor_size_atomic,
			   Kokkos::View<lno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace> > verts_to_send_view,
			   Kokkos::View<size_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace>,
					Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_size_atomic,
			   Kokkos::View<size_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > recoloringSize,
			   Kokkos::View<int*,
			                Kokkos::Device<ExecutionSpace, MemorySpace> > rand,
			   Kokkos::View<gno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace> > gid,
			   Kokkos::View<gno_t*,
			                Kokkos::Device<ExecutionSpace, MemorySpace> > ghost_degrees,
			   bool recolor_degrees){
      size_t local_recoloring_size;
      Kokkos::RangePolicy<ExecutionSpace> policy(n_local,rand.size());
      Kokkos::parallel_reduce("D1-2GL Conflict Detection",policy, KOKKOS_LAMBDA (const int& i, size_t& recoloring_size){
        lno_t localIdx = i;
        int currColor = femv_colors(localIdx);
        int currDegree = ghost_degrees(i-n_local);
        for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
          int nborColor = femv_colors(dist_adjs(j));
          int nborDegree = 0;
          if((size_t)dist_adjs(j) < n_local) nborDegree = dist_offsets(dist_adjs(j)+1) - dist_offsets(dist_adjs(j));
          else nborDegree = ghost_degrees(dist_adjs(j) - n_local);
          if(currColor == nborColor ){
            if(currDegree < nborDegree && recolor_degrees){
              femv_colors(localIdx) = 0;
              recoloring_size++;
              break;
            }else if(nborDegree < currDegree && recolor_degrees){
              femv_colors(dist_adjs(j)) = 0;
              recoloring_size++;
            }else if(rand(localIdx) > rand(dist_adjs(j))){
              recoloring_size++;
              femv_colors(localIdx) = 0;
              break;
            }else if(rand(dist_adjs(j)) > rand(localIdx)){
                recoloring_size++;
                femv_colors(dist_adjs(j)) = 0;
            } else {
              if (gid(localIdx) >= gid(dist_adjs(j))){
                femv_colors(localIdx) = 0;
                recoloring_size++;
                break;
              } else {
                  femv_colors(dist_adjs(j)) = 0;
                  recoloring_size++;
              }
            }
          }
        }
      },local_recoloring_size);
      Kokkos::deep_copy(recoloringSize, local_recoloring_size);
      Kokkos::fence();
      Kokkos::parallel_for("rebuild verts_to_send and verts_to_recolor",
		           Kokkos::RangePolicy<ExecutionSpace>(0,femv_colors.size()), 
			   KOKKOS_LAMBDA (const size_t& i){
        if(femv_colors(i) == 0){
          if(i < n_local){
            verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
          }
          verts_to_recolor_view(verts_to_recolor_size_atomic(0)++) = i;
        }
      });
      Kokkos::fence();
      			   
    }
			   
    virtual void detectConflicts(const size_t n_local,
                                 Kokkos::View<offset_t*, device_type > dist_offsets_dev,
                                 Kokkos::View<lno_t*, device_type > dist_adjs_dev,
                                 Kokkos::View<int*,device_type > femv_colors, 
                                 Kokkos::View<lno_t*, device_type > boundary_verts_view,
                                 Kokkos::View<lno_t*,
                                              device_type > verts_to_recolor_view,
                                 Kokkos::View<int*,
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_recolor_size_atomic,
                                 Kokkos::View<lno_t*,
                                              device_type > verts_to_send_view,
                                 Kokkos::View<size_t*,
                                              device_type,
                                              Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic,
				 Kokkos::View<size_t*, device_type> recoloringSize,
                                 Kokkos::View<int*,
                                              device_type> rand,
                                 Kokkos::View<gno_t*,
                                              device_type> gid,
                                 Kokkos::View<gno_t*,
                                              device_type> ghost_degrees,
				 bool recolor_degrees){
      this->detectD1Conflicts<execution_space, memory_space>(n_local,
		                                             dist_offsets_dev,
							     dist_adjs_dev,
							     femv_colors,
							     boundary_verts_view,
							     verts_to_recolor_view,
							     verts_to_recolor_size_atomic,
							     verts_to_send_view,
							     verts_to_send_size_atomic,
							     recoloringSize,
							     rand,
							     gid,
							     ghost_degrees,
							     recolor_degrees);
    }

    virtual void detectConflicts_serial(const size_t n_local,
                                 typename Kokkos::View<offset_t*, device_type >::HostMirror dist_offsets_host,
                                 typename Kokkos::View<lno_t*, device_type >::HostMirror dist_adjs_host,
                                 typename Kokkos::View<int*,device_type >::HostMirror femv_colors,
                                 typename Kokkos::View<lno_t*, device_type >::HostMirror boundary_verts_view,
                                 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_recolor,
                                 typename Kokkos::View<int*,device_type>::HostMirror verts_to_recolor_size,
                                 typename Kokkos::View<lno_t*,device_type>::HostMirror verts_to_send,
                                 typename Kokkos::View<size_t*,device_type>::HostMirror verts_to_send_size,
				 typename Kokkos::View<size_t*, device_type>::HostMirror recoloringSize,
                                 typename Kokkos::View<int*,  device_type>::HostMirror rand,
                                 typename Kokkos::View<gno_t*,device_type>::HostMirror gid,
                                 typename Kokkos::View<gno_t*,device_type>::HostMirror ghost_degrees,
				 bool recolor_degrees) {
      this->detectD1Conflicts<host_exec, host_mem >(n_local,
		                                    dist_offsets_host,
						    dist_adjs_host,
						    femv_colors,
						    boundary_verts_view,
						    verts_to_recolor,
						    verts_to_recolor_size,
						    verts_to_send,
						    verts_to_send_size,
						    recoloringSize,
						    rand,
						    gid,
						    ghost_degrees,
						    recolor_degrees);

    }

    virtual void constructBoundary(const size_t n_local,
                                   Kokkos::View<offset_t*, device_type> dist_offsets_dev,
                                   Kokkos::View<lno_t*, device_type> dist_adjs_dev,
				   typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host,
				   typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host,
                                   Kokkos::View<lno_t*, device_type>& boundary_verts,
                                   Kokkos::View<lno_t*,
                                                device_type > verts_to_send_view,
                                   Kokkos::View<size_t*,
                                                device_type,
                                                Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_size_atomic){

      Kokkos::parallel_for("constructBoundary",
        Kokkos::RangePolicy<execution_space, int>(0,n_local), 
        KOKKOS_LAMBDA(const int& i){
        for(offset_t j = dist_offsets_dev(i); j < dist_offsets_dev(i+1); j++){
          if((size_t)dist_adjs_dev(j) >= n_local){
            verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
            break;
          }
          bool found = false;
          for(offset_t k = dist_offsets_dev(dist_adjs_dev(j)); k < dist_offsets_dev(dist_adjs_dev(j)+1); k++){
            if((size_t)dist_adjs_dev(k) >= n_local){
              verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
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
