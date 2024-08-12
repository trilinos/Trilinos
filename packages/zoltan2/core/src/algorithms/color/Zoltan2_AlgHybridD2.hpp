// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_DISTANCE2_HPP_
#define _ZOLTAN2_DISTANCE2_HPP_

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
    using device_type = typename femv_t::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using host_exec = typename femv_t::host_view_type::device_type::execution_space;
    using host_mem = typename femv_t::host_view_type::device_type::memory_space;

  private:

    //This is both the serial and parallel local coloring.
    template <class ExecutionSpace, typename MemorySpace>
    void localColoring(const size_t nVtx,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> adjs_view,
		       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> offset_view,
		       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > vertex_list,
		       size_t vertex_list_size = 0,
		       bool use_vertex_based_coloring = false){
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
	  <offset_t, lno_t, lno_t, ExecutionSpace, MemorySpace, MemorySpace>;
      KernelHandle kh;

      //Instead of switching between vertex-based and net-based algorithms,
      //we only use the net-based algorithm, as it is faster than its
      //vertex-based counterpart.
      kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_NB_BIT);

      //vertex_list_size indicates whether we have provided a list of vertices to recolor
      //NB_BIT does not make use of this argument currently.
      if(vertex_list_size != 0){
        kh.get_distance2_graph_coloring_handle()->set_vertex_list(vertex_list, vertex_list_size);
      }

      //the verbose argument should carry through the local coloring
      kh.set_verbose(this->verbose);
      
      //set initial colors to be the colors from femv
      auto femvColors = femv->template getLocalView<Kokkos::Device<ExecutionSpace,MemorySpace> >(Tpetra::Access::ReadWrite);
      auto sv = subview(femvColors, Kokkos::ALL, 0);
      kh.get_distance2_graph_coloring_handle()->set_vertex_colors(sv);

      //call coloring
      KokkosGraph::Experimental::graph_color_distance2(&kh, nVtx, offset_view, adjs_view);
      
      
      //output total time
      if(this->verbose){
        std::cout<<"\nKokkosKernels Coloring: "
		 <<kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time()
		 <<"\n";
      }
    }

    //Entry point for device-based coloring
    virtual void colorInterior(const size_t nVtx, 
                       Kokkos::View<lno_t*,device_type> adjs_view,
                       Kokkos::View<offset_t*, device_type> offset_view,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, device_type> vertex_list,
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

    //Entry point for serial coloring
    virtual void colorInterior_serial(const size_t nVtx,
		       typename Kokkos::View<lno_t*, device_type>::HostMirror adjs_view,
		       typename Kokkos::View<offset_t*, device_type>::HostMirror offset_view,
		       Teuchos::RCP<femv_t> femv,
		       typename Kokkos::View<lno_t*, device_type>::HostMirror vertex_list,
		       size_t vertex_list_size = 0,
		       bool recolor=false){
      this->localColoring<host_exec, host_mem>(nVtx,
		                               adjs_view,
					       offset_view,
					       femv,
					       vertex_list,
					       vertex_list_size,
					       recolor);
    }
  public:
    //this function must be public due to Cuda Lambda restrictions.
    //It is both the serial and parallel conflict detection function.
    template< class ExecutionSpace, typename MemorySpace>
    void detectD2Conflicts(const size_t n_local,
		           Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace,MemorySpace>> dist_offsets,
			   Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> dist_adjs,
			   Kokkos::View<int*, Kokkos::Device<ExecutionSpace, MemorySpace>> femv_colors,
			   Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> boundary_verts_view,
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
                           Kokkos::View<size_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> recoloringSize,
			   Kokkos::View<int*, Kokkos::Device<ExecutionSpace, MemorySpace>> rand,
			   Kokkos::View<gno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> gid,
			   Kokkos::View<gno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>> ghost_degrees,
			   bool recolor_degrees){
      Kokkos::RangePolicy<ExecutionSpace> policy(0, boundary_verts_view.extent(0));
      size_t local_recoloring_size;
      Kokkos::parallel_reduce("D2 conflict detection",policy, KOKKOS_LAMBDA(const uint64_t& i,size_t& recoloring_size){
	//we only detect conflicts for vertices in the boundary
        const size_t curr_lid = boundary_verts_view(i);
        const int curr_color = femv_colors(curr_lid);
        const size_t vid_d1_adj_begin = dist_offsets(curr_lid);
        const size_t vid_d1_adj_end   = dist_offsets(curr_lid+1);
        const size_t curr_degree = vid_d1_adj_end - vid_d1_adj_begin;
        for(size_t vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++){
	  //check all distance-1 neighbors for conflicts
          size_t vid_d1 = dist_adjs(vid_d1_adj);
          size_t vid_d1_degree = 0;
	  //calculate the degree for degree-base recoloring
          if(vid_d1 < n_local){
            vid_d1_degree = dist_offsets(vid_d1+1) - dist_offsets(vid_d1);
          } else {
            vid_d1_degree = ghost_degrees(vid_d1-n_local);
          }
          if( vid_d1 != curr_lid && femv_colors(vid_d1) == curr_color){
            if(curr_degree < vid_d1_degree && recolor_degrees){
              femv_colors(curr_lid) = 0;
              recoloring_size++;
              break;//----------------------------------------------------                             
            } else if (vid_d1_degree < curr_degree && recolor_degrees){//|
              femv_colors(vid_d1) = 0;                                 //|
              recoloring_size++;                                       //|
            } else if(rand(curr_lid) < rand(vid_d1)){                  //|
              femv_colors(curr_lid) = 0;                               //|
              recoloring_size++;                                       //|
              break;//---------------------------------------------------|
            } else if(rand(vid_d1) < rand(curr_lid)){                  //|
              femv_colors(vid_d1) = 0;                                 //|
              recoloring_size++;                                       //|
            } else{                                                    //|
              if(gid(curr_lid) >= gid(vid_d1)){                        //|
                femv_colors(curr_lid) = 0;                             //|
                recoloring_size++;                                     //|
                break;//-------------------------------------------------|
              } else {                        //                         v
                femv_colors(vid_d1) = 0;      //If we uncolor the vertex whose
                recoloring_size++;            //neighbors we're checking, each         
              }                               //subsquent conflict check will
            }                                 //not do anything productive.
          }
          size_t d2_adj_begin = 0;
          size_t d2_adj_end   = 0;
          d2_adj_begin = dist_offsets(vid_d1);
          d2_adj_end   = dist_offsets(vid_d1+1);

	  //If we find a conflict that uncolors curr_lid, then we can safely stop
	  //detecting further conflicts. Since this is a nested loop, we need to 
	  //break twice, using the found boolean.
          bool found = false;
          for(size_t vid_d2_adj = d2_adj_begin; vid_d2_adj < d2_adj_end; vid_d2_adj++){
	    //check all distance-2 neighbors for conflicts
            const size_t vid_d2 = dist_adjs(vid_d2_adj);
            size_t vid_d2_degree = 0;
	    //calculate the degree for degree-based recoloring
            if(vid_d2 < n_local){
              vid_d2_degree = dist_offsets(vid_d2+1) - dist_offsets(vid_d2);
            } else {
              vid_d2_degree = ghost_degrees(vid_d2-n_local);
            }
            if(curr_lid != vid_d2 && femv_colors(vid_d2) == curr_color){
              if(curr_degree < vid_d2_degree && recolor_degrees){
                found = true;
                femv_colors(curr_lid) = 0;
                recoloring_size++;
                break;//---------------------------------------------------
              } else if(vid_d2_degree < curr_degree && recolor_degrees){//|
                femv_colors(vid_d2) = 0;                                //|
                recoloring_size++;                                      //|
              } else if(rand(curr_lid) < rand(vid_d2)){                 //|
                found = true;                                           //|
                femv_colors(curr_lid) = 0;                              //|
                recoloring_size++;                                      //|
                break;//--------------------------------------------------|
              } else if(rand(vid_d2) < rand(curr_lid)){                 //|
                femv_colors(vid_d2) = 0;                                //|
                recoloring_size++;                                      //|
              } else {                                                  //|
                if(gid(curr_lid) >= gid(vid_d2)){                       //|
                  found = true;                                         //|
                  femv_colors(curr_lid) = 0;                            //|
                  recoloring_size++;                                    //|
                  break;//------------------------------------------------|
                } else {                                                //|
                  femv_colors(vid_d2) = 0;                              //|
                  recoloring_size++;//                                    v
                }//              If we uncolor the vertex whose neighbors we're
              }  //              checking, each subsequent conflict check will 
            }    //              not do anything productive. We need this------
          }      //              to completely move on to the next vertex.    |
          if(found) break;//<--------------------------------------------------
        }
      },local_recoloring_size);
      Kokkos::deep_copy(recoloringSize,local_recoloring_size);
      Kokkos::fence();

      //update the verts_to_send and verts_to_recolor views.
      Kokkos::parallel_for("rebuild verts_to_send and verts_to_recolor",
		           Kokkos::RangePolicy<ExecutionSpace>(0,femv_colors.size()), 
			   KOKKOS_LAMBDA(const uint64_t& i){
        if(femv_colors(i) == 0){
	  //we only send vertices owned by the current process
          if(i < n_local){
            verts_to_send_view(verts_to_send_size_atomic(0)++) = i;
          }
	  //we need to recolor all vertices, for consistency.
          verts_to_recolor_view(verts_to_recolor_size_atomic(0)++) = i;
        }
      });
      Kokkos::fence();
      			   
    }
    
    //Entry point for parallel conflict detection
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

      this->detectD2Conflicts<execution_space, memory_space>(n_local,
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
  //Entry point for serial conflict detection
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
      this->detectD2Conflicts<host_exec, host_mem>(n_local,
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
    //count the number of boundary vertices to correctly allocate
    //the boundary vertex view on device.
    gno_t boundary_size_temp = 0; 
    for(size_t i = 0; i < n_local; i++){
      for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
        if((size_t)dist_adjs_host(j) >= n_local){
          boundary_size_temp++;
          break;
        }
        bool found = false;
        for(offset_t k = dist_offsets_host(dist_adjs_host(j)); k < dist_offsets_host(dist_adjs_host(j)+1); k++){
          if((size_t)dist_adjs_host(k) >= n_local){
            boundary_size_temp++;
            found = true;
            break;
          }
        }
        if(found) break;
      }
    }
        
    //create a host mirror to fill in the list of boundary vertices
    boundary_verts = Kokkos::View<lno_t*, device_type>("boundary verts",boundary_size_temp);
    typename Kokkos::View<lno_t*, device_type>::HostMirror boundary_verts_host = Kokkos::create_mirror_view(boundary_verts);
    
    //reset the boundary size count to use as an index to construct the view
    boundary_size_temp = 0;

    //a boundary vertex is any vertex within two edges of a ghost vertex.
    for(size_t i = 0; i < n_local; i++){
      for(offset_t j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
        if((size_t)dist_adjs_host(j) >= n_local){
          boundary_verts_host(boundary_size_temp++) = i;
          break;
        }
        bool found = false;
        for(offset_t k = dist_offsets_host(dist_adjs_host(j)); k < dist_offsets_host(dist_adjs_host(j)+1); k++){
          if((size_t)dist_adjs_host(k) >= n_local){
            boundary_verts_host(boundary_size_temp++) = i;
            found = true;
            break;
          }
        }
        if(found) break;
      }
    }
    //copy the boundary to the device view
    Kokkos::deep_copy(boundary_verts, boundary_verts_host);

    //initialize the list of verts to send
    Kokkos::parallel_for("init verts to send", 
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
    AlgDistance2(
      const RCP<const base_adapter_t> &adapter_,
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : AlgTwoGhostLayer<Adapter>(adapter_,pl_,env_,comm_){}

}; //end class


}//end namespace Zoltan2

#endif
