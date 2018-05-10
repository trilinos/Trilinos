/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <KokkosKernels_ExecSpaceUtils.hpp>
#include "KokkosGraph_GraphColorHandle.hpp"
#include "KokkosSparse_gauss_seidel_handle.hpp"
#include "KokkosSparse_spgemm_handle.hpp"
#include "KokkosSparse_spadd_handle.hpp"
#ifndef _KOKKOSKERNELHANDLE_HPP
#define _KOKKOSKERNELHANDLE_HPP

namespace KokkosKernels{

namespace Experimental{

template <class size_type_, class lno_t_, class scalar_t_,
          class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class KokkosKernelsHandle{
public:

  typedef typename ExecutionSpace::execution_space HandleExecSpace;
  typedef typename TemporaryMemorySpace::memory_space HandleTempMemorySpace;
  typedef typename PersistentMemorySpace::memory_space HandlePersistentMemorySpace;
  //typedef Kokkos::Device<ExecutionSpace::execution_space,ExecutionSpace::memory_space>  HandleExecSpace;
  //typedef Kokkos::Device<TemporaryMemorySpace::execution_space,TemporaryMemorySpace::memory_space> HandleTempMemorySpace;
  //typedef Kokkos::Device<PersistentMemorySpace::execution_space,PersistentMemorySpace::memory_space> HandlePersistentMemorySpace;


  typedef typename std::remove_const<size_type_>::type  size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type  nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type  nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  template <  typename right_size_type_, typename right_lno_t_, typename right_scalar_t_,
  	  	  	  typename right_ExecutionSpace, typename right_TemporaryMemorySpace, typename right_PersistentMemorySpace>
  //KokkosKernelsHandle<const_size_type,const_nnz_lno_t, const_nnz_scalar_t, HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace> &operator=
  KokkosKernelsHandle
  	  (KokkosKernelsHandle<right_size_type_, right_lno_t_, right_scalar_t_,
	  	  	  right_ExecutionSpace, right_TemporaryMemorySpace, right_PersistentMemorySpace> & right_side_handle){

	  static_assert (std::is_same<size_type_,  const_size_type>::value,
			  "Kernel handle left hand side should have const size type in assignment");
	  static_assert (std::is_same<lno_t_,  const_nnz_lno_t>::value,
			  "Kernel handle left hand side should have const lno type in assignment");
	  static_assert (std::is_same<scalar_t_,  const_nnz_scalar_t>::value,
			  "Kernel handle left hand side should have const scalar type in assignment");

	  static_assert (std::is_same<ExecutionSpace,  HandleExecSpace>::value,
			  "Kernel handle left hand side should have execution space in assignment");
	  static_assert (std::is_same<TemporaryMemorySpace,  HandleTempMemorySpace>::value,
			  "Kernel handle left hand side should have temp memory space in assignment");
	  static_assert (std::is_same<PersistentMemorySpace,  HandlePersistentMemorySpace>::value,
			  "Kernel handle left hand side should have persistent memory space in assignment");

	  typedef typename  std::remove_const<right_size_type_>::type  nonconst_right_size_type;
	  typedef const nonconst_right_size_type const_right_size_type;

	  typedef typename  std::remove_const<right_lno_t_>::type  nonconst_right_lno_t;
	  typedef const nonconst_right_lno_t const_right_lno_t;

	  typedef typename  std::remove_const<right_scalar_t_>::type  nonconst_right_scalar_t;
	  typedef const nonconst_right_scalar_t const_right_scalar_t;




	  static_assert (std::is_same<size_type_,  const_right_size_type>::value,
			  "Kernel handle left and right sides should have same size type in assignment");
	  static_assert (std::is_same<lno_t_,  const_right_lno_t>::value,
			  "Kernel handle left and right sides should have same lno type in assignment");
	  static_assert (std::is_same<scalar_t_,  const_right_scalar_t>::value,
			  "Kernel handle left and right sides should have same scalar type in assignment");



	  static_assert (std::is_same<typename ExecutionSpace::execution_space,  typename right_ExecutionSpace::execution_space>::value,
			  "Kernel handle left and right sides should have same execution_space in assignment");
	  /*
	  static_assert (std::is_same<typename TemporaryMemorySpace::execution_space, typename right_TemporaryMemorySpace::execution_space>::value,
			  "Kernel handle left and right sides should have same TemporaryMemorySpace in assignment");
	  static_assert (std::is_same<typename PersistentMemorySpace::execution_space, typename right_PersistentMemorySpace::execution_space>::value,
			  "Kernel handle left and right sides should have same PersistentMemorySpace in assignment");
      */
	  static_assert (std::is_same<typename ExecutionSpace::memory_space, typename right_ExecutionSpace::memory_space>::value,
			  "Kernel handle left and right sides should have same ExecutionSpace in assignment");
	  static_assert (std::is_same<typename TemporaryMemorySpace::memory_space, typename right_TemporaryMemorySpace::memory_space>::value,
			  "Kernel handle left and right sides should have same TemporaryMemorySpace in assignment");
	  static_assert (std::is_same<typename PersistentMemorySpace::memory_space, typename right_PersistentMemorySpace::memory_space>::value,
			  "Kernel handle left and right sides should have same PersistentMemorySpace in assignment");


	  this->gcHandle = right_side_handle.get_graph_coloring_handle();
	  this->gsHandle = right_side_handle.get_gs_handle();
	  this->spgemmHandle = right_side_handle.get_spgemm_handle();


	  this->team_work_size = right_side_handle.get_set_team_work_size();
	  this->shared_memory_size = right_side_handle.get_shmem_size();
	  this->suggested_team_size = right_side_handle.get_set_suggested_team_size();

	  this->my_exec_space = right_side_handle.get_handle_exec_space();
	  this->use_dynamic_scheduling = right_side_handle.is_dynamic_scheduling();
	  this->KKVERBOSE = right_side_handle.get_verbose();
	  this->vector_size = right_side_handle.get_set_suggested_vector_size();

	  is_owner_of_the_gc_handle = false;
	  is_owner_of_the_gs_handle = false;
	  is_owner_of_the_spgemm_handle = false;
	  is_owner_of_the_spadd_handle = false;
	  //return *this;
  }


  typedef typename KokkosGraph::GraphColoringHandle
      <const_size_type, const_nnz_lno_t, const_nnz_lno_t,
	  HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace> GraphColoringHandleType;
  typedef typename KokkosSparse::GaussSeidelHandle
      <const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
	  HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace> GaussSeidelHandleType;

  typedef typename KokkosSparse::SPGEMMHandle
      <const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
	  HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace> SPGEMMHandleType;

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> in_scalar_nnz_view_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> size_type_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> size_type_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type
  typedef typename size_type_persistent_work_view_t::HostMirror size_type_persistent_work_host_view_t; //Host view type
  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type
  typedef typename Kokkos::View<bool *, HandlePersistentMemorySpace> bool_persistent_view_t;
  typedef typename Kokkos::View<bool *, HandleTempMemorySpace> bool_temp_view_t;

  typedef typename KokkosSparse::SPADDHandle
         <row_lno_temp_work_view_t,
          nnz_lno_temp_work_view_t,
          scalar_temp_work_view_t,
          HandleExecSpace,
          HandleTempMemorySpace> SPADDHandleType;

private:

  GraphColoringHandleType *gcHandle;
  GaussSeidelHandleType *gsHandle;
  SPGEMMHandleType *spgemmHandle;
  SPADDHandleType *spaddHandle;

  int team_work_size;
  size_t shared_memory_size;
  int suggested_team_size;

  KokkosKernels::Impl::ExecSpaceType my_exec_space;
  bool use_dynamic_scheduling;
  bool KKVERBOSE;
  int vector_size;

  bool is_owner_of_the_gc_handle;
  bool is_owner_of_the_gs_handle;
  bool is_owner_of_the_spgemm_handle;
  bool is_owner_of_the_spadd_handle;


public:



  KokkosKernelsHandle():
      gcHandle(NULL), gsHandle(NULL),spgemmHandle(NULL),spaddHandle(NULL),
      team_work_size (-1), shared_memory_size(16128),
      suggested_team_size(-1),
      my_exec_space(KokkosKernels::Impl::kk_get_exec_space_type<HandleExecSpace>()),
      use_dynamic_scheduling(true), KKVERBOSE(false),vector_size(-1),
	  is_owner_of_the_gc_handle(true), is_owner_of_the_gs_handle(true), is_owner_of_the_spgemm_handle(true),
    is_owner_of_the_spadd_handle(true) {}

  ~KokkosKernelsHandle(){
    this->destroy_gs_handle();
    this->destroy_graph_coloring_handle();
    this->destroy_spgemm_handle();
    this->destroy_spadd_handle();
  }


  void set_verbose(bool verbose_){
    this->KKVERBOSE = verbose_;
  }
  bool get_verbose (){
    return this->KKVERBOSE;
  }
  /**
   * \brief Sets the chunk size to be handled by each team.
   * Big chunks are good, they reduce the overhead that occurs at the every iteration.
   * However, too big chunks causes load-imbalances.
   * If this is not set, the algorithm will use teamsize as chunk sizes.
   * On GPUs, this is usually okay. On CPUs, we get better performance for bigger chunksizes.
   * This depends on how big is your work. On SPGEMM a single work is the multiplication for a row.
   * For laplace we had 50-500 FLOPs per row for AP and RxAP. Best was to use 256 chunks. ~12800 - 128K flops per team was fine.
   * For brick we had 250 - 2500 FLOPs per  row for AP and RxAP. Best was to use 256 chunks. 64K to 640K flops per team was fine.
   * But on EMPIRE 1200 - 20800 FLOPs per row for AP and RxAP. Best was to use 16 chunks. 19K to 330K was flops per team was fine.
   * This bases on the load balancing issues introduced. A general way is to have 100K flops per team for SPGEMM.
   * \param team_work_size_: input, the size of the chunks.
   */
  void set_team_work_size(const int team_work_size_){
    this->team_work_size = team_work_size_;
  }
  int get_set_team_work_size(){
    return this->team_work_size;
  }
  /**
   * \brief Returns the enum type for the execution space.
   */
  KokkosKernels::Impl::ExecSpaceType get_handle_exec_space(){
    return this->my_exec_space;
  }

  /**
   * \brief Returns the suggested team work size. If set with set_team_work_size,
   * it will return the set value. Otherwise it will return the teamsize.
   * \param team_size: input, team size used by the kernel.
   * \param concurrency: input, the number of threads overall. Not used currently.
   * \param overall_work_size: The overall work size.
   */
  int get_team_work_size(const int team_size, const int concurrency, const nnz_lno_t overall_work_size){
    if (this->team_work_size != -1){
      return this->team_work_size;
    }
    else {
      if (my_exec_space == KokkosKernels::Impl::Exec_CUDA){
        return team_size;
      }
      else {
        return 16;
      }

    }
  }

  /**
   * \brief Sets whether to use dynamic scheduling or
   * not for those kernels where load-imbalances might occur.
   * \input is_dynamic: true or false -> dynamic or static scheduling.
   */
  void set_dynamic_scheduling(const bool is_dynamic){
    this->use_dynamic_scheduling = is_dynamic;
  }

  /**
   * \brief Returns true or false, use dynamic scheduling or not.
   */
  bool is_dynamic_scheduling(){
    return this->use_dynamic_scheduling;
  }



  /**
   * \brief sets the shared memory size to be used by the kernels using shared memory on GPUs.
   * \param shared_memory_size: input, shared memory size to be used by the kernel.   *
   */
  void set_shmem_size(const size_t shared_memory_size_){
    this->shared_memory_size = shared_memory_size_;
  }

  /**
   * \brief Returns the shared memory size suggested by the handle.
   */
  size_t get_shmem_size(){
    return shared_memory_size;
  }


  /**
   * \brief Returns the suggested vector size based on the execution space.
   * Basically, it calculates the average degree, and rounds it up to the
   * closes power of 2 on GPUs. on CPUs, it returns 1.
   * \param nr: number of rows, of vertices.
   * \param nnz: number of nonzeroes, or edges.
   */
  int get_suggested_vector_size(const size_t nr,const size_t nnz){
    if (vector_size == -1){
      return KokkosKernels::Impl::kk_get_suggested_vector_size(nr, nnz, my_exec_space);
    }
    else {
      return vector_size;
    }
  }

  void set_suggested_vector_size(int vector_size_){
    this->vector_size = vector_size_;
  }

  int get_set_suggested_vector_size(){
    return this->vector_size;
  }
  /**
   * \brief Sets the team size to be used by the kernels. On GPUs and CPUs
   * usually the defaults are fine. But on CPUs with hyperthreads it might be
   * worth of trying different team sizes.
   * \param suggested_team_size_: team size to set.
   */
  void set_suggested_team_size(const int suggested_team_size_){
    this->suggested_team_size = suggested_team_size_;
  }

  int get_set_suggested_team_size(){
    return this->suggested_team_size;
  }


  /**
   * \brief Returns the team size, either set by the user or suggested by the handle.
   * \param vector_size: suggested vector size by the handle.
   */
  int get_suggested_team_size(const int vector_size_){
    if (this->suggested_team_size != -1){
      return this->suggested_team_size;
    }
    else {
      return KokkosKernels::Impl::kk_get_suggested_team_size(vector_size_, my_exec_space);
    }
  }




  SPGEMMHandleType *get_spgemm_handle(){
    return this->spgemmHandle;
  }

  void create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm spgemm_algo = KokkosSparse::SPGEMM_DEFAULT){
    this->destroy_spgemm_handle();
    this->is_owner_of_the_spgemm_handle = true;
    this->spgemmHandle = new SPGEMMHandleType(spgemm_algo);
  }

  void destroy_spgemm_handle(){
    if (is_owner_of_the_spgemm_handle && this->spgemmHandle != NULL){
      delete this->spgemmHandle;
      this->spgemmHandle = NULL;
    }
  }

  GraphColoringHandleType *get_graph_coloring_handle(){
    return this->gcHandle;
  }
  void create_graph_coloring_handle(KokkosGraph::ColoringAlgorithm coloring_type = KokkosGraph::COLORING_DEFAULT){
    this->destroy_graph_coloring_handle();
    this->is_owner_of_the_gc_handle = true;
    this->gcHandle = new GraphColoringHandleType();
    this->gcHandle->set_algorithm(coloring_type, true);
    this->gcHandle->set_tictoc(KKVERBOSE);

  }
  void destroy_graph_coloring_handle(){
    if (is_owner_of_the_gc_handle &&  this->gcHandle != NULL){
      delete this->gcHandle;
      this->gcHandle = NULL;
    }
  }


  GaussSeidelHandleType *get_gs_handle(){
    return this->gsHandle;
  }
  void create_gs_handle(
    KokkosSparse::GSAlgorithm gs_algorithm = KokkosSparse::GS_DEFAULT){
    this->destroy_gs_handle();
    this->is_owner_of_the_gs_handle = true;
    this->gsHandle = new GaussSeidelHandleType(gs_algorithm);
  }
  void destroy_gs_handle(){
    if (is_owner_of_the_gs_handle && this->gsHandle != NULL){
      if (this->gsHandle->is_owner_of_coloring()){
        this->destroy_graph_coloring_handle();
      }
      delete this->gsHandle;
      this->gsHandle = NULL;
    }
  }


  SPADDHandleType *get_spadd_handle(){
    return this->spaddHandle;
  }

  void create_spadd_handle(bool input_sorted) {
    this->destroy_spadd_handle();
    this->is_owner_of_the_spadd_handle = true;
    this->spaddHandle = new SPADDHandleType(input_sorted);
  }

  void destroy_spadd_handle(){
    if (is_owner_of_the_spadd_handle && this->spaddHandle != NULL)
    {
      delete this->spaddHandle;
      this->spaddHandle = NULL;
    }
  }

};

}
}

#endif //_KOKKOSKERNELHANDLE_HPP
