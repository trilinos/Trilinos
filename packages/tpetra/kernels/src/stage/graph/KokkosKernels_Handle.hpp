#include <KokkosKernels_GraphColorHandle.hpp>
#include <KokkosKernels_GaussSeidelHandle.hpp>
#include <KokkosKernels_SPGEMMHandle.hpp>
#include "utils/KokkosKernels_ExecSpaceUtils.hpp"
#ifndef _KOKKOSKERNELHANDLE_HPP
#define _KOKKOSKERNELHANDLE_HPP

namespace KokkosKernels{

namespace Experimental{



template <class lno_row_view_t_, class lno_nnz_view_t_, class scalar_nnz_view_t_,
          class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class KokkosKernelsHandle{
public:

  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename in_lno_row_view_t::non_const_value_type size_type;
  typedef typename in_lno_row_view_t::array_layout row_lno_view_array_layout;
  typedef typename in_lno_row_view_t::device_type row_lno_view_device_t;
  typedef typename in_lno_row_view_t::memory_traits row_lno_view_memory_traits;
  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type
  typedef typename in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename in_lno_nnz_view_t::array_layout nnz_lno_view_array_layout;
  typedef typename in_lno_nnz_view_t::device_type nnz_lno_view_device_t;
  typedef typename in_lno_nnz_view_t::memory_traits nnz_lno_view_memory_traits;
  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type
  typedef typename in_scalar_nnz_view_t::non_const_value_type nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::array_layout nnz_scalar_view_array_layout;
  typedef typename in_scalar_nnz_view_t::device_type nnz_scalar_view_device_t;
  typedef typename in_scalar_nnz_view_t::memory_traits nnz_scalar_view_memory_traits;
  typedef typename in_scalar_nnz_view_t::HostMirror nnz_scalar_view_t; //Host view type
  typedef typename in_lno_row_view_t::const_value_type const_row_lno_t;
  typedef typename in_lno_row_view_t::const_value_type const_size_type;
  typedef typename in_lno_row_view_t::non_const_value_type non_const_row_lno_t;
  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;
  typedef typename in_lno_nnz_view_t::const_value_type const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;
  typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;
  typedef typename in_scalar_nnz_view_t::const_data_type const_nnz_scalar_t; //nnz_scalar_t
  typedef typename in_scalar_nnz_view_t::non_const_data_type non_const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::const_type const_scalar_nnz_view_t;
  typedef typename in_scalar_nnz_view_t::non_const_type non_const_scalar_nnz_view_t;
  typedef typename Graph::GraphColoringHandle
      <in_lno_row_view_t, non_const_lno_nnz_view_t, in_lno_nnz_view_t,
      ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> GraphColoringHandleType;
  typedef typename Graph::GaussSeidelHandle
      <in_lno_row_view_t, in_lno_nnz_view_t, in_scalar_nnz_view_t,
      ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> GaussSeidelHandleType;

  typedef typename Graph::SPGEMMHandle
      <in_lno_row_view_t, in_lno_nnz_view_t, in_scalar_nnz_view_t,
      ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> SPGEMMHandleType;

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

private:
  GraphColoringHandleType *gcHandle;
  GaussSeidelHandleType *gsHandle;
  SPGEMMHandleType *spgemmHandle;
  int team_work_size;
  size_t shared_memory_size;
  int suggested_team_size;
  KokkosKernels::Experimental::Util::ExecSpaceType my_exec_space;
  bool use_dynamic_scheduling;
public:



  KokkosKernelsHandle():
      gcHandle(NULL), gsHandle(NULL),spgemmHandle(NULL),
      team_work_size (-1), shared_memory_size(16128),
      suggested_team_size(-1),
      my_exec_space(KokkosKernels::Experimental::Util::kk_get_exec_space_type<HandleExecSpace>()),
      use_dynamic_scheduling(false){}

  ~KokkosKernelsHandle(){
    this->destroy_gs_handle();
    this->destroy_graph_coloring_handle();
    this->destroy_spgemm_handle();
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

  /**
   * \brief Returns the enum type for the execution space.
   */
  KokkosKernels::Experimental::Util::ExecSpaceType get_handle_exec_space(){
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
      return team_size;
      /*
      if (my_exec_space == KokkosKernels::Experimental::Util::Exec_CUDA){
        return team_size;
      }
      else {
        return overall_work_size / (concurrency * team_size) ;
      }
      */
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
    return KokkosKernels::Experimental::Util::kk_get_suggested_vector_size(nr, nnz, my_exec_space);
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



  /**
   * \brief Returns the team size, either set by the user or suggested by the handle.
   * \param vector_size: suggested vector size by the handle.
   */
  int get_suggested_team_size(const int vector_size){
    if (this->suggested_team_size != -1){
      return this->suggested_team_size;
    }
    else {
      return KokkosKernels::Experimental::Util::kk_get_suggested_team_size(vector_size, my_exec_space);
    }
  }




  SPGEMMHandleType *get_spgemm_handle(){
    return this->spgemmHandle;
  }

  void create_spgemm_handle(Graph::SPGEMMAlgorithm spgemm_algo = Graph::SPGEMM_DEFAULT){
    this->destroy_spgemm_handle();
    this->spgemmHandle = new SPGEMMHandleType(spgemm_algo);

  }
  void destroy_spgemm_handle(){
    if (this->spgemmHandle != NULL){
      delete this->spgemmHandle;
      this->spgemmHandle = NULL;
    }
  }

  GraphColoringHandleType *get_graph_coloring_handle(){
    return this->gcHandle;
  }
  void create_graph_coloring_handle(Graph::ColoringAlgorithm coloring_type = Graph::COLORING_DEFAULT){
    this->destroy_graph_coloring_handle();
    this->gcHandle = new GraphColoringHandleType();
    this->gcHandle->set_algorithm(coloring_type, true);
  }
  void destroy_graph_coloring_handle(){
    if (this->gcHandle != NULL){
      delete this->gcHandle;
      this->gcHandle = NULL;
    }
  }


  GaussSeidelHandleType *get_gs_handle(){
    return this->gsHandle;
  }
  void create_gs_handle(
    Graph::GSAlgorithm gs_algorithm = Graph::GS_DEFAULT){
    this->destroy_gs_handle();
    this->gsHandle = new GaussSeidelHandleType(gs_algorithm);
  }
  void destroy_gs_handle(){
    if (this->gsHandle != NULL){
      if (this->gsHandle->is_owner_of_coloring()){
        this->destroy_graph_coloring_handle();
      }
      delete this->gsHandle;
      this->gsHandle = NULL;
    }
  }



};

}
}

#endif //_KOKKOSKERNELHANDLE_HPP
