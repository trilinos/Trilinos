#include <GraphColoringHandle.hpp>
#include <GaussSeidelHandle.hpp>
#include <SPGEMMHandle.hpp>
#ifndef _KOKKOSKERNELHANDLE_HPP
#define _KOKKOSKERNELHANDLE_HPP

namespace KokkosKernels{

namespace Experimental{



template <class idx_array_type_, class idx_edge_array_type_, class value_array_type_,
          class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class KokkosKernelsHandle{
public:

  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;



  typedef idx_array_type_ idx_array_type;
  typedef idx_edge_array_type_ idx_edge_array_type;
  typedef value_array_type_ value_array_type;

  typedef typename idx_array_type::value_type idx;
  typedef typename idx_array_type::array_layout idx_array_layout;
  typedef typename idx_array_type::device_type idx_device_type;
  typedef typename idx_array_type::memory_traits idx_memory_traits;
  typedef typename idx_array_type::HostMirror host_view_type; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename idx_edge_array_type::value_type idx_edge;
  typedef typename idx_edge_array_type::array_layout idx_edge_array_layout;
  typedef typename idx_edge_array_type::device_type idx_edge_device_type;
  typedef typename idx_edge_array_type::memory_traits idx_edge_memory_traits;
  typedef typename idx_edge_array_type::HostMirror host_edge_view_type; //Host view type
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename value_array_type::value_type value_type;
  typedef typename value_array_type::array_layout value_type_array_layout;
  typedef typename value_array_type::device_type value_type_device_type;
  typedef typename value_array_type::memory_traits value_type_memory_traits;
  typedef typename value_array_type::HostMirror host_value_view_type; //Host view type

  typedef typename Graph::GraphColoringHandle
      <idx_array_type, idx_array_type, idx_edge_array_type, ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> GraphColoringHandleType;

  typedef typename Graph::GaussSeidelHandle
      <idx_array_type, idx_edge_array_type, value_array_type, ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> GaussSeidelHandleType;

  typedef typename Graph::SPGEMMHandle
      <idx_array_type, idx_edge_array_type, value_array_type, ExecutionSpace, TemporaryMemorySpace, PersistentMemorySpace> SPGEMMHandleType;


private:
  GraphColoringHandleType *gcHandle;
  GaussSeidelHandleType *gsHandle;
  SPGEMMHandleType *spgemmHandle;
  //idx_array_type row_map;
  //idx_edge_array_type entries;
  //value_array_type values;

public:



  KokkosKernelsHandle():gcHandle(NULL), gsHandle(NULL),spgemmHandle(NULL){}
  ~KokkosKernelsHandle(){

    this->destroy_gs_handle();
    this->destroy_graph_coloring_handle();
    this->destroy_spgemm_handle();
  }

  SPGEMMHandleType *get_spgemm_handle(){
    return this->spgemmHandle;
  }

  void create_spgemm_handle(Graph::SPGEMMAlgorithm spgemm_algo = Graph::SPGEMM_DEFAULT){
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


  //idx_array_type get_row_map(){return this->row_map;}
  //idx get_num_rows(){return this->row_map.dimension_0() - 1;}
  //idx get_num_nonzeroes(){return this->entries.dimension_0();}

  //idx_edge_array_type get_entries(){return this->entries;}
  //value_array_type get_values(){return this->values;}
  //void set_row_map(const idx_array_type &rm){this->row_map = rm;}
  //void set_entries(const idx_edge_array_type &e){this->entries = e;}
  //void set_values(const value_array_type &v){this->values = v;}



};

}
}

#endif //_KOKKOSKERNELHANDLE_HPP
