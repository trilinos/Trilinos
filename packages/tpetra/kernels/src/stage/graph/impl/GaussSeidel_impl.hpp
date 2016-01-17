#include "GraphColor.hpp"
#include "KokkosKernelsUtils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_MemoryTraits.hpp>

#ifndef _KOKKOSGSIMP_HPP
#define _KOKKOSGSIMP_HPP

namespace KokkosKernels{

namespace Experimental{

namespace Graph{


namespace Impl{


template <typename HandleType>
class GaussSeidel{

public:


  typedef typename HandleType::idx_array_type idx_array_type;
  typedef typename HandleType::idx_edge_array_type idx_edge_array_type;
  typedef typename HandleType::value_array_type value_array_type;

  typedef typename idx_array_type::value_type idx;
  typedef typename idx_array_type::array_layout idx_array_layout;
  typedef typename idx_array_type::device_type idx_device_type;
  typedef typename idx_array_type::memory_traits idx_memory_traits;
  typedef typename idx_array_type::HostMirror host_view_type; //Host view type



  typedef typename idx_edge_array_type::value_type idx_edge;
  typedef typename idx_edge_array_type::array_layout idx_edge_array_layout;
  typedef typename idx_edge_array_type::device_type idx_edge_device_type;
  typedef typename idx_edge_array_type::memory_traits idx_edge_memory_traits;
  typedef typename idx_edge_array_type::HostMirror host_edge_view_type; //Host view type


  typedef typename value_array_type::value_type value_type;
  typedef typename value_array_type::array_layout value_type_array_layout;
  typedef typename value_array_type::device_type value_type_device_type;
  typedef typename value_array_type::memory_traits value_type_memory_traits;
  typedef typename value_array_type::HostMirror host_value_view_type; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef Kokkos::View<idx *, MyTempMemorySpace> idx_temp_work_array_type;
  typedef Kokkos::View<idx *, MyPersistentMemorySpace> idx_persistent_work_array_type;
  typedef typename idx_persistent_work_array_type::HostMirror host_idx_persistent_view_type; //Host view type


  typedef Kokkos::View<value_type *, MyTempMemorySpace> value_temp_work_array_type;
  typedef Kokkos::View<value_type *, MyPersistentMemorySpace> value_persistent_work_array_type;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef idx color_type;
  typedef Kokkos::View<value_type *, MyTempMemorySpace> color_array_type;

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy ;
  typedef typename team_policy::member_type team_member ;

private:
  HandleType *handle;
  idx num_rows, num_cols;
  idx_array_type row_map;
  idx_edge_array_type entries;
  value_array_type values;
  bool is_symmetric;
public:

  struct PSGS{

    idx_persistent_work_array_type _xadj;
    idx_persistent_work_array_type _adj; // CSR storage of the graph.
    value_persistent_work_array_type _adj_vals; // CSR storage of the graph.

    value_persistent_work_array_type _Xvector /*output*/;
    value_persistent_work_array_type _Yvector;

    PSGS(idx_persistent_work_array_type xadj_, idx_persistent_work_array_type adj_, value_persistent_work_array_type adj_vals_,
        value_persistent_work_array_type Xvector_, value_persistent_work_array_type Yvector_, idx_persistent_work_array_type color_adj_):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_),
          _Xvector( Xvector_),
          _Yvector( Yvector_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {

      idx row_begin = _xadj[ii];
      idx row_end = _xadj[ii + 1];

      value_type sum = _Yvector[ii];
      value_type diagonalVal = 1;

      for (idx adjind = row_begin; adjind < row_end; ++adjind){
        idx colIndex = _adj[adjind];
        value_type val = _adj_vals[adjind];

        if (colIndex == ii){
          diagonalVal = val;
        }
        else {
          sum -= val * _Xvector[colIndex];
        }
      }
      _Xvector[ii] = sum / diagonalVal;
    }
  };

  struct Team_PSGS{

    idx_persistent_work_array_type _xadj;
    idx_persistent_work_array_type _adj; // CSR storage of the graph.
    value_persistent_work_array_type _adj_vals; // CSR storage of the graph.

    value_persistent_work_array_type _Xvector /*output*/;
    value_persistent_work_array_type _Yvector;
    idx _color_set_begin;
    idx _color_set_end;



    Team_PSGS(idx_persistent_work_array_type xadj_, idx_persistent_work_array_type adj_, value_persistent_work_array_type adj_vals_,
        value_persistent_work_array_type Xvector_, value_persistent_work_array_type Yvector_,
        idx color_set_begin, idx color_set_end):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_),
          _Xvector( Xvector_),
          _Yvector( Yvector_),
          _color_set_begin(color_set_begin),
          _color_set_end(color_set_end){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member & teamMember) const {
      //idx ii = _color_adj[i];
      //int ii = teamMember.league_rank()  + _shift_index;

      idx ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank() + _color_set_begin;
      //check ii is out of range. if it is, just return.
      if (ii >= _color_set_end) return;



      idx row_begin = _xadj[ii];
      idx row_end = _xadj[ii + 1];

      bool am_i_the_diagonal = false;
      value_type diagonal = 1;
      value_type product = 0 ;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          //Kokkos::TeamThreadRange(teamMember, row_end - row_begin),
          [&] (idx i, value_type & valueToUpdate) {
        idx adjind = i + row_begin;
        idx colIndex = _adj[adjind];
        value_type val = _adj_vals[adjind];
        if (colIndex == ii){
          diagonal = val;
          am_i_the_diagonal = true;
        }
        else {
          valueToUpdate += val * _Xvector[colIndex];
        }
      },
      product);

      if (am_i_the_diagonal) {
        _Xvector[ii] = (_Yvector[ii] - product) / diagonal;
      }
     }
  };



  /**
   * \brief constructor
   */

  GaussSeidel(HandleType *handle_,
      idx num_rows_,
      idx num_cols_,
      idx_array_type row_map_,
      idx_edge_array_type entries_,
      value_array_type values_):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_), row_map(row_map_), entries(entries_), values(values_), is_symmetric(false){}


  GaussSeidel(HandleType *handle_,
      idx num_rows_,
      idx num_cols_,
      idx_array_type row_map_,
      idx_edge_array_type entries_):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_),
        row_map(row_map_),
        entries(entries_),
        values(), is_symmetric(false){}



  /**
   * \brief constructor
   */
  GaussSeidel(HandleType *handle_,
      idx num_rows_,
      idx num_cols_,
      idx_array_type row_map_,
      idx_edge_array_type entries_,
      value_array_type values_,
      bool is_symmetric_):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_), row_map(row_map_), entries(entries_), values(values_), is_symmetric(is_symmetric_){}




  void initialize_symbolic(){
    typename HandleType::GraphColoringHandleType *gchandle = this->handle->get_graph_coloring_handle();
    if (gchandle == NULL){

      this->handle->create_graph_coloring_handle();
      //this->handle->create_gs_handle();
      this->handle->get_gs_handle()->set_owner_of_coloring();
      gchandle = this->handle->get_graph_coloring_handle();
    }

    idx_array_type xadj = this->row_map;
    idx_edge_array_type adj = this->entries;
    idx nnz = adj.dimension_0();

    {
      idx_array_type tmp_xadj = xadj;
      idx_array_type tmp_adj = adj;
      if (!is_symmetric || num_rows < num_cols){
        KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < idx_array_type, idx_edge_array_type,
        idx_array_type, idx_edge_array_type,
        MyExecSpace>
        (num_rows, xadj, adj, tmp_xadj,tmp_adj );

      }
      graph_color_symbolic <HandleType> (this->handle, num_rows, num_rows, tmp_xadj , tmp_adj);
    }


    idx numColors = gchandle->get_num_colors();

    typename HandleType::GraphColoringHandleType::color_array_type colors =  gchandle->get_vertex_colors();

    idx_persistent_work_array_type color_xadj;

    idx_persistent_work_array_type color_adj;



    KokkosKernels::Experimental::Util::create_reverse_map
      <typename HandleType::GraphColoringHandleType::color_array_type, idx_persistent_work_array_type, MyExecSpace>
          (num_rows, numColors, colors, color_xadj, color_adj);
    MyExecSpace::fence();


    host_idx_persistent_view_type  h_color_xadj = Kokkos::create_mirror_view (color_xadj);
    Kokkos::deep_copy (h_color_xadj , color_xadj);

    //TODO: Change this to 1 sort for all matrix.
    for (idx i = 0; i < numColors; ++i){
      idx color_index_begin = h_color_xadj(i);
      idx color_index_end = h_color_xadj(i + 1);
      if (color_index_begin + 1 >= color_index_end ) continue;
      idx_persistent_work_array_type colorsubset =
          subview(color_adj, Kokkos::pair<idx, idx> (color_index_begin, color_index_end));
      Kokkos::sort (colorsubset);
    }

    idx_persistent_work_array_type permuted_xadj ("new xadj", num_rows + 1);
    idx_persistent_work_array_type old_to_new_map ("old_to_new_index_", num_rows );
    idx_persistent_work_array_type permuted_adj ("newadj_", nnz );
    Kokkos::parallel_for( my_exec_space(0,num_rows),
        create_permuted_xadj(
            color_adj,
            xadj,
            permuted_xadj,
            old_to_new_map));


    KokkosKernels::Experimental::Util::inclusive_parallel_prefix_sum
        <idx_array_type, MyExecSpace>(num_rows + 1, permuted_xadj);

    //value_persistent_work_array_type newvals_ ("newvals_", nnz );

    Kokkos::parallel_for( my_exec_space(0,num_rows),
        fill_matrix_symbolic(
            num_rows,
            color_adj,
            xadj,
            adj,
            //adj_vals,
            permuted_xadj,
            permuted_adj,
            //newvals_,
            old_to_new_map));


    typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();
    gsHandler->set_color_set_xadj(h_color_xadj);
    gsHandler->set_color_set_adj(color_adj);
    gsHandler->set_num_colors(numColors);
    gsHandler->set_new_xadj(permuted_xadj);
    gsHandler->set_new_adj(permuted_adj);
    //gsHandler->set_new_adj_val(newvals_);
    gsHandler->set_old_to_new_map(old_to_new_map);
    if (this->handle->get_gs_handle()->is_owner_of_coloring()){
      this->handle->destroy_graph_coloring_handle();
      this->handle->get_gs_handle()->set_owner_of_coloring(false);
    }
    this->handle->get_gs_handle()->set_call_symbolic(true);
    this->handle->get_gs_handle()->allocate_x_y_vectors(this->num_rows, this->num_cols);
  }

  struct create_permuted_xadj{
    idx_persistent_work_array_type color_adj;
    idx_array_type oldxadj;
    idx_persistent_work_array_type newxadj;
    idx_persistent_work_array_type old_to_new_index;
    create_permuted_xadj(
        idx_persistent_work_array_type color_adj_,
        idx_array_type oldxadj_,
        idx_persistent_work_array_type newxadj_,
        idx_persistent_work_array_type old_to_new_index_):
          color_adj(color_adj_), oldxadj(oldxadj_),
          newxadj(newxadj_),old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i) const{
      idx index = color_adj(i);
      newxadj(i + 1) = oldxadj[index + 1] - oldxadj[index];
      old_to_new_index[index] = i;
    }
  };

  struct fill_matrix_symbolic{
    idx num_rows;
    idx_persistent_work_array_type color_adj;
    idx_array_type oldxadj;
    idx_edge_array_type oldadj;
    //value_array_type oldadjvals;
    idx_persistent_work_array_type newxadj;
    idx_persistent_work_array_type newadj;
    //value_persistent_work_array_type newadjvals;
    idx_persistent_work_array_type old_to_new_index;
    fill_matrix_symbolic(
        idx num_rows_,
        idx_persistent_work_array_type color_adj_,
        idx_array_type oldxadj_,
        idx_edge_array_type oldadj_,
        //value_array_type oldadjvals_,
        idx_persistent_work_array_type newxadj_,
        idx_persistent_work_array_type newadj_,
        //value_persistent_work_array_type newadjvals_,
        idx_persistent_work_array_type old_to_new_index_):
          num_rows(num_rows_),
          color_adj(color_adj_), oldxadj(oldxadj_), oldadj(oldadj_), //oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadj(newadj_), //newadjvals(newadjvals_),
          old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i) const{
      idx index = color_adj(i);
      idx xadj_begin = newxadj(i);

      idx old_xadj_end = oldxadj[index + 1];
      for (idx j = oldxadj[index]; j < old_xadj_end; ++j){
        idx neighbor = oldadj[j];
        if(neighbor < num_rows) neighbor = old_to_new_index[neighbor];
        newadj[xadj_begin++] = neighbor;
        //newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };


  struct fill_matrix_numeric{
    idx_persistent_work_array_type color_adj;
    idx_array_type oldxadj;
    value_array_type oldadjvals;
    idx_persistent_work_array_type newxadj;
    value_persistent_work_array_type newadjvals;
    fill_matrix_numeric(
        idx_persistent_work_array_type color_adj_,
        idx_array_type oldxadj_,
        value_array_type oldadjvals_,
        idx_persistent_work_array_type newxadj_,
        value_persistent_work_array_type newadjvals_):
          color_adj(color_adj_), oldxadj(oldxadj_),  oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadjvals(newadjvals_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i) const{
      idx index = color_adj(i);
      idx xadj_begin = newxadj(i);

      idx old_xadj_end = oldxadj[index + 1];
      for (idx j = oldxadj[index]; j < old_xadj_end; ++j){
        newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };

  void initialize_numeric(){

    if (this->handle->get_gs_handle()->is_symbolic_called() == false){

      this->initialize_symbolic();

    }
    //else
    {


      idx_array_type xadj = this->row_map;
      idx_edge_array_type adj = this->entries;

      idx nnz = adj.dimension_0();
      value_array_type adj_vals = this->values;

      typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();



      idx_persistent_work_array_type newxadj_ = gsHandler->get_new_xadj();
      idx_persistent_work_array_type old_to_new_map = gsHandler->get_old_to_new_map();
      idx_persistent_work_array_type newadj_ = gsHandler->get_new_adj();

      idx_persistent_work_array_type color_adj = gsHandler->get_color_adj();
      value_persistent_work_array_type permuted_adj_vals ("newvals_", nnz );

      Kokkos::parallel_for( my_exec_space(0,num_rows),
          fill_matrix_numeric(
              color_adj,
              xadj,
              //adj,
              adj_vals,
              newxadj_,
              //newadj_,
              permuted_adj_vals
              //,old_to_new_map
              ));
      gsHandler->set_new_adj_val(permuted_adj_vals);
      this->handle->get_gs_handle()->set_call_numeric(true);

    }
  }

  void apply(
      value_array_type x_lhs_output_vec,
      value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1,
      bool apply_forward = true,
      bool apply_backward = true){
    if (this->handle->get_gs_handle()->is_numeric_called() == false){
      this->initialize_numeric();
    }

    typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();
    value_persistent_work_array_type Permuted_Yvector = gsHandler->get_permuted_y_vector();
    value_persistent_work_array_type Permuted_Xvector = gsHandler->get_permuted_x_vector();


    idx_persistent_work_array_type newxadj_ = gsHandler->get_new_xadj();
    idx_persistent_work_array_type old_to_new_map = gsHandler->get_old_to_new_map();
    idx_persistent_work_array_type newadj_ = gsHandler->get_new_adj();
    idx_persistent_work_array_type color_adj = gsHandler->get_color_adj();

    idx numColors = gsHandler->get_num_colors();



    //Kokkos::parallel_for( my_exec_space(0,nr), PermuteVector(y_rhs_input_vec, Permuted_Yvector, old_to_new_map));
    KokkosKernels::Experimental::Util::permute_vector
      <value_persistent_work_array_type, idx_persistent_work_array_type, MyExecSpace>(
        num_rows,
        old_to_new_map,
        y_rhs_input_vec,
        Permuted_Yvector
        );
    MyExecSpace::fence();
    if(init_zero_x_vector){
      KokkosKernels::Experimental::Util::zero_vector<value_persistent_work_array_type, MyExecSpace>(num_cols, Permuted_Xvector);
    }
    else{
      KokkosKernels::Experimental::Util::permute_vector<value_persistent_work_array_type, idx_persistent_work_array_type, MyExecSpace>(
          num_cols,
          old_to_new_map,
          x_lhs_output_vec,
          Permuted_Xvector
          );
    }
    MyExecSpace::fence();

    idx_persistent_work_array_type permuted_xadj = gsHandler->get_new_xadj();
    idx_persistent_work_array_type permuted_adj = gsHandler->get_new_adj();
    value_persistent_work_array_type permuted_adj_vals = gsHandler->get_new_adj_val();


    host_idx_persistent_view_type h_color_xadj = gsHandler->get_color_xadj();

    if (gsHandler->get_algorithm_type()== GS_PERMUTED){
      PSGS gs(permuted_xadj, permuted_adj, permuted_adj_vals,
          Permuted_Xvector, Permuted_Yvector, color_adj);

      this->IterativePSGS(
          gs,
          numColors,
          h_color_xadj,
          numIter,
          apply_forward,
          apply_backward);
    }
    else{

      Team_PSGS gs(permuted_xadj, permuted_adj, permuted_adj_vals,
          Permuted_Xvector, Permuted_Yvector,0,0);

      this->IterativePSGS(
          gs,
          numColors,
          h_color_xadj,
          numIter,
          apply_forward,
          apply_backward);
    }

    //Kokkos::parallel_for( my_exec_space(0,nr), PermuteVector(x_lhs_output_vec, Permuted_Xvector, color_adj));


    KokkosKernels::Experimental::Util::permute_vector<value_persistent_work_array_type, idx_persistent_work_array_type, MyExecSpace>(
        num_cols,
        color_adj,
        Permuted_Xvector,
        x_lhs_output_vec
        );
    MyExecSpace::fence();

  }

  void IterativePSGS(
      Team_PSGS &gs,
      idx numColors,
      host_view_type h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }

  void DoPSGS(Team_PSGS &gs, idx numColors, host_view_type h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    int teamSizeMax = 0;
    int vector_size = 0;
    int max_allowed_team_size = team_policy::team_size_max(gs);


    idx nnz = this->entries.dimension_0();


    this->handle->get_gs_handle()->vector_team_size(max_allowed_team_size, vector_size, teamSizeMax, num_rows, nnz);
    /*std::cout
        << "max_allowed_team_size"  << max_allowed_team_size
        << " vector_size:" << vector_size
        << " teamSizeMax:" << teamSizeMax << std::endl;
    */
    if (apply_forward){
      for (idx i = 0; i < numColors; ++i){
        idx color_index_begin = h_color_xadj(i);
        idx color_index_end = h_color_xadj(i + 1);

        int overall_work = color_index_end - color_index_begin;// /256 + 1;


        gs._color_set_begin = color_index_begin;
        gs._color_set_end = color_index_end;

        Kokkos::parallel_for(
            team_policy(overall_work / teamSizeMax + 1 , teamSizeMax, vector_size),
            gs );
        MyExecSpace::fence();
      }
    }
    if (apply_backward){
      for (idx i = numColors - 1; ; --i){
        idx color_index_begin = h_color_xadj(i);
        idx color_index_end = h_color_xadj(i + 1);

        int numberOfTeams = color_index_end - color_index_begin;// /256 + 1;
        gs._color_set_begin = color_index_begin;
        gs._color_set_end = color_index_end;

        Kokkos::parallel_for(
            team_policy(numberOfTeams / teamSizeMax + 1 , teamSizeMax, vector_size),
            gs );
        MyExecSpace::fence();
        if (i == 0){
          break;
        }
      }
    }
  }

  void IterativePSGS(
      PSGS &gs,
      idx numColors,
      host_view_type h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }



  void DoPSGS(PSGS &gs, idx numColors, host_view_type h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    if (apply_forward){
      for (idx i = 0; i < numColors; ++i){
        idx color_index_begin = h_color_xadj(i);
        idx color_index_end = h_color_xadj(i + 1);
        Kokkos::parallel_for (my_exec_space (color_index_begin, color_index_end) , gs);
        MyExecSpace::fence();
      }
    }
    if (apply_backward){
      for (idx i = numColors - 1; ; --i){
        idx color_index_begin = h_color_xadj(i);
        idx color_index_end = h_color_xadj(i + 1);
        Kokkos::parallel_for (my_exec_space (color_index_begin, color_index_end) , gs);
        MyExecSpace::fence();
        if (i == 0){
          break;
        }
      }
    }
  }
};

}
}
}
}
#endif
