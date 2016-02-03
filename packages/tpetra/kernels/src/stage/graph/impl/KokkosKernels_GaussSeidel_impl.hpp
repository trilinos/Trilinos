#include "KokkosKernels_GraphColor.hpp"
#include "KokkosKernels_Utils.hpp"
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


template <typename HandleType, typename in_row_index_view_type_, typename in_nonzero_index_view_type_, typename in_nonzero_value_view_type_>
class GaussSeidel{

public:


  typedef in_row_index_view_type_ in_row_index_view_type;
  typedef in_nonzero_index_view_type_ in_nonzero_index_view_type;
  typedef in_nonzero_value_view_type_ in_nonzero_value_view_type;

  typedef typename in_row_index_view_type::non_const_value_type row_index_type;
  typedef typename in_row_index_view_type::array_layout row_view_array_layout;
  typedef typename in_row_index_view_type::device_type row_view_device_type;
  typedef typename in_row_index_view_type::memory_traits row_view_memory_traits;
  typedef typename in_row_index_view_type::HostMirror row_host_view_type; //Host view type

  typedef typename in_nonzero_index_view_type::non_const_value_type nonzero_index_type;
  typedef typename in_nonzero_index_view_type::array_layout nonzero_index_view_array_layout;
  typedef typename in_nonzero_index_view_type::device_type nonzero_index_view_device_type;
  typedef typename in_nonzero_index_view_type::memory_traits nonzero_index_view_memory_traits;
  typedef typename in_nonzero_index_view_type::HostMirror nonzero_index_host_view_type; //Host view type


  typedef typename in_nonzero_value_view_type::non_const_value_type nonzero_value_type;
  typedef typename in_nonzero_value_view_type::array_layout nonzero_value_view_array_layout;
  typedef typename in_nonzero_value_view_type::device_type nonzero_value_view_device_type;
  typedef typename in_nonzero_value_view_type::memory_traits nonzero_value_view_memory_traits;
  typedef typename in_nonzero_value_view_type::HostMirror nonzero_value_host_view_type; //Host view type


  typedef typename in_row_index_view_type::const_data_type const_row_data_type;
  typedef typename in_row_index_view_type::non_const_data_type non_const_row_data_type;
  typedef typename in_row_index_view_type::memory_space row_view_memory_space;
  typedef typename Kokkos::View<const_row_data_type, row_view_array_layout,
      row_view_memory_space, row_view_memory_traits> const_row_index_view_type;
  typedef typename Kokkos::View<non_const_row_data_type, row_view_array_layout,
      row_view_memory_space, row_view_memory_traits> non_const_row_index_view_type;

  typedef typename in_nonzero_index_view_type::const_data_type const_nonzero_index_data_type;
  typedef typename in_nonzero_index_view_type::non_const_data_type non_const_nonzero_index_data_type;
  typedef typename in_nonzero_index_view_type::memory_space nonzero_index_view_memory_space;
  typedef typename Kokkos::View<const_nonzero_index_data_type, nonzero_index_view_array_layout,
      nonzero_index_view_memory_space, nonzero_index_view_memory_traits> const_nonzero_index_view_type;
  typedef typename Kokkos::View<non_const_nonzero_index_data_type, nonzero_index_view_array_layout,
      nonzero_index_view_memory_space, nonzero_index_view_memory_traits> non_const_nonzero_index_view_type;

  typedef typename in_nonzero_value_view_type::const_data_type const_nonzero_value_data_type;
  typedef typename in_nonzero_value_view_type::non_const_data_type non_const_nonzero_value_data_type;
  typedef typename in_nonzero_value_view_type::memory_space nonzero_value_view_memory_space;
  typedef typename Kokkos::View<const_nonzero_value_data_type, nonzero_value_view_array_layout,
      nonzero_value_view_memory_space, nonzero_value_view_memory_traits> const_nonzero_value_view_type;
  typedef typename Kokkos::View<non_const_nonzero_value_data_type, nonzero_value_view_array_layout,
      nonzero_value_view_memory_space, nonzero_value_view_memory_traits> non_const_nonzero_value_view_type;


  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef typename HandleType::row_index_temp_work_view_type row_index_temp_work_view_type;
  typedef typename HandleType::row_index_persistent_work_view_type row_index_persistent_work_view_type;
  typedef typename HandleType::row_index_persistent_host_view_type row_index_persistent_host_view_type; //Host view type


  typedef typename HandleType::nonzero_value_temp_work_view_type nonzero_value_temp_work_view_type;
  typedef typename HandleType::nonzero_value_persistent_work_view_type nonzero_value_persistent_work_view_type;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef row_index_type color_type;
  typedef Kokkos::View<row_index_type *, MyTempMemorySpace> color_view_type;

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy ;
  typedef typename team_policy::member_type team_member ;

private:
  HandleType *handle;
  row_index_type num_rows, num_cols;

  const_row_index_view_type row_map;
  const_nonzero_index_view_type entries;
  const_nonzero_value_view_type values;
  bool is_symmetric;
public:

  struct PSGS{

    row_index_persistent_work_view_type _xadj;
    row_index_persistent_work_view_type _adj; // CSR storage of the graph.
    nonzero_value_persistent_work_view_type _adj_vals; // CSR storage of the graph.

    nonzero_value_persistent_work_view_type _Xvector /*output*/;
    nonzero_value_persistent_work_view_type _Yvector;

    PSGS(row_index_persistent_work_view_type xadj_, row_index_persistent_work_view_type adj_, nonzero_value_persistent_work_view_type adj_vals_,
        nonzero_value_persistent_work_view_type Xvector_, nonzero_value_persistent_work_view_type Yvector_, row_index_persistent_work_view_type color_adj_):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_),
          _Xvector( Xvector_),
          _Yvector( Yvector_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &ii) const {

      row_index_type row_begin = _xadj[ii];
      row_index_type row_end = _xadj[ii + 1];

      nonzero_value_type sum = _Yvector[ii];
      nonzero_value_type diagonalVal = 1;

      for (row_index_type adjind = row_begin; adjind < row_end; ++adjind){
        row_index_type colIndex = _adj[adjind];
        nonzero_value_type val = _adj_vals[adjind];

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

    row_index_persistent_work_view_type _xadj;
    row_index_persistent_work_view_type _adj; // CSR storage of the graph.
    nonzero_value_persistent_work_view_type _adj_vals; // CSR storage of the graph.

    nonzero_value_persistent_work_view_type _Xvector /*output*/;
    nonzero_value_persistent_work_view_type _Yvector;
    row_index_type _color_set_begin;
    row_index_type _color_set_end;



    Team_PSGS(row_index_persistent_work_view_type xadj_, row_index_persistent_work_view_type adj_, nonzero_value_persistent_work_view_type adj_vals_,
        nonzero_value_persistent_work_view_type Xvector_, nonzero_value_persistent_work_view_type Yvector_,
        row_index_type color_set_begin, row_index_type color_set_end):
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

      row_index_type ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank() + _color_set_begin;
      //check ii is out of range. if it is, just return.
      if (ii >= _color_set_end) return;



      row_index_type row_begin = _xadj[ii];
      row_index_type row_end = _xadj[ii + 1];

      bool am_i_the_diagonal = false;
      nonzero_value_type diagonal = 1;
      nonzero_value_type product = 0 ;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          //Kokkos::TeamThreadRange(teamMember, row_end - row_begin),
          [&] (row_index_type i, nonzero_value_type & valueToUpdate) {
        row_index_type adjind = i + row_begin;
        row_index_type colIndex = _adj[adjind];
        nonzero_value_type val = _adj_vals[adjind];
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
      row_index_type num_rows_,
      row_index_type num_cols_,
      const_row_index_view_type row_map_,
      const_nonzero_index_view_type entries_,
      const_nonzero_value_view_type values_):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_), row_map(row_map_), entries(entries_), values(values_), is_symmetric(true){}


  GaussSeidel(HandleType *handle_,
      row_index_type num_rows_,
      row_index_type num_cols_,
      const_row_index_view_type row_map_,
      const_nonzero_index_view_type entries_,
      bool is_symmetric_ = true):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_),
        row_map(row_map_),
        entries(entries_),
        values(), is_symmetric(is_symmetric_){}



  /**
   * \brief constructor
   */
  GaussSeidel(HandleType *handle_,
      row_index_type num_rows_,
      row_index_type num_cols_,
      const_row_index_view_type row_map_,
      const_nonzero_index_view_type entries_,
      const_nonzero_value_view_type values_,
      bool is_symmetric_):
        num_rows(num_rows_), num_cols(num_cols_),
        handle(handle_), row_map(row_map_), entries(entries_), values(values_), is_symmetric(is_symmetric_){}




  void initialize_symbolic(){
    //std::cout << std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
    typename HandleType::GraphColoringHandleType *gchandle = this->handle->get_graph_coloring_handle();
    if (gchandle == NULL){

      this->handle->create_graph_coloring_handle();
      //this->handle->create_gs_handle();
      this->handle->get_gs_handle()->set_owner_of_coloring();
      gchandle = this->handle->get_graph_coloring_handle();
    }



    const_row_index_view_type xadj = this->row_map;
    const_nonzero_index_view_type adj = this->entries;
    row_index_type nnz = adj.dimension_0();
    //std::cout << "initialize_symbolic start" << std::endl;
    {
      if (!is_symmetric){

        if (gchandle->get_coloring_type() == KokkosKernels::Experimental::Graph::COLORING_EB){

          gchandle->symmetrize_and_calculate_lower_diagonal_edge_list(num_rows, xadj, adj);
          graph_color_symbolic <HandleType, const_row_index_view_type, const_nonzero_index_view_type>
              (this->handle, num_rows, num_rows, xadj , adj);
        }
        else {

          row_index_temp_work_view_type tmp_xadj;
          row_index_temp_work_view_type tmp_adj;

          KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap
          < const_row_index_view_type, const_nonzero_index_view_type,
          row_index_temp_work_view_type, row_index_temp_work_view_type,
          MyExecSpace>
          (num_rows, xadj, adj, tmp_xadj, tmp_adj );
          //std::cout << "symmetrize_graph_symbolic " << std::endl;

          graph_color_symbolic <HandleType, row_index_temp_work_view_type, row_index_temp_work_view_type> (this->handle, num_rows, num_rows, tmp_xadj , tmp_adj);
        }
      }
      //std::cout << "graph_color_symbolic STARTTTTTT num_rows:" << num_rows  << " tmp_adj:" << tmp_adj.dimension_0() << " adj:" << adj.dimension_0()<< std::endl;

      else {

        graph_color_symbolic <HandleType, const_row_index_view_type, const_nonzero_index_view_type> (this->handle, num_rows, num_rows, xadj , adj);
      }
      //std::cout << "graph_color_symbolic ENDDDDDDDDDD num_rows:" << num_rows << std::endl;
    }

    row_index_type numColors = gchandle->get_num_colors();


    typename HandleType::GraphColoringHandleType::color_view_type colors =  gchandle->get_vertex_colors();

    row_index_persistent_work_view_type color_xadj;

    row_index_persistent_work_view_type color_adj;



    KokkosKernels::Experimental::Util::create_reverse_map
      <typename HandleType::GraphColoringHandleType::color_view_type,
        row_index_persistent_work_view_type, MyExecSpace>
          (num_rows, numColors, colors, color_xadj, color_adj);
    MyExecSpace::fence();
    //std::cout << "create_reverse_map" << std::endl;

    row_index_persistent_host_view_type  h_color_xadj = Kokkos::create_mirror_view (color_xadj);
    Kokkos::deep_copy (h_color_xadj , color_xadj);

    //TODO: Change this to 1 sort for all matrix.
    for (row_index_type i = 0; i < numColors; ++i){
      row_index_type color_index_begin = h_color_xadj(i);
      row_index_type color_index_end = h_color_xadj(i + 1);
      if (color_index_begin + 1 >= color_index_end ) continue;
      row_index_persistent_work_view_type colorsubset =
          subview(color_adj, Kokkos::pair<row_index_type, row_index_type> (color_index_begin, color_index_end));
      Kokkos::sort (colorsubset);
    }
    //std::cout << "sort" << std::endl;

    row_index_persistent_work_view_type permuted_xadj ("new xadj", num_rows + 1);
    row_index_persistent_work_view_type old_to_new_map ("old_to_new_index_", num_rows );
    row_index_persistent_work_view_type permuted_adj ("newadj_", nnz );
    Kokkos::parallel_for( my_exec_space(0,num_rows),
        create_permuted_xadj(
            color_adj,
            xadj,
            permuted_xadj,
            old_to_new_map));
    //std::cout << "create_permuted_xadj" << std::endl;


    KokkosKernels::Experimental::Util::inclusive_parallel_prefix_sum
        <row_index_persistent_work_view_type, MyExecSpace>
        (num_rows + 1, permuted_xadj);

    //std::cout << "inclusive_parallel_prefix_sum" << std::endl;

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

    //std::cout << "fill_matrix_symbolic" << std::endl;


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
    //std::cout << "all end" << std::endl;
  }

  struct create_permuted_xadj{
    row_index_persistent_work_view_type color_adj;
    const_row_index_view_type oldxadj;
    row_index_persistent_work_view_type newxadj;
    row_index_persistent_work_view_type old_to_new_index;
    create_permuted_xadj(
        row_index_persistent_work_view_type color_adj_,
        const_row_index_view_type oldxadj_,
        row_index_persistent_work_view_type newxadj_,
        row_index_persistent_work_view_type old_to_new_index_):
          color_adj(color_adj_), oldxadj(oldxadj_),
          newxadj(newxadj_),old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &i) const{
      row_index_type index = color_adj(i);
      newxadj(i + 1) = oldxadj[index + 1] - oldxadj[index];
      old_to_new_index[index] = i;
    }
  };

  struct fill_matrix_symbolic{
    row_index_type num_rows;
    row_index_persistent_work_view_type color_adj;
    const_row_index_view_type oldxadj;
    const_nonzero_index_view_type oldadj;
    //value_array_type oldadjvals;
    row_index_persistent_work_view_type newxadj;
    row_index_persistent_work_view_type newadj;
    //value_persistent_work_array_type newadjvals;
    row_index_persistent_work_view_type old_to_new_index;
    fill_matrix_symbolic(
        row_index_type num_rows_,
        row_index_persistent_work_view_type color_adj_,
        const_row_index_view_type oldxadj_,
        const_nonzero_index_view_type oldadj_,
        //value_array_type oldadjvals_,
        row_index_persistent_work_view_type newxadj_,
        row_index_persistent_work_view_type newadj_,
        //value_persistent_work_array_type newadjvals_,
        row_index_persistent_work_view_type old_to_new_index_):
          num_rows(num_rows_),
          color_adj(color_adj_), oldxadj(oldxadj_), oldadj(oldadj_), //oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadj(newadj_), //newadjvals(newadjvals_),
          old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &i) const{
      row_index_type index = color_adj(i);
      row_index_type xadj_begin = newxadj(i);

      row_index_type old_xadj_end = oldxadj[index + 1];
      for (row_index_type j = oldxadj[index]; j < old_xadj_end; ++j){
        row_index_type neighbor = oldadj[j];
        if(neighbor < num_rows) neighbor = old_to_new_index[neighbor];
        newadj[xadj_begin++] = neighbor;
        //newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };


  struct fill_matrix_numeric{
    row_index_persistent_work_view_type color_adj;
    const_row_index_view_type oldxadj;
    const_nonzero_value_view_type oldadjvals;
    row_index_persistent_work_view_type newxadj;
    nonzero_value_persistent_work_view_type newadjvals;
    fill_matrix_numeric(
        row_index_persistent_work_view_type color_adj_,
        const_row_index_view_type oldxadj_,
        const_nonzero_value_view_type oldadjvals_,
        row_index_persistent_work_view_type newxadj_,
        nonzero_value_persistent_work_view_type newadjvals_):
          color_adj(color_adj_), oldxadj(oldxadj_),  oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadjvals(newadjvals_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &i) const{
      row_index_type index = color_adj(i);
      row_index_type xadj_begin = newxadj(i);

      row_index_type old_xadj_end = oldxadj[index + 1];
      for (row_index_type j = oldxadj[index]; j < old_xadj_end; ++j){
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


      const_row_index_view_type xadj = this->row_map;
      const_nonzero_index_view_type adj = this->entries;

      nonzero_index_type nnz = adj.dimension_0();
      const_nonzero_value_view_type adj_vals = this->values;

      typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();



      row_index_persistent_work_view_type newxadj_ = gsHandler->get_new_xadj();
      row_index_persistent_work_view_type old_to_new_map = gsHandler->get_old_to_new_map();
      row_index_persistent_work_view_type newadj_ = gsHandler->get_new_adj();

      row_index_persistent_work_view_type color_adj = gsHandler->get_color_adj();
      nonzero_value_persistent_work_view_type permuted_adj_vals ("newvals_", nnz );

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

  template <typename x_value_array_type, typename y_value_array_type>
  void apply(
      x_value_array_type x_lhs_output_vec,
      y_value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1,
      bool apply_forward = true,
      bool apply_backward = true,
      bool update_y_vector = true){
    if (this->handle->get_gs_handle()->is_numeric_called() == false){
      this->initialize_numeric();
    }

    typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();
    nonzero_value_persistent_work_view_type Permuted_Yvector = gsHandler->get_permuted_y_vector();
    nonzero_value_persistent_work_view_type Permuted_Xvector = gsHandler->get_permuted_x_vector();


    row_index_persistent_work_view_type newxadj_ = gsHandler->get_new_xadj();
    row_index_persistent_work_view_type old_to_new_map = gsHandler->get_old_to_new_map();
    row_index_persistent_work_view_type newadj_ = gsHandler->get_new_adj();
    row_index_persistent_work_view_type color_adj = gsHandler->get_color_adj();

    row_index_type numColors = gsHandler->get_num_colors();



    if (update_y_vector){
      KokkosKernels::Experimental::Util::permute_vector
        <y_value_array_type,
        nonzero_value_persistent_work_view_type,
        row_index_persistent_work_view_type, MyExecSpace>(
          num_rows,
          old_to_new_map,
          y_rhs_input_vec,
          Permuted_Yvector
      );
    }
    MyExecSpace::fence();
    if(init_zero_x_vector){
      KokkosKernels::Experimental::Util::zero_vector<
          nonzero_value_persistent_work_view_type, MyExecSpace>(num_cols, Permuted_Xvector);
    }
    else{
      KokkosKernels::Experimental::Util::permute_vector
        <x_value_array_type, nonzero_value_persistent_work_view_type, row_index_persistent_work_view_type, MyExecSpace>(
          num_cols,
          old_to_new_map,
          x_lhs_output_vec,
          Permuted_Xvector
          );
    }
    MyExecSpace::fence();

    row_index_persistent_work_view_type permuted_xadj = gsHandler->get_new_xadj();
    row_index_persistent_work_view_type permuted_adj = gsHandler->get_new_adj();
    nonzero_value_persistent_work_view_type permuted_adj_vals = gsHandler->get_new_adj_val();


    row_index_persistent_host_view_type h_color_xadj = gsHandler->get_color_xadj();

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


    KokkosKernels::Experimental::Util::permute_vector
    <nonzero_value_persistent_work_view_type,x_value_array_type,  row_index_persistent_work_view_type, MyExecSpace>(
        num_cols,
        color_adj,
        Permuted_Xvector,
        x_lhs_output_vec
        );
    MyExecSpace::fence();

  }

  void IterativePSGS(
      Team_PSGS &gs,
      row_index_type numColors,
      row_index_persistent_host_view_type h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }

  void DoPSGS(Team_PSGS &gs, row_index_type numColors, row_index_persistent_host_view_type h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    int teamSizeMax = 0;
    int vector_size = 0;
    int max_allowed_team_size = team_policy::team_size_max(gs);


    row_index_type nnz = this->entries.dimension_0();


    this->handle->get_gs_handle()->vector_team_size(max_allowed_team_size, vector_size, teamSizeMax, num_rows, nnz);
    /*std::cout
        << "max_allowed_team_size"  << max_allowed_team_size
        << " vector_size:" << vector_size
        << " teamSizeMax:" << teamSizeMax << std::endl;
    */
    if (apply_forward){
      for (row_index_type i = 0; i < numColors; ++i){
        row_index_type color_index_begin = h_color_xadj(i);
        row_index_type color_index_end = h_color_xadj(i + 1);

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
      if (numColors > 0)
      for (row_index_type i = numColors - 1;  ; --i){
        row_index_type color_index_begin = h_color_xadj(i);
        row_index_type color_index_end = h_color_xadj(i + 1);

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
      row_index_type numColors,
      row_index_persistent_host_view_type h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }



  void DoPSGS(PSGS &gs, row_index_type numColors, row_index_persistent_host_view_type h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    if (apply_forward){
      for (row_index_type i = 0; i < numColors; ++i){
        row_index_type color_index_begin = h_color_xadj(i);
        row_index_type color_index_end = h_color_xadj(i + 1);
        Kokkos::parallel_for (my_exec_space (color_index_begin, color_index_end) , gs);
        MyExecSpace::fence();
      }
    }
    if (apply_backward){
      for (row_index_type i = numColors - 1; ; --i){
        row_index_type color_index_begin = h_color_xadj(i);
        row_index_type color_index_end = h_color_xadj(i + 1);
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
