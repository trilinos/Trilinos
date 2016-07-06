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


template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
class GaussSeidel{

public:


  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename in_lno_row_view_t::non_const_value_type row_lno_t;
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


  typedef typename in_lno_row_view_t::const_data_type const_row_lno_t;
  typedef typename in_lno_row_view_t::non_const_data_type non_const_row_lno_t;
  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;

  typedef typename in_lno_nnz_view_t::const_data_type const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::non_const_data_type non_const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;
  typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;

  typedef typename in_scalar_nnz_view_t::const_data_type const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::non_const_data_type non_const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::const_type const_scalar_nnz_view_t;
  typedef typename in_scalar_nnz_view_t::non_const_type non_const_scalar_nnz_view_t;


  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t; //Host view type


  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef row_lno_t color_t;
  typedef Kokkos::View<row_lno_t *, MyTempMemorySpace> color_view_t;

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t ;
  typedef typename team_policy_t::member_type team_member_t ;

private:
  HandleType *handle;
  row_lno_t num_rows, num_cols;

  const_lno_row_view_t row_map;
  const_lno_nnz_view_t entries;
  const_scalar_nnz_view_t values;
  bool is_symmetric;
public:

  struct PSGS{

    row_lno_persistent_work_view_t _xadj;
    row_lno_persistent_work_view_t _adj; // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals; // CSR storage of the graph.

    scalar_persistent_work_view_t _Xvector /*output*/;
    scalar_persistent_work_view_t _Yvector;

    scalar_persistent_work_view_t _permuted_diagonals;

    PSGS(row_lno_persistent_work_view_t xadj_, row_lno_persistent_work_view_t adj_, scalar_persistent_work_view_t adj_vals_,
        scalar_persistent_work_view_t Xvector_, scalar_persistent_work_view_t Yvector_, row_lno_persistent_work_view_t color_adj_,
        scalar_persistent_work_view_t permuted_diagonals_):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_),
          _Xvector( Xvector_),
          _Yvector( Yvector_), _permuted_diagonals(permuted_diagonals_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &ii) const {

      row_lno_t row_begin = _xadj[ii];
      row_lno_t row_end = _xadj[ii + 1];

      nnz_scalar_t sum = _Yvector[ii];

      for (row_lno_t adjind = row_begin; adjind < row_end; ++adjind){
        row_lno_t colIndex = _adj[adjind];
        nnz_scalar_t val = _adj_vals[adjind];
        sum -= val * _Xvector[colIndex];
      }
      nnz_scalar_t diagonalVal = _permuted_diagonals[ii];
      _Xvector[ii] = (sum + diagonalVal * _Xvector[ii])/ diagonalVal;
    }
  };

  struct Team_PSGS{

    row_lno_persistent_work_view_t _xadj;
    row_lno_persistent_work_view_t _adj; // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals; // CSR storage of the graph.

    scalar_persistent_work_view_t _Xvector /*output*/;
    scalar_persistent_work_view_t _Yvector;
    row_lno_t _color_set_begin;
    row_lno_t _color_set_end;

    scalar_persistent_work_view_t _permuted_diagonals;


    Team_PSGS(row_lno_persistent_work_view_t xadj_, row_lno_persistent_work_view_t adj_, scalar_persistent_work_view_t adj_vals_,
        scalar_persistent_work_view_t Xvector_, scalar_persistent_work_view_t Yvector_,
        row_lno_t color_set_begin, row_lno_t color_set_end,
        scalar_persistent_work_view_t permuted_diagonals_):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_),
          _Xvector( Xvector_),
          _Yvector( Yvector_),
          _color_set_begin(color_set_begin),
          _color_set_end(color_set_end), _permuted_diagonals(permuted_diagonals_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {
      //idx ii = _color_adj[i];
      //int ii = teamMember.league_rank()  + _shift_index;

      row_lno_t ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank() + _color_set_begin;
      //check ii is out of range. if it is, just return.
      if (ii >= _color_set_end)
        return;



      row_lno_t row_begin = _xadj[ii];
      row_lno_t row_end = _xadj[ii + 1];

      //bool am_i_the_diagonal = false;
      //nnz_scalar_t diagonal = 1;
      nnz_scalar_t product = 0 ;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          //Kokkos::TeamThreadRange(teamMember, row_end - row_begin),
          [&] (row_lno_t i, nnz_scalar_t & valueToUpdate) {
        row_lno_t adjind = i + row_begin;
        row_lno_t colIndex = _adj[adjind];
        nnz_scalar_t val = _adj_vals[adjind];
        /*
        if (colIndex == ii){
          diagonal = val;
          am_i_the_diagonal = true;
        }
        else {
          valueToUpdate += val * _Xvector[colIndex];
        }
        */
        valueToUpdate += val * _Xvector[colIndex];
      },
      product);

      /*
      if (am_i_the_diagonal) {
        _Xvector[ii] = (_Yvector[ii] - product) / diagonal;
      }*/

      Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
        nnz_scalar_t diagonalVal = _permuted_diagonals[ii];
        _Xvector[ii] = (_Yvector[ii] - product + diagonalVal * _Xvector[ii])/ diagonalVal;
      });

     }
  };



  /**
   * \brief constructor
   */

  GaussSeidel(HandleType *handle_,
      row_lno_t num_rows_,
      row_lno_t num_cols_,
      const_lno_row_view_t row_map_,
      const_lno_nnz_view_t entries_,
      const_scalar_nnz_view_t values_):
        handle(handle_), num_rows(num_rows_), num_cols(num_cols_),
        row_map(row_map_), entries(entries_), values(values_), is_symmetric(true){}


  GaussSeidel(HandleType *handle_,
      row_lno_t num_rows_,
      row_lno_t num_cols_,
      const_lno_row_view_t row_map_,
      const_lno_nnz_view_t entries_,
      bool is_symmetric_ = true):
        handle(handle_),
        num_rows(num_rows_), num_cols(num_cols_),

        row_map(row_map_),
        entries(entries_),
        values(), is_symmetric(is_symmetric_){}



  /**
   * \brief constructor
   */
  GaussSeidel(HandleType *handle_,
      row_lno_t num_rows_,
      row_lno_t num_cols_,
      const_lno_row_view_t row_map_,
      const_lno_nnz_view_t entries_,
      const_scalar_nnz_view_t values_,
      bool is_symmetric_):
        handle(handle_),
        num_rows(num_rows_), num_cols(num_cols_),
        row_map(row_map_), entries(entries_), values(values_), is_symmetric(is_symmetric_){}




  void initialize_symbolic(){
    //std::cout << std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
    typename HandleType::GraphColoringHandleType *gchandle = this->handle->get_graph_coloring_handle();


    if (gchandle == NULL){

      this->handle->create_graph_coloring_handle();
      //this->handle->create_gs_handle();
      this->handle->get_gs_handle()->set_owner_of_coloring();
      gchandle = this->handle->get_graph_coloring_handle();
    }



    const_lno_row_view_t xadj = this->row_map;
    const_lno_nnz_view_t adj = this->entries;
    row_lno_t nnz = adj.dimension_0();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    Kokkos::Impl::Timer timer;
#endif
    {
      if (!is_symmetric){

        if (gchandle->get_coloring_algo_type() == KokkosKernels::Experimental::Graph::COLORING_EB){
	 
          gchandle->symmetrize_and_calculate_lower_diagonal_edge_list(num_rows, xadj, adj);
          graph_color_symbolic <HandleType, const_lno_row_view_t, const_lno_nnz_view_t>
              (this->handle, num_rows, num_rows, xadj , adj);
        }
        else {
          row_lno_temp_work_view_t tmp_xadj;
          row_lno_temp_work_view_t tmp_adj;
          KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap
          < const_lno_row_view_t, const_lno_nnz_view_t,
          row_lno_temp_work_view_t, row_lno_temp_work_view_t,
          MyExecSpace>
          (num_rows, xadj, adj, tmp_xadj, tmp_adj );


          graph_color_symbolic <HandleType, row_lno_temp_work_view_t, row_lno_temp_work_view_t> (this->handle, num_rows, num_rows, tmp_xadj , tmp_adj);
        }
      }
      else {
        graph_color_symbolic <HandleType, const_lno_row_view_t, const_lno_nnz_view_t> (this->handle, num_rows, num_rows, xadj , adj);
      }
    }
    row_lno_t numColors = gchandle->get_num_colors();
   //std::cout << "numCol:" << numColors << " numRows:" << num_rows << " cols:" << num_cols << " nnz:" << adj.dimension_0() <<  std::endl;

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "COLORING_TIME:" << timer.seconds() << std::endl;
#endif


    typename HandleType::GraphColoringHandleType::color_view_t colors =  gchandle->get_vertex_colors();

    row_lno_persistent_work_view_t color_xadj;

    row_lno_persistent_work_view_t color_adj;


#ifdef KOKKOSKERNELS_TIME_REVERSE
    timer.reset();
#endif

    KokkosKernels::Experimental::Util::create_reverse_map
      <typename HandleType::GraphColoringHandleType::color_view_t,
        row_lno_persistent_work_view_t, MyExecSpace>
          (num_rows, numColors, colors, color_xadj, color_adj);
    MyExecSpace::fence();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "CREATE_REVERSE_MAP:" << timer.seconds() << std::endl;
    timer.reset();
#endif

    row_lno_persistent_work_host_view_t  h_color_xadj = Kokkos::create_mirror_view (color_xadj);
    Kokkos::deep_copy (h_color_xadj , color_xadj);
    MyExecSpace::fence();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "DEEP_COPY:" << timer.seconds() << std::endl;
    timer.reset();
#endif


#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, MyExecSpace >::value){
      for (row_lno_t i = 0; i < numColors; ++i){
        row_lno_t color_index_begin = h_color_xadj(i);
        row_lno_t color_index_end = h_color_xadj(i + 1);
        if (color_index_begin + 1 >= color_index_end ) continue;
        row_lno_persistent_work_view_t colorsubset =
            subview(color_adj, Kokkos::pair<row_lno_t, row_lno_t> (color_index_begin, color_index_end));
        Kokkos::sort (colorsubset);
      }
    }
#endif

    MyExecSpace::fence();
#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "SORT_TIME:" << timer.seconds() << std::endl;
    timer.reset();
    //std::cout << "sort" << std::endl;
#endif

    row_lno_persistent_work_view_t permuted_xadj ("new xadj", num_rows + 1);
    row_lno_persistent_work_view_t old_to_new_map ("old_to_new_index_", num_rows );
    row_lno_persistent_work_view_t permuted_adj ("newadj_", nnz );
    Kokkos::parallel_for( my_exec_space(0,num_rows),
        create_permuted_xadj(
            color_adj,
            xadj,
            permuted_xadj,
            old_to_new_map));
    //std::cout << "create_permuted_xadj" << std::endl;
    MyExecSpace::fence();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "CREATE_PERMUTED_XADJ:" << timer.seconds() << std::endl;

    timer.reset();
#endif


    KokkosKernels::Experimental::Util::inclusive_parallel_prefix_sum
        <row_lno_persistent_work_view_t, MyExecSpace>
        (num_rows + 1, permuted_xadj);
    MyExecSpace::fence();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "INCLUSIVE_PPS:" << timer.seconds() << std::endl;
    timer.reset();
#endif


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
    MyExecSpace::fence();

#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "SYMBOLIC_FILL:" << timer.seconds() << std::endl;
    timer.reset();
#endif



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
#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "ALLOC:" << timer.seconds() << std::endl;
#endif
  }

  struct create_permuted_xadj{
    row_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    row_lno_persistent_work_view_t newxadj;
    row_lno_persistent_work_view_t old_to_new_index;
    create_permuted_xadj(
        row_lno_persistent_work_view_t color_adj_,
        const_lno_row_view_t oldxadj_,
        row_lno_persistent_work_view_t newxadj_,
        row_lno_persistent_work_view_t old_to_new_index_):
          color_adj(color_adj_), oldxadj(oldxadj_),
          newxadj(newxadj_),old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i) const{
      row_lno_t index = color_adj(i);
      newxadj(i + 1) = oldxadj[index + 1] - oldxadj[index];
      old_to_new_index[index] = i;
    }
  };

  struct fill_matrix_symbolic{
    row_lno_t num_rows;
    row_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    const_lno_nnz_view_t oldadj;
    //value_array_type oldadjvals;
    row_lno_persistent_work_view_t newxadj;
    row_lno_persistent_work_view_t newadj;
    //value_persistent_work_array_type newadjvals;
    row_lno_persistent_work_view_t old_to_new_index;
    fill_matrix_symbolic(
        row_lno_t num_rows_,
        row_lno_persistent_work_view_t color_adj_,
        const_lno_row_view_t oldxadj_,
        const_lno_nnz_view_t oldadj_,
        //value_array_type oldadjvals_,
        row_lno_persistent_work_view_t newxadj_,
        row_lno_persistent_work_view_t newadj_,
        //value_persistent_work_array_type newadjvals_,
        row_lno_persistent_work_view_t old_to_new_index_):
          num_rows(num_rows_),
          color_adj(color_adj_), oldxadj(oldxadj_), oldadj(oldadj_), //oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadj(newadj_), //newadjvals(newadjvals_),
          old_to_new_index(old_to_new_index_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i) const{
      row_lno_t index = color_adj(i);
      row_lno_t xadj_begin = newxadj(i);

      row_lno_t old_xadj_end = oldxadj[index + 1];
      for (row_lno_t j = oldxadj[index]; j < old_xadj_end; ++j){
        row_lno_t neighbor = oldadj[j];
        if(neighbor < num_rows) neighbor = old_to_new_index[neighbor];
        newadj[xadj_begin++] = neighbor;
        //newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };


  struct fill_matrix_numeric{
    row_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    const_scalar_nnz_view_t oldadjvals;
    row_lno_persistent_work_view_t newxadj;
    scalar_persistent_work_view_t newadjvals;
    fill_matrix_numeric(
        row_lno_persistent_work_view_t color_adj_,
        const_lno_row_view_t oldxadj_,
        const_scalar_nnz_view_t oldadjvals_,
        row_lno_persistent_work_view_t newxadj_,
        scalar_persistent_work_view_t newadjvals_):
          color_adj(color_adj_), oldxadj(oldxadj_),  oldadjvals(oldadjvals_),
          newxadj(newxadj_), newadjvals(newadjvals_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i) const{
      row_lno_t index = color_adj(i);
      row_lno_t xadj_begin = newxadj(i);

      row_lno_t old_xadj_end = oldxadj[index + 1];
      for (row_lno_t j = oldxadj[index]; j < old_xadj_end; ++j){
        newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };


  struct Get_Matrix_Diagonals{

    row_lno_persistent_work_view_t _xadj;
    row_lno_persistent_work_view_t _adj; // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals; // CSR storage of the graph.
    scalar_persistent_work_view_t _diagonals;
    row_lno_t nr;

    Get_Matrix_Diagonals(
        row_lno_persistent_work_view_t xadj_,
        row_lno_persistent_work_view_t adj_,
        scalar_persistent_work_view_t adj_vals_,
        scalar_persistent_work_view_t diagonals_):
          _xadj( xadj_),
          _adj( adj_),
          _adj_vals( adj_vals_), _diagonals(diagonals_),
          nr(xadj_.dimension_0() - 1){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t & ii) const {
      row_lno_t row_begin = _xadj[ii];
      row_lno_t row_end = _xadj[ii + 1];
      for (row_lno_t c = row_begin; c < row_end; ++c){
        row_lno_t colIndex = _adj[c];
        if (colIndex == ii){
          nnz_scalar_t val = _adj_vals[c];
          _diagonals[ii] = val;
        }
      }
    }
  };

  void initialize_numeric(){

    if (this->handle->get_gs_handle()->is_symbolic_called() == false){
      this->initialize_symbolic();
    }
    //else
#ifdef KOKKOSKERNELS_TIME_REVERSE
    Kokkos::Impl::Timer timer;
#endif
    {


      const_lno_row_view_t xadj = this->row_map;
      const_lno_nnz_view_t adj = this->entries;

      nnz_lno_t nnz = adj.dimension_0();
      const_scalar_nnz_view_t adj_vals = this->values;

      typename HandleType::GaussSeidelHandleType *gsHandler = this->handle->get_gs_handle();



      row_lno_persistent_work_view_t newxadj_ = gsHandler->get_new_xadj();
      row_lno_persistent_work_view_t old_to_new_map = gsHandler->get_old_to_new_map();
      row_lno_persistent_work_view_t newadj_ = gsHandler->get_new_adj();

      row_lno_persistent_work_view_t color_adj = gsHandler->get_color_adj();
      scalar_persistent_work_view_t permuted_adj_vals (Kokkos::ViewAllocateWithoutInitializing("newvals_"), nnz );

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
      MyExecSpace::fence();
      gsHandler->set_new_adj_val(permuted_adj_vals);



      scalar_persistent_work_view_t permuted_diagonals (Kokkos::ViewAllocateWithoutInitializing("permuted_diagonals"), num_rows );

      Get_Matrix_Diagonals gmd(newxadj_, newadj_, permuted_adj_vals, permuted_diagonals);
      /*
      int teamSizeMax = 0;
      int vector_size = 0;
      int max_allowed_team_size = team_policy_t::team_size_max(gmd);

      this->handle->get_gs_handle()->vector_team_size(max_allowed_team_size, vector_size, teamSizeMax, num_rows, nnz);
      Kokkos::parallel_for(
          team_policy_t(num_rows / teamSizeMax + 1 , teamSizeMax, vector_size),
          gmd );
          */
      Kokkos::parallel_for(
                my_exec_space(0,num_rows),
                gmd );
      MyExecSpace::fence();
      this->handle->get_gs_handle()->set_permuted_diagonals(permuted_diagonals);


      this->handle->get_gs_handle()->set_call_numeric(true);

    }
#ifdef KOKKOSKERNELS_TIME_REVERSE
    std::cout << "NUMERIC:" << timer.seconds() << std::endl;
#endif
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
    scalar_persistent_work_view_t Permuted_Yvector = gsHandler->get_permuted_y_vector();
    scalar_persistent_work_view_t Permuted_Xvector = gsHandler->get_permuted_x_vector();


    row_lno_persistent_work_view_t newxadj_ = gsHandler->get_new_xadj();
    row_lno_persistent_work_view_t old_to_new_map = gsHandler->get_old_to_new_map();
    row_lno_persistent_work_view_t newadj_ = gsHandler->get_new_adj();
    row_lno_persistent_work_view_t color_adj = gsHandler->get_color_adj();

    row_lno_t numColors = gsHandler->get_num_colors();



    if (update_y_vector){
      KokkosKernels::Experimental::Util::permute_vector
        <y_value_array_type,
        scalar_persistent_work_view_t,
        row_lno_persistent_work_view_t, MyExecSpace>(
          num_rows,
          old_to_new_map,
          y_rhs_input_vec,
          Permuted_Yvector
      );
    }
    MyExecSpace::fence();
    if(init_zero_x_vector){
      KokkosKernels::Experimental::Util::zero_vector<
          scalar_persistent_work_view_t, MyExecSpace>(num_cols, Permuted_Xvector);
    }
    else{
      KokkosKernels::Experimental::Util::permute_vector
        <x_value_array_type, scalar_persistent_work_view_t, row_lno_persistent_work_view_t, MyExecSpace>(
          num_cols,
          old_to_new_map,
          x_lhs_output_vec,
          Permuted_Xvector
          );
    }
    MyExecSpace::fence();

    row_lno_persistent_work_view_t permuted_xadj = gsHandler->get_new_xadj();
    row_lno_persistent_work_view_t permuted_adj = gsHandler->get_new_adj();
    scalar_persistent_work_view_t permuted_adj_vals = gsHandler->get_new_adj_val();
    scalar_persistent_work_view_t permuted_diagonals = gsHandler->get_permuted_diagonals();

    row_lno_persistent_work_host_view_t h_color_xadj = gsHandler->get_color_xadj();



    if (gsHandler->get_algorithm_type()== GS_PERMUTED){
      PSGS gs(permuted_xadj, permuted_adj, permuted_adj_vals,
          Permuted_Xvector, Permuted_Yvector, color_adj, permuted_diagonals);

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
          Permuted_Xvector, Permuted_Yvector,0,0, permuted_diagonals);

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
    <scalar_persistent_work_view_t,x_value_array_type,  row_lno_persistent_work_view_t, MyExecSpace>(
        num_cols,
        color_adj,
        Permuted_Xvector,
        x_lhs_output_vec
        );
    MyExecSpace::fence();

  }

  void IterativePSGS(
      Team_PSGS &gs,
      row_lno_t numColors,
      row_lno_persistent_work_host_view_t h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }

  void DoPSGS(Team_PSGS &gs, row_lno_t numColors, row_lno_persistent_work_host_view_t h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    int teamSizeMax = 0;
    int vector_size = 0;
    int max_allowed_team_size = team_policy_t::team_size_max(gs);


    row_lno_t nnz = this->entries.dimension_0();


    this->handle->get_gs_handle()->vector_team_size(max_allowed_team_size, vector_size, teamSizeMax, num_rows, nnz);
    /*std::cout
        << "max_allowed_team_size"  << max_allowed_team_size
        << " vector_size:" << vector_size
        << " teamSizeMax:" << teamSizeMax << std::endl;
    */
    if (apply_forward){
      for (row_lno_t i = 0; i < numColors; ++i){
        row_lno_t color_index_begin = h_color_xadj(i);
        row_lno_t color_index_end = h_color_xadj(i + 1);

        int overall_work = color_index_end - color_index_begin;// /256 + 1;


        gs._color_set_begin = color_index_begin;
        gs._color_set_end = color_index_end;

        Kokkos::parallel_for(
            team_policy_t(overall_work / teamSizeMax + 1 , teamSizeMax, vector_size),
            gs );
        MyExecSpace::fence();
      }
    }
    if (apply_backward){
      if (numColors > 0)
      for (row_lno_t i = numColors - 1;  ; --i){
        row_lno_t color_index_begin = h_color_xadj(i);
        row_lno_t color_index_end = h_color_xadj(i + 1);

        int numberOfTeams = color_index_end - color_index_begin;// /256 + 1;
        gs._color_set_begin = color_index_begin;
        gs._color_set_end = color_index_end;

        Kokkos::parallel_for(
            team_policy_t(numberOfTeams / teamSizeMax + 1 , teamSizeMax, vector_size),
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
      row_lno_t numColors,
      row_lno_persistent_work_host_view_t h_color_xadj,
      int num_iteration,
      bool apply_forward,
      bool apply_backward){

    for (int i = 0; i < num_iteration; ++i){
      this->DoPSGS(gs, numColors, h_color_xadj, apply_forward, apply_backward);
    }
  }



  void DoPSGS(PSGS &gs, row_lno_t numColors, row_lno_persistent_work_host_view_t h_color_xadj,
      bool apply_forward,
      bool apply_backward){
    if (apply_forward){
      for (row_lno_t i = 0; i < numColors; ++i){
        row_lno_t color_index_begin = h_color_xadj(i);
        row_lno_t color_index_end = h_color_xadj(i + 1);
        Kokkos::parallel_for (my_exec_space (color_index_begin, color_index_end) , gs);
        MyExecSpace::fence();
      }
    }
    if (apply_backward){
      for (row_lno_t i = numColors - 1; ; --i){
        row_lno_t color_index_begin = h_color_xadj(i);
        row_lno_t color_index_end = h_color_xadj(i + 1);
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
