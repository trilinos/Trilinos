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

#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Utils.hpp>
#ifndef _GAUSSSEIDELHANDLE_HPP
#define _GAUSSSEIDELHANDLE_HPP
//#define VERBOSE

namespace KokkosSparse{

enum GSAlgorithm{GS_DEFAULT, GS_PERMUTED, GS_TEAM};

template <class size_type_, class lno_t_, class scalar_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GaussSeidelHandle{
public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type  size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type  nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type  nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;


  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename scalar_persistent_work_view_t::HostMirror scalar_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type


private:
  bool owner_of_coloring;
  GSAlgorithm algorithm_type;

  nnz_lno_persistent_work_host_view_t color_set_xadj;
  nnz_lno_persistent_work_view_t color_sets;
  nnz_lno_t numColors;

  row_lno_persistent_work_view_t permuted_xadj;
  nnz_lno_persistent_work_view_t permuted_adj;
  scalar_persistent_work_view_t permuted_adj_vals;
  nnz_lno_persistent_work_view_t old_to_new_map;

  bool called_symbolic;
  bool called_numeric;


  scalar_persistent_work_view_t permuted_y_vector;
  scalar_persistent_work_view_t permuted_x_vector;

  int suggested_vector_size;
  int suggested_team_size;

  scalar_persistent_work_view_t permuted_diagonals;
  nnz_lno_t block_size; //this is for block sgs

  nnz_lno_t max_nnz_input_row, num_values_in_l1, num_values_in_l2, num_big_rows;
  size_t level_1_mem, level_2_mem;
  public:

  /**
   * \brief Default constructor.
   */
  GaussSeidelHandle(GSAlgorithm gs = GS_DEFAULT):
    owner_of_coloring(false),
    algorithm_type(gs),
    color_set_xadj(), color_sets(), numColors(0),
    permuted_xadj(),  permuted_adj(), permuted_adj_vals(), old_to_new_map(),
    called_symbolic(false), called_numeric(false), permuted_y_vector(), permuted_x_vector(),
    suggested_vector_size(0), suggested_team_size(0), permuted_diagonals(), block_size(1), max_nnz_input_row(-1),
	num_values_in_l1(-1), num_values_in_l2(-1),num_big_rows(0), level_1_mem(0), level_2_mem(0)
    {
    if (gs == GS_DEFAULT){
      this->choose_default_algorithm();
    }


  }

  void set_block_size(nnz_lno_t bs){this->block_size = bs; }
  nnz_lno_t get_block_size(){return this->block_size;}

    /** \brief Chooses best algorithm based on the execution space. COLORING_EB if cuda, COLORING_VB otherwise.
   */
  void choose_default_algorithm(){
#if defined( KOKKOS_ENABLE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_THREADS )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->algorithm_type = GS_TEAM;
#ifdef VERBOSE
      std::cout << "Cuda Execution Space, Default Algorithm: GS_TEAM" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif
  }

  virtual ~GaussSeidelHandle(){};

  //getters
  GSAlgorithm get_algorithm_type() const {return this->algorithm_type;}
  bool is_owner_of_coloring() const {return this->owner_of_coloring;}

  nnz_lno_persistent_work_host_view_t get_color_xadj() {
    return this->color_set_xadj;
  }
  nnz_lno_persistent_work_view_t get_color_adj() {
    return this->color_sets;
  }
  nnz_lno_t get_num_colors() {
    return this->numColors;
  }

  row_lno_persistent_work_view_t get_new_xadj() {
    return this->permuted_xadj;
  }
  nnz_lno_persistent_work_view_t get_new_adj() {
    return this->permuted_adj;
  }
  scalar_persistent_work_view_t get_new_adj_val() {
    return this->permuted_adj_vals;
  }
  nnz_lno_persistent_work_view_t get_old_to_new_map() {
    return this->old_to_new_map;
  }

  bool is_symbolic_called(){return this->called_symbolic;}
  bool is_numeric_called(){return this->called_numeric;}

  //setters
  void set_algorithm_type(const GSAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_owner_of_coloring(bool owner = true){this->owner_of_coloring = owner;}

  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

  void set_color_set_xadj(const nnz_lno_persistent_work_host_view_t &color_set_xadj_) {
    this->color_set_xadj = color_set_xadj_;
  }
  void set_color_set_adj(const nnz_lno_persistent_work_view_t &color_sets_) {
    this->color_sets = color_sets_;
  }
  void set_num_colors(const nnz_lno_t &numColors_) {
    this->numColors = numColors_;
  }

  void set_new_xadj(const row_lno_persistent_work_view_t &xadj_) {
    this->permuted_xadj = xadj_;
  }
  void set_new_adj(const nnz_lno_persistent_work_view_t &adj_) {
    this->permuted_adj = adj_;
  }
  void set_new_adj_val(const scalar_persistent_work_view_t &adj_vals_) {
    this->permuted_adj_vals = adj_vals_;
  }
  void set_old_to_new_map(const nnz_lno_persistent_work_view_t &old_to_new_map_) {
    this->old_to_new_map = old_to_new_map_;
  }
  void set_permuted_diagonals (const scalar_persistent_work_view_t permuted_diagonals_){
    this->permuted_diagonals = permuted_diagonals_;
  }

  scalar_persistent_work_view_t get_permuted_diagonals (){
    return this->permuted_diagonals;
  }


  void set_level_1_mem(size_t _level_1_mem){
	  this->level_1_mem = _level_1_mem;
  }
  void set_level_2_mem(size_t _level_2_mem){
	  this->level_2_mem = _level_2_mem;
  }

  void set_num_values_in_l1(nnz_lno_t _num_values_in_l1){
	  this->num_values_in_l1 = _num_values_in_l1;
  }
  void set_num_values_in_l2(nnz_lno_t _num_values_in_l2){
	  this->num_values_in_l2 = _num_values_in_l2;
  }

  void set_num_big_rows(nnz_lno_t _big_rows){
	  this->num_big_rows = _big_rows;
  }

  size_t get_level_1_mem() const {
	  return this->level_1_mem;
  }
  size_t get_level_2_mem()const {
	  return this->level_2_mem;
  }

  nnz_lno_t get_num_values_in_l1()const {
	  return this->num_values_in_l1 ;
  }
  nnz_lno_t get_num_values_in_l2()const {
	  return this->num_values_in_l2 ;
  }
  nnz_lno_t get_num_big_rows()const {
	  return this->num_big_rows;
  }


  nnz_lno_t get_max_nnz() const{
    return this->max_nnz_input_row ;
  }


  void set_max_nnz(nnz_lno_t num_result_nnz_){
    this->max_nnz_input_row = num_result_nnz_;
  }


  void allocate_x_y_vectors(nnz_lno_t num_rows, nnz_lno_t num_cols){
    if(permuted_y_vector.extent(0) != size_t(num_rows)){
      permuted_y_vector = scalar_persistent_work_view_t("PERMUTED Y VECTOR", num_rows);
    }
    if(permuted_x_vector.extent(0) != size_t(num_cols)){
      permuted_x_vector = scalar_persistent_work_view_t("PERMUTED X VECTOR", num_cols);
    }
  }

  scalar_persistent_work_view_t get_permuted_y_vector (){return this->permuted_y_vector;}
  scalar_persistent_work_view_t get_permuted_x_vector (){return this->permuted_x_vector;}
  void vector_team_size(
      int max_allowed_team_size,
      int &suggested_vector_size_,
      int &suggested_team_size_,
      size_type nr, size_type nnz){
    //suggested_team_size_ =  this->suggested_team_size = 1;
    //suggested_vector_size_=this->suggested_vector_size = 1;
    //return;
    if (this->suggested_team_size && this->suggested_vector_size) {
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_ = this->suggested_team_size;
      return;
    }
    else {
      KokkosKernels::Impl::get_suggested_vector_team_size<size_type, ExecutionSpace>(
          max_allowed_team_size, suggested_vector_size_, suggested_team_size_, nr, nnz);
      this->suggested_team_size = suggested_vector_size_;
      this->suggested_vector_size = suggested_vector_size_;

    }
  }

};
}

#endif
