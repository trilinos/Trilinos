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
#include <iostream>
#include <string>

#ifndef _SPTRSVHANDLE_HPP
#define _SPTRSVHANDLE_HPP

namespace KokkosSparse {
namespace Experimental {

// TP2 algorithm has issues with some offset-ordinal combo to be addressed
enum class SPTRSVAlgorithm { SEQLVLSCHD_RP, SEQLVLSCHD_TP1/*, SEQLVLSCHED_TP2*/ };

template <class size_type_, class lno_t_, class scalar_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPTRSVHandle {
public:

  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef ExecutionSpace execution_space;
  typedef HandlePersistentMemorySpace memory_space;


  typedef typename std::remove_const<size_type_>::type  size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type  nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type  nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;


  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> nnz_row_view_temp_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> nnz_row_view_t;
 // typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> nnz_scalar_view_temp_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> nnz_scalar_view_t;


  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_view_temp_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_view_t;
 // typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type


  typedef typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type signed_integral_t;
  typedef Kokkos::View< signed_integral_t*, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type, typename nnz_row_view_t::memory_traits > signed_nnz_lno_view_t;


private:

  signed_nnz_lno_view_t level_list;
  nnz_lno_view_t nodes_per_level;
  nnz_lno_view_t nodes_grouped_by_level;

  size_type nrows;
  size_type nlevel;

  bool lower_tri;

  bool symbolic_complete;

  SPTRSVAlgorithm algm;

  int team_size;
  int vector_size;

public:

  SPTRSVHandle ( SPTRSVAlgorithm choice, const size_type nrows_, bool lower_tri_, bool symbolic_complete_ = false ) :
    level_list(),
    nodes_per_level(),
    nodes_grouped_by_level(),
    nrows(nrows_),
    nlevel(0),
    lower_tri( lower_tri_ ),
    symbolic_complete( symbolic_complete_ ),
    algm(choice),
    team_size(-1),
    vector_size(-1)
  {}

#if 0
  SPTRSVHandle ( SPTRSVAlgorithm choice, const size_type nrows_, bool lower_tri_, bool symbolic_complete_ = false ) :
    level_list( Kokkos::ViewAllocateWithoutInitializing("level_list"), nrows),
    nodes_per_level("nodes_per_level", nrows),
    nodes_grouped_by_level("nodes_grouped_by_level", nrows),
    nrows(nrows_),
    nlevel(0),
    lower_tri( lower_tri_ ),
    symbolic_complete( symbolic_complete_ ),
    algm(choice)
  {
    // WithoutInitializing
    Kokkos::deep_copy( level_list, signed_integral_t(-1) );
  }

/*
  template <class rhslno_row_view_t_,
          class rhslno_nnz_view_t_,
          class rhsscalar_nnz_view_t_,
          class rhsExecutionSpace,
          class rhsMemorySpace>
  SPTRSVHandle ( SPTRSVHandle< rhslno_row_view_t_, rhslno_nnz_view_t_, rhsscalar_nnz_view_t_, rhsExecutionSpace, rhsMemorySpace > & rhs ) {

    this->level_list = rhs.level_list;
    this->nodes_per_level = rhs.nodes_per_level;
    this->nodes_grouped_by_level = rhs.nodes_grouped_by_level;
    this->nrows = rhs.nrows;
    this->nlevel = rhs.nlevel;
    this->lower_tri = rhs.lower_tri;
    this->symbolic_complete = rhs.symbolic_complete;
    this->algm = rhs.algm;
  }

  template <class rhslno_row_view_t_,
          class rhslno_nnz_view_t_,
          class rhsscalar_nnz_view_t_,
          class rhsExecutionSpace,
          class rhsMemorySpace>
  SPTRSVHandle & operator= ( SPTRSVHandle< rhslno_row_view_t_, rhslno_nnz_view_t_, rhsscalar_nnz_view_t_, rhsExecutionSpace, rhsMemorySpace > & rhs ) {

    this->level_list = rhs.level_list;
    this->nodes_per_level = rhs.nodes_per_level;
    this->nodes_grouped_by_level = rhs.nodes_grouped_by_level;
    this->nrows = rhs.nrows;
    this->nlevel = rhs.nlevel;
    this->lower_tri = rhs.lower_tri;
    this->symbolic_complete = rhs.symbolic_complete;
    this->algm = rhs.algm;
    return *this;
  }
*/

#endif

  void reset_handle( const size_type nrows_ ) {
    set_nrows(nrows_);
    set_num_levels(0);
    level_list = signed_nnz_lno_view_t( Kokkos::ViewAllocateWithoutInitializing("level_list"), nrows_),
    Kokkos::deep_copy( level_list, signed_integral_t(-1) );
    nodes_per_level =  nnz_lno_view_t("nodes_per_level", nrows_),
    nodes_grouped_by_level = nnz_lno_view_t("nodes_grouped_by_level", nrows_),
    reset_symbolic_complete();
  }

  virtual ~SPTRSVHandle() {};


  void set_algorithm(SPTRSVAlgorithm choice) { algm = choice; }

  SPTRSVAlgorithm get_algorithm() { return algm; }

  KOKKOS_INLINE_FUNCTION
  signed_nnz_lno_view_t get_level_list() const { return level_list; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_nodes_per_level() const { return nodes_per_level; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_nodes_grouped_by_level() const { return nodes_grouped_by_level; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nrows() const { return nrows; }

  KOKKOS_INLINE_FUNCTION
  void set_nrows(const size_type nrows_) { this->nrows = nrows_; }

  bool is_lower_tri() const { return lower_tri; }
  bool is_upper_tri() const { return !lower_tri; }

  bool is_symbolic_complete() const { return symbolic_complete; }

  size_type get_num_levels() const { return nlevel; }
  void set_num_levels(size_type nlevels_) { this->nlevel = nlevels_; }

  void set_symbolic_complete() { this->symbolic_complete = true; }
  void reset_symbolic_complete() { this->symbolic_complete = false; }

  void set_team_size(const int ts) {this->team_size = ts;}
  int get_team_size() const {return this->team_size;}

  void set_vector_size(const int vs) {this->vector_size = vs;}
  int get_vector_size() const {return this->vector_size;}

  void print_algorithm() { 
    if ( algm == SPTRSVAlgorithm::SEQLVLSCHD_RP )
      std::cout << "SEQLVLSCHD_RP" << std::endl;;

    if ( algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1 )
      std::cout << "SEQLVLSCHD_TP1" << std::endl;;

    /*
    if ( algm == SPTRSVAlgorithm::SEQLVLSCHED_TP2 ) {
      std::cout << "SEQLVLSCHED_TP2" << std::endl;;
      std::cout << "WARNING: With CUDA this is currently only reliable with int-int ordinal-offset pair" << std::endl;
    }
    */
  }

  inline SPTRSVAlgorithm StringToSPTRSVAlgorithm(std::string & name) {
    if(name=="SPTRSV_DEFAULT")             return SPTRSVAlgorithm::SEQLVLSCHD_RP;
    else if(name=="SPTRSV_RANGEPOLICY")    return SPTRSVAlgorithm::SEQLVLSCHD_RP;
    else if(name=="SPTRSV_TEAMPOLICY1")    return SPTRSVAlgorithm::SEQLVLSCHD_TP1;
    /*else if(name=="SPTRSV_TEAMPOLICY2")    return SPTRSVAlgorithm::SEQLVLSCHED_TP2;*/
    else
      throw std::runtime_error("Invalid SPTRSVAlgorithm name");
  }

};

} // namespace Experimental
} // namespace Kokkos

#endif
