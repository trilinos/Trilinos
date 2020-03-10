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

#ifndef _SPILUKHANDLE_HPP
#define _SPILUKHANDLE_HPP

//#define EXPAND_FACT 3
#define KEEP_DIAG

namespace KokkosSparse {
namespace Experimental {

// TP2 algorithm has issues with some offset-ordinal combo to be addressed
enum class SPILUKAlgorithm { SEQLVLSCHD_RP, SEQLVLSCHD_TP1/*, SEQLVLSCHED_TP2*/ };

template <class size_type_, class lno_t_, class scalar_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPILUKHandle {
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

  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> nnz_row_view_t;
  
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_view_t;

  typedef typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type signed_integral_t;
  typedef Kokkos::View< signed_integral_t*, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type, typename nnz_row_view_t::memory_traits > signed_nnz_lno_view_t;


private:

  nnz_row_view_t level_list;//level IDs which the rows belong to
  nnz_lno_view_t level_idx; //the list of rows in each level
  nnz_lno_view_t level_ptr; //the starting index (into the view level_idx) of each level

  size_type nrows;
  size_type nlevel;
  size_type nnzL;
  size_type nnzU;
  size_type level_maxrows;//maximum number of rows of levels

  bool symbolic_complete;

  SPILUKAlgorithm algm;

  int team_size;
  int vector_size;

public:

  SPILUKHandle ( SPILUKAlgorithm choice, const size_type nrows_, const size_type nnzL_, const size_type nnzU_, bool symbolic_complete_ = false ) :
    level_list(),
    level_idx(),
    level_ptr(),
    nrows(nrows_),
    nlevel(0),
    nnzL(nnzL_),
    nnzU(nnzU_),
    level_maxrows(0),
    symbolic_complete( symbolic_complete_ ),
    algm(choice),
    team_size(-1),
    vector_size(-1)
  {}

  void reset_handle( const size_type nrows_, const size_type nnzL_, const size_type nnzU_ ) {
    set_nrows(nrows_);
    set_num_levels(0);
    set_nnzL(nnzL_);
    set_nnzU(nnzU_);
    set_level_maxrows(0);
    level_list = nnz_row_view_t("level_list", nrows_),
    level_idx  = nnz_lno_view_t("level_idx", nrows_),
    level_ptr  = nnz_lno_view_t("level_ptr", nrows_),
    reset_symbolic_complete();
  }

  virtual ~SPILUKHandle() {};


  void set_algorithm(SPILUKAlgorithm choice) { algm = choice; }

  SPILUKAlgorithm get_algorithm() { return algm; }

  KOKKOS_INLINE_FUNCTION
  nnz_row_view_t get_level_list() const { return level_list; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_level_idx() const { return level_idx; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_level_ptr() const { return level_ptr; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nrows() const { return nrows; }

  KOKKOS_INLINE_FUNCTION
  void set_nrows(const size_type nrows_) { this->nrows = nrows_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzL() const { return nnzL; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzL(const size_type nnzL_) { this->nnzL = nnzL_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzU() const { return nnzU; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzU(const size_type nnzU_) { this->nnzU = nnzU_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_level_maxrows() const { return level_maxrows; }

  KOKKOS_INLINE_FUNCTION
  void set_level_maxrows(const size_type level_maxrows_) { this->level_maxrows = level_maxrows_; }

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
    if ( algm == SPILUKAlgorithm::SEQLVLSCHD_RP )
      std::cout << "SEQLVLSCHD_RP" << std::endl;;

    if ( algm == SPILUKAlgorithm::SEQLVLSCHD_TP1 )
      std::cout << "SEQLVLSCHD_TP1" << std::endl;;

    /*
    if ( algm == SPILUKAlgorithm::SEQLVLSCHED_TP2 ) {
      std::cout << "SEQLVLSCHED_TP2" << std::endl;;
      std::cout << "WARNING: With CUDA this is currently only reliable with int-int ordinal-offset pair" << std::endl;
    }
    */
  }

  inline SPILUKAlgorithm StringToSPILUKAlgorithm(std::string & name) {
    if(name=="SPILUK_DEFAULT")             return SPILUKAlgorithm::SEQLVLSCHD_RP;
    else if(name=="SPILUK_RANGEPOLICY")    return SPILUKAlgorithm::SEQLVLSCHD_RP;
    else if(name=="SPILUK_TEAMPOLICY1")    return SPILUKAlgorithm::SEQLVLSCHD_TP1;
    /*else if(name=="SPILUK_TEAMPOLICY2")    return SPILUKAlgorithm::SEQLVLSCHED_TP2;*/
    else
      throw std::runtime_error("Invalid SPILUKAlgorithm name");
  }

};

} // namespace Experimental
} // namespace Kokkos

#endif
