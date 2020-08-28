/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef _SPADDHANDLE_HPP
#define _SPADDHANDLE_HPP

namespace KokkosSparse{

template <class lno_row_view_t_,
          class lno_nnz_view_t_,
          class scalar_nnz_view_t_,
          class ExecutionSpace,
          class MemorySpace>
class SPADDHandle {
public:
  typedef typename lno_nnz_view_t_::non_const_type nnz_lno_view_t;
  typedef typename lno_row_view_t_::non_const_type nnz_row_view_t;
  typedef typename lno_row_view_t_::non_const_value_type size_type;
  typedef ExecutionSpace execution_space;
private:
  bool input_sorted;

  size_type result_nnz_size;

  bool called_symbolic;
  bool called_numeric;

  //a_pos and b_pos are used by the unsorted version of the kernel
  //both have same length as a_entries and b_entries
  //each entry provides the index in C row where the corresponding entry is added
  nnz_lno_view_t a_pos;
  nnz_lno_view_t b_pos;

public:
  /**
   * \brief sets the result nnz size.
   * \param result_nnz_size: size of the output matrix.
   */

  void set_a_b_pos(nnz_lno_view_t a_pos_in, const nnz_lno_view_t b_pos_in)
  {
    a_pos = a_pos_in;
    b_pos = b_pos_in;
  }

  nnz_lno_view_t get_a_pos()
  {
    return a_pos;
  }

  nnz_lno_view_t get_b_pos()
  {
    return b_pos;
  }

  /**
   * \brief sets the result nnz size.
   * \param result_nnz_size: size of the output matrix.
   */
  void set_c_nnz(size_type result_nnz_size_){
    this->result_nnz_size = result_nnz_size_;
  }

  /**
   * \brief returns the result nnz size.
   */
  size_type get_c_nnz(){
    return this->result_nnz_size;
  }

  void set_sort_option(int option){
    this->sort_option = option;
  }

  int get_sort_option(){
    return this->sort_option;
  }

  /**
   * \brief Default constructor.
   */
  SPADDHandle(bool input_is_sorted) :
    input_sorted(input_is_sorted), result_nnz_size(0),
    called_symbolic(false), called_numeric(false)
    {}

  virtual ~SPADDHandle() {};

  bool is_symbolic_called(){return this->called_symbolic;}
  bool is_numeric_called(){return this->called_numeric;}

  //setters
  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

  bool is_input_sorted()
  {
    return input_sorted;
  }
};

}

#endif
