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
#ifndef _KOKKOS_SPGEMM_HPP
#define _KOKKOS_SPGEMM_HPP

#include "KokkosSparse_spgemm_numeric.hpp"
#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_jacobi.hpp"

namespace KokkosSparse {

template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
void spgemm_symbolic(KernelHandle& kh, const AMatrix& A, const bool Amode,
                     const BMatrix& B, const bool Bmode, CMatrix& C) {
  using graph_type   = typename CMatrix::staticcrsgraph_type;
  using row_map_type = typename CMatrix::row_map_type::non_const_type;
  using entries_type = typename CMatrix::index_type::non_const_type;
  using values_type  = typename CMatrix::values_type::non_const_type;

  row_map_type row_mapC(
      Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"),
      A.numRows() + 1);
  entries_type entriesC;
  values_type valuesC;

  KokkosSparse::Experimental::spgemm_symbolic(
      &kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
      A.graph.entries, Amode, B.graph.row_map, B.graph.entries, Bmode,
      row_mapC);

  const size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size) {
    entriesC = entries_type(Kokkos::ViewAllocateWithoutInitializing("entriesC"),
                            c_nnz_size);
    valuesC  = values_type(Kokkos::ViewAllocateWithoutInitializing("valuesC"),
                          c_nnz_size);
  }

  graph_type graphC(entriesC, row_mapC);
  C = CMatrix("matrix", graphC);
}

template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
void spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode,
                    const BMatrix& B, const bool Bmode, CMatrix& C) {
  // using row_map_type = typename CMatrix::index_type::non_const_type;
  // using entries_type = typename CMatrix::row_map_type::non_const_type;
  // using values_type  = typename CMatrix::values_type::non_const_type;

  KokkosSparse::Experimental::spgemm_numeric(
      &kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
      A.graph.entries, A.values, Amode, B.graph.row_map, B.graph.entries,
      B.values, Bmode, C.graph.row_map, C.graph.entries, C.values);
  kh.destroy_spgemm_handle();
}

}  // namespace KokkosSparse

#endif
