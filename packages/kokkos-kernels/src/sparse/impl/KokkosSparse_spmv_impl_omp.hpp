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
namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOS_ENABLE_OPENMP
template<typename AMatrix, typename XVector, typename YVector>
void spmv_raw_openmp_no_transpose(typename YVector::const_value_type& s_a, AMatrix A, XVector x, typename YVector::const_value_type& s_b, YVector y) {

  typedef typename YVector::non_const_value_type value_type;
  typedef typename AMatrix::ordinal_type         ordinal_type;
  typedef typename AMatrix::non_const_size_type            size_type;

  typename XVector::const_value_type* KOKKOS_RESTRICT x_ptr = x.data();
  typename YVector::non_const_value_type* KOKKOS_RESTRICT y_ptr = y.data();

  const typename AMatrix::value_type* KOKKOS_RESTRICT matrixCoeffs = A.values.data();
  const ordinal_type* KOKKOS_RESTRICT matrixCols     = A.graph.entries.data();
  const size_type* KOKKOS_RESTRICT matrixRowOffsets  = A.graph.row_map.data();
  const size_type* KOKKOS_RESTRICT threadStarts     = A.graph.row_block_offsets.data();

#if defined(KOKKOS_ENABLE_PROFILING)
    uint64_t kpID = 0;
     if(Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::beginParallelFor("KokkosSparse::spmv<RawOpenMP,NoTranspose>", 0, &kpID);
     }
#endif

  typename YVector::const_value_type zero = 0;
  #pragma omp parallel
  {
#ifdef KOKKOS_COMPILER_INTEL
    __assume_aligned(x_ptr, 64);
    __assume_aligned(y_ptr, 64);
#endif

    const int myID    = omp_get_thread_num();
    const size_type myStart = threadStarts[myID];
    const size_type myEnd   = threadStarts[myID + 1];

    for(size_type row = myStart; row < myEnd; ++row) {
      const size_type rowStart = matrixRowOffsets[row];
      const size_type rowEnd   = matrixRowOffsets[row + 1];

      value_type sum = 0.0;

      for(size_type i = rowStart; i < rowEnd; ++i) {
        const ordinal_type x_entry =  matrixCols[i];
        const value_type alpha_MC  =  s_a * matrixCoeffs[i];
        sum                    += alpha_MC * x_ptr[x_entry];
      }

      if(zero == s_b) {
        y_ptr[row] = sum;
      } else {
        y_ptr[row] = s_b * y_ptr[row] + sum;
      }
   }
  }
#if defined(KOKKOS_ENABLE_PROFILING)
     if(Kokkos::Profiling::profileLibraryLoaded()) {
        Kokkos::Profiling::endParallelFor(kpID);
     }
#endif

}

#endif
}
}
