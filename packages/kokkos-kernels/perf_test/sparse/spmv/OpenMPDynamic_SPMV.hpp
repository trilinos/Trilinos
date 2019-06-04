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

#ifndef OPENMP_DYNAMIC_SPMV_HPP_
#define OPENMP_DYNAMIC_SPMV_HPP_

template<typename AType, typename XType, typename YType, typename LocalOrdinal, typename Scalar>
void openmp_dynamic_matvec(AType A, XType x, YType y, int rows_per_thread, int team_size, int vector_length) {

  #define OMP_BENCH_RESTRICT __restrict__

  const double s_a                                = 1.0;
  const double s_b                                = 0.0;

  const int rowCount                              = A.numRows();
  const double* OMP_BENCH_RESTRICT x_ptr          = (double*) x.data();
  double* OMP_BENCH_RESTRICT y_ptr                = (double*) y.data();
  const double* OMP_BENCH_RESTRICT matrixCoeffs   = A.values.data();
  const int* OMP_BENCH_RESTRICT matrixCols        = A.graph.entries.data();
  const int* OMP_BENCH_RESTRICT matrixRowOffsets  = &A.graph.row_map(0);
  
  #pragma omp parallel for schedule(dynamic)
  for(int row = 0; row < rowCount; ++row) {
  
  	const int rowStart = matrixRowOffsets[row];
  	const int rowEnd   = matrixRowOffsets[row + 1];
  	
  	double sum = 0.0;
  	
  	for(int i = rowStart; i < rowEnd; ++i) {
  		const int x_entry = matrixCols[i];
  		const double alpha_MC = s_a * matrixCoeffs[i];
  		sum += alpha_MC * x_ptr[x_entry];
  	}
  	
  	if(0.0 == s_b) {
  		y_ptr[row] = sum;
  	} else {
  		y_ptr[row] = s_b * y_ptr[row] + sum;
  	}
  
  }

  #undef OMP_BENCH_RESTRICT

}

#endif /* OPENMP_DYNAMIC_SPMV_HPP_ */
