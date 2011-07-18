// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef AMESOS2_SUPERLU_MATRIXHELPER_HPP
#define AMESOS2_SUPERLU_MATRIXHELPER_HPP

#include <Teuchos_ArrayRCP.hpp>

#include "Amesos2_MatrixHelper.hpp"
#include "Amesos2_Superlu_FunctionMap.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace SLU {
  extern "C" {
    typedef int int_t;
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"
  }
}

template <class Matrix, class Vector> class Superlu;


namespace Amesos2 {


  template <>
  struct MatrixHelper<Superlu>
  {

    /** \brief Creates a Superlu compressed-column Matrix from the given Matrix
     *
     * \tparam Matrix A matrix type conforming to the interface, in particular
     *         an Amesos2::MatrixAdapter<>.
     *
     * \param [in]     mat    The matrix which will be converted to SuperLU format
     * \param [in,out] nzval  A user-provided persisting store for the nonzero
     *                        values of the matrix.  The SuperLU matrix will expect
     *                        the array to persist throughout its existence.
     * \param [in,out] rowind A user-provided persisting store for the column
     *                        indices of the matrix.
     * \param [in,out] colptr User-provided persisting store for row pointers
     * \param [out]    A      Pointer to the SuperLU SuperMatrix which is to be constructed
     * \param [out]    mtxRedistTime Will have additional time added to it for the
     *                        time to redistribute the \c mat matrix.
     *
     * \callgraph
     */
    template <class Matrix>
    static void createCCSMatrix(const Teuchos::Ptr<Matrix>& mat,
                                const Teuchos::ArrayView<typename TypeMap<Superlu,typename Matrix::scalar_t>::type>& nzval,
                                const Teuchos::ArrayView<int>& rowind,
                                const Teuchos::ArrayView<int>& colptr,
                                const Teuchos::Ptr<SLU::SuperMatrix>& A,
                                Teuchos::Time& mtxRedistTime
                                )
    {
      typedef typename Matrix::scalar_t                         scalar_type;
      typedef typename Matrix::global_ordinal_t                     go_type;
      typedef typename Matrix::global_size_t                        gs_type;
      typedef typename TypeMap<Amesos2::Superlu,scalar_type>::type slu_type;
      // Get the SLU data type for this type of matrix
      SLU::Dtype_t dtype = TypeMap<Amesos2::Superlu,scalar_type>::dtype;

      // Extract the necessary information from mat and call SLU function
      using Teuchos::Array;
      using Teuchos::ArrayView;

      int nnz, rows, cols;
      nnz  = Teuchos::as<int>(mat->getGlobalNNZ());
      rows = Teuchos::as<int>(mat->getGlobalNumRows());
      cols = Teuchos::as<int>(mat->getGlobalNumCols());

      TEST_FOR_EXCEPTION( Teuchos::as<int>(nzval.size()) < nnz,
                          std::runtime_error,
                          "nzval array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<int>(rowind.size()) < nnz,
                          std::runtime_error,
                          "rowind array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<int>(colptr.size()) < rows + 1,
                          std::runtime_error,
                          "colptr array not large enough to hold data");

      int nnz_ret = 0;

      // Actually uses the compressed-row-store format, and tells Superlu this
      // while creating a compressed-column store.
      {
        Teuchos::TimeMonitor mtxRedistTimer( mtxRedistTime );

	Util::get_ccs_helper<Matrix,slu_type,int,int>::do_get(mat, nzval, rowind,
							      colptr, nnz_ret,
							      Util::Rooted,
							      Util::Arbitrary);
      }

      TEST_FOR_EXCEPTION( mat->getComm()->getRank() == 0 && nnz_ret != nnz,
			  std::runtime_error,
			  "Root rank failed to get all non-zero values in getCrs()");

      typedef FunctionMap<Superlu,scalar_type> function_map;
      function_map::create_CompCol_Matrix(A.getRawPtr(), rows, cols, nnz,
					  nzval.getRawPtr(), rowind.getRawPtr(),
					  colptr.getRawPtr(), SLU::SLU_NC,
					  dtype, SLU::SLU_GE);
    }


    /** \brief Creates a Superlu Dense Matrix from the given MultiVector
     *
     * \tparam MV A multi-vector type conforming to the interface, in particular
     *         an Amesos2::MultiVecAdapter<>.
     *
     * \param [in]     mv   The MultiVector which will be converted to SuperLU format
     * \param [in,out] vals A user-provided persisting store for the values of
     *                      the multi-vector.  The SuperLU matrix will expect the array
     *                      to persist throughout its existence.
     * \param [out]    ldx  Leading dimension of \c vals
     * \param [out]    X    Pointer to the SuperLU Dense SuperMatrix which is to be
     *                      constructed
     * \param [out]    vecRedistTime Will have time added for redistribution of \c \mv
     */
    template <class MV>
    static void createMVDenseMatrix(
        const Teuchos::Ptr<MV>& mv,
        const Teuchos::ArrayView<typename TypeMap<Superlu,typename MV::scalar_t>::type>& vals,
        size_t& ldx,
        const Teuchos::Ptr<SLU::SuperMatrix>& X,
        Teuchos::Time& vecRedistTime)
    {
      typedef typename MV::scalar_t scalar_type;
      typedef typename TypeMap<Superlu,scalar_type>::type slu_type;
      SLU::Dtype_t dtype = TypeMap<Superlu,scalar_type>::dtype;

      int rows, cols;
      rows = Teuchos::as<int>(mv->getGlobalLength());
      cols = Teuchos::as<int>(mv->getGlobalNumVectors());
      ldx  = Teuchos::as<size_t>(rows);

      {
        Teuchos::TimeMonitor redistTimer( vecRedistTime );

	Util::get_1d_copy_helper<MV,slu_type>::do_get(mv, vals, ldx, Util::Rooted);
    }

      FunctionMap<Superlu,scalar_type>::create_Dense_Matrix(
        X.getRawPtr(), rows, cols, vals.getRawPtr(), Teuchos::as<int>(ldx),
	SLU::SLU_DN, dtype, SLU::SLU_GE);
    }
    
  };				// end struct MatrixHelper

} // end namespace Amesos2

#endif  // end AMESOS2_SUPERLU_MATRIXHELPER_HPP
