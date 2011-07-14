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

#ifndef AMESOS2_SUPERLUDIST_MATRIXHELPER_HPP
#define AMESOS2_SUPERLUDIST_MATRIXHELPER_HPP

#include <Teuchos_ArrayRCP.hpp>

#include "Amesos2_MatrixHelper.hpp"
#include "Amesos2_Superludist_FunctionMap.hpp"
#include "Amesos2_Util.hpp"	// for get_ccs_helper class
#include "Amesos2_Util_is_same.hpp"

namespace SLUD {
  extern "C" {
    typedef int int_t;


  }
}

template <class Matrix, class Vector> class Superludist;

namespace Amesos {


  template <>
  struct MatrixHelper<Superludist>
  {

    /** \brief Creates a Superludist compressed-column Matrix from the given Matrix
     *
     * \tparam Matrix A matrix type conforming to the interface, in particular
     *         an Amesos::MatrixAdapter<>.
     *
     * \param [in]     mat    The matrix which will be converted to Superludist format
     * \param [in,out] nzval  A user-provided persisting store for the nonzero
     *                        values of the matrix.  The Superludist matrix will expect
     *                        the array to persist throughout its existence.
     * \param [in,out] rowind A user-provided persisting store for the column
     *                        indices of the matrix.
     * \param [in,out] colptr User-provided persisting store for row pointers
     * \param [out]    A      Pointer to the Superludist SuperMatrix which is to be constructed
     * \param [out]    mtxRedistTime Will have additional time added to it for the
     *                        time to redistribute the \c mat matrix.
     *
     * \callgraph
     */
    template <class Matrix>
    static void createCcsMatrix(const Teuchos::Ptr<const Matrix>& mat,
                                const Teuchos::ArrayView<typename TypeMap<Superludist,typename Matrix::scalar_t>::type>& nzval,
                                const Teuchos::ArrayView<int>& rowind,
                                const Teuchos::ArrayView<int>& colptr,
                                const Teuchos::Ptr<SLUD::SuperMatrix>& A,
                                Teuchos::Time& mtxRedistTime)
    {
      typedef typename Matrix::scalar_t                          scalar_type;
      typedef typename Matrix::global_ordinal_t                      go_type;
      typedef typename Matrix::global_size_t                         gs_type;
      typedef typename TypeMap<Superludist,scalar_type>::type slu_type;
      // Get the SLU data type for this type of matrix
      SLUD::Dtype_t dtype = TypeMap<Superludist,scalar_type>::dtype;

      // Extract the necessary information from mat and call SLU function
      using Teuchos::Array;
      using Teuchos::ArrayView;

      SLUD::int_t nnz, rows, cols;
      nnz  = Teuchos::as<SLUD::int_t>(mat->getGlobalNNZ());
      rows = Teuchos::as<SLUD::int_t>(mat->getGlobalNumRows());
      cols = Teuchos::as<SLUD::int_t>(mat->getGlobalNumCols());

#ifdef HAVE_AMESOS2_DEBUG
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(nzval.size()) < nnz,
                          std::runtime_error,
                          "nzval array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(rowind.size()) < nnz,
                          std::runtime_error,
                          "rowind array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(colptr.size()) < rows + 1,
                          std::runtime_error,
                          "colptr array not large enough to hold data");
#endif // HAVE_AMESOS2_DEBUG

      SLUD::int_t nnz_ret = 0;
      {
        Teuchos::TimeMonitor mtxRedistTimer( mtxRedistTime );

	Util::get_ccs_helper<
	  Matrix,
	  slu_type,
	  SLUD::int_t,
	  SLUD::int_t >::do_get(mat,
				nzval, rowind, colptr, nnz_ret,
				Util::Globally_Replicated,
				Util::Arbitrary);
      }

      TEST_FOR_EXCEPTION( nnz_ret != nnz,
                          std::runtime_error,
                          "Did not get the expected number of non-zero vals");

      typedef FunctionMap<Superludist,scalar_type> func_map;
      func_map::create_CompCol_Matrix(A.getRawPtr(),
				      rows, cols, nnz,
				      nzval.getRawPtr(),
				      rowind.getRawPtr(),
				      colptr.getRawPtr(),
				      SLUD::SLU_NC,
				      dtype, SLUD::SLU_GE);
    }

    /**
     * \brief Creates a SuperLU_DIST distributed Matrix from the
     * given Matrix.  The matrix is stored locally in compressed-row
     * format.
     *
     * \tparam Matrix A matrix type conforming to the interface, in particular
     *         an Amesos::MatrixAdapter<>.
     *
     * \param [in]     mat    The matrix which will be converted to SuperLU_DIST
     *                        NR_col format.  It is expected that the matrix has
     *                        a distribution matching that which SuperLU_DIST expects
     *                        for its symbolic and numeric factorization routines.
     *                        Namely, that only ranks < grid.nprow * grid.npcol have
     *                        rows of the matrix.
     * \param [in,out] nzval  A user-provided persisting store for the nonzero
     *                        values of the matrix.  The Superludist matrix will
     *                        expect the array to persist throughout its existence.
     * \param [in,out] rowind A user-provided persisting store for the column
     *                        indices of the matrix.
     * \param [in,out] colptr User-provided persisting store for row pointers
     * \param [out]    A      Pointer to the Superludist SuperMatrix which is to
     *                        be constructed
     * \param [out]    mtxRedistTime Will have additional time added to it for the
     *                        time to redistribute the \c mat matrix (if necessary).
     *
     * \note All processors should call this function.
     *
     * \callgraph
     */
    template <class Matrix>
    static void createNRlocMatrix(const Teuchos::Ptr<const Matrix>& mat,
				  const Teuchos::ArrayView<typename TypeMap<Superludist,typename Matrix::scalar_t>::type>& nzval,
				  const Teuchos::ArrayView<int>& colind,
				  const Teuchos::ArrayView<int>& rowptr,
				  const Teuchos::Ptr<SLUD::SuperMatrix>& A,
				  const Teuchos::Ptr<const Tpetra::Map<typename Matrix::local_ordinal_t, typename Matrix::global_ordinal_t, typename Matrix::node_t> > slu_dist_map,
				  Teuchos::Time& mtxRedistTime)
    {
      typedef typename Matrix::scalar_t                          scalar_type;
      typedef typename Matrix::global_ordinal_t                      go_type;
      typedef typename Matrix::global_size_t                         gs_type;
      typedef typename TypeMap<Superludist,scalar_type>::type slu_type;
      // Get the SLU data type for this type of matrix
      SLUD::Dtype_t dtype = TypeMap<Superludist,scalar_type>::dtype;

      // Extract the necessary information from mat and call SLU function
      using Teuchos::Array;
      using Teuchos::ArrayView;

#ifdef HAVE_AMESOS2_DEBUG
      slu_dist_map.assert_not_null();
#endif

      Teuchos::RCP<const Matrix> redist_mat = mat->get(slu_dist_map);

      SLUD::int_t l_nnz, l_rows, g_rows, g_cols, fst_global_row;
      l_nnz  = Teuchos::as<SLUD::int_t>(redist_mat->getLocalNNZ());
      l_rows = Teuchos::as<SLUD::int_t>(redist_mat->getLocalNumRows());
      g_rows = Teuchos::as<SLUD::int_t>(redist_mat->getGlobalNumRows());
      // g_cols = Teuchos::as<SLUD::int_t>(redist_mat->getGlobalNumCols());
      g_cols = g_rows;		// should be a square matrix anyhow
      fst_global_row = Teuchos::as<SLUD::int_t>(slu_dist_map->getMinGlobalIndex());

#ifdef HAVE_AMESOS2_DEBUG
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(nzval.size()) < l_nnz,
                          std::runtime_error,
                          "nzval array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(colind.size()) < l_nnz,
                          std::runtime_error,
                          "colind array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<SLUD::int_t>(rowptr.size()) < l_rows + 1,
                          std::runtime_error,
                          "rowptr array not large enough to hold data");
#endif

      SLUD::int_t nnz_ret = 0;
      {
        Teuchos::TimeMonitor mtxRedistTimer( mtxRedistTime );

	Util::get_crs_helper<
	  Matrix,
	  slu_type,
	  SLUD::int_t,
	  SLUD::int_t >::do_get(redist_mat.ptr(), nzval, colind,
				rowptr, nnz_ret, slu_dist_map,
				Util::Arbitrary);
      }

      TEST_FOR_EXCEPTION( nnz_ret != l_nnz,
                          std::runtime_error,
                          "Did not get the expected number of non-zero vals");

      typedef FunctionMap<Superludist,slu_type> func_map;
      func_map::create_CompRowLoc_Matrix(A.getRawPtr(),
					 g_rows, g_cols,
					 l_nnz, l_rows, fst_global_row,
					 nzval.getRawPtr(),
					 colind.getRawPtr(),
					 rowptr.getRawPtr(),
					 SLUD::SLU_NR_loc,
					 dtype, SLUD::SLU_GE);
    }

};                              // end struct MatrixHelper


} // end namespace Amesos

#endif  // end AMESOS2_SUPERLUDIST_MATRIXHELPER_HPP
