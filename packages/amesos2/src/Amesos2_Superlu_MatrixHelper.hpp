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

namespace Amesos {


  template <>
  struct MatrixHelper<Superlu>
  {

    /** \brief Creates a Superlu compressed-row Matrix from the given Matrix
     *
     * \tparam Matrix A matrix type conforming to the interface, in particular
     *         an Amesos::MatrixAdapter<>.
     *
     * \param [in]     mat    The matrix which will be converted to SuperLU format
     * \param [in,out] nzval  A user-provided persisting store for the nonzero
     *                        values of the matrix.  The SuperLU matrix will expect
     *                        the array to persist throughout its existence.
     * \param [in,out] colind A user-provided persisting store for the column
     *                        indices of the matrix.
     * \param [in,out] rowptr User-provided persisting store for row pointers
     * \param [out]    A      Pointer to the SuperLU SuperMatrix which is to be constructed
     * \param [out]    mtxRedistTime Will have additional time added to it for the
     *                        time to redistribute the \c mat matrix.
     *
     * \callgraph
     */
    template <class Matrix>
    static void createCRSMatrix(const Teuchos::Ptr<Matrix>& mat,
                                const Teuchos::ArrayView<typename TypeMap<Superlu,typename Matrix::scalar_t>::type>& nzval,
                                const Teuchos::ArrayView<int>& colind,
                                const Teuchos::ArrayView<int>& rowptr,
                                const Teuchos::Ptr<SLU::SuperMatrix>& A,
                                Teuchos::Time& mtxRedistTime
                                )
    {
      typedef typename Matrix::scalar_t                        scalar_type;
      typedef typename Matrix::global_ordinal_t                    go_type;
      typedef typename Matrix::global_size_t                       gs_type;
      typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;
      // Get the SLU data type for this type of matrix
      SLU::Dtype_t dtype = TypeMap<Amesos::Superlu,scalar_type>::dtype;

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
      TEST_FOR_EXCEPTION( Teuchos::as<int>(colind.size()) < nnz,
                          std::runtime_error,
                          "colind array not large enough to hold data");
      TEST_FOR_EXCEPTION( Teuchos::as<int>(rowptr.size()) < rows + 1,
                          std::runtime_error,
                          "rowptr array not large enough to hold data");

      int nnz_ret = 0;

      // Actually uses the compressed-row-store format, and tells Superlu this
      // while creating a compressed-column store.
      {
        Teuchos::TimeMonitor mtxRedistTimer( mtxRedistTime );

	Util::get_crs_helper<Matrix,slu_type,int,int>::do_get(mat, nzval, colind,
							      rowptr, nnz_ret,
							      Util::Rooted,
							      Util::Arbitrary);
      }

      TEST_FOR_EXCEPTION( nnz_ret != nnz,
                          std::runtime_error,
                          "Number of nonzeros returned by getCrs() different from getGlobalNNZ()");

      FunctionMap<Superlu,scalar_type>::create_CompRow_Matrix(A.getRawPtr(), rows, cols, nnz, nzval.getRawPtr(),
                                                              colind.getRawPtr(), rowptr.getRawPtr(), SLU::SLU_NR,
                                                              dtype, SLU::SLU_GE);
    }


    /*
     * If the scalar and SuperLU types are the same, then we can do
     * a simple straight copy.
     */
    template <typename MV>
    struct same_type_get_copy {
      static void apply(const Teuchos::Ptr<MV>& mv,
			const Teuchos::ArrayView<typename MV::scalar_type>& v,
			const size_t& ldx)
      {
	bool get_root_global_copy = true;
	mv->get1dCopy(v, ldx, get_root_global_copy);
      }
    };

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding SuperLU type are different, then we need to
     * first get a copy of the scalar values, then convert each one
     * into the SuperLU type before inserting into the vals array.
     */
    template <typename MV>
    struct diff_type_get_copy {
      static void apply(const Teuchos::Ptr<MV>& mv,
			const Teuchos::ArrayView<typename TypeMap<Superlu,typename MV::scalar_type>::type>& v,
			const size_t& ldx)
      {
	typedef typename MV::scalar_type scalar_type;
	typedef typename TypeMap<Superlu,scalar_type>::type slu_type;
	  
	bool get_root_global_copy = true;
	int vals_length = v.size();
	Teuchos::Array<scalar_type> vals_tmp(vals_length);
	  
	mv->get1dCopy(vals_tmp(), ldx, get_root_global_copy);

	for ( int i = 0; i < vals_length; ++i ){
	  v[i] = Teuchos::as<slu_type>(vals_tmp[i]);
	}
      }
    };


    /** \brief Creates a Superlu Dense Matrix from the given MultiVector
     *
     * \tparam MV A multi-vector type conforming to the interface, in particular
     *         an Amesos::MultiVecAdapter<>.
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
        const Teuchos::ArrayView<typename TypeMap<Superlu,typename MV::scalar_type>::type>& vals,
        size_t& ldx,
        const Teuchos::Ptr<SLU::SuperMatrix>& X,
        Teuchos::Time& vecRedistTime)
    {
      typedef typename MV::scalar_type scalar_type;
      typedef typename TypeMap<Superlu,scalar_type>::type slu_type;
      SLU::Dtype_t dtype = TypeMap<Superlu,scalar_type>::dtype;

      int rows, cols;
      rows = Teuchos::as<int>(mv->getGlobalLength());
      cols = Teuchos::as<int>(mv->getGlobalNumVectors());
      ldx  = Teuchos::as<size_t>(rows);

      {
        Teuchos::TimeMonitor redistTimer( vecRedistTime );

	// Dispatch to the copy function appropriate for the type
        Util::if_then_else<Util::is_same<scalar_type,slu_type>::value,
          same_type_get_copy<MV>,
          diff_type_get_copy<MV> >::type::apply(mv, vals, ldx);
      }

      FunctionMap<Superlu,scalar_type>::create_Dense_Matrix(
        X.getRawPtr(), rows, cols, vals.getRawPtr(), Teuchos::as<int>(ldx),
	SLU::SLU_DN, dtype, SLU::SLU_GE);
    }


    /** \brief Creates a Superlu Dense Matrix from the given MultiVector
     *
     * \tparam MV A multi-vector type conforming to the interface, in
     *         particular an Amesos::MultiVecAdapter<>.
     *
     * \param [in]     mv   The MultiVector which will be converted to SuperLU format
     * \param [out]    X    Pointer to the SuperLU Dense SuperMatrix which is to be
     *                      constructed
     * \param [out]    vecRedistTime Will have time added for redistribution of \c \mv
     *
     * \return A Teuchos::ArrayRCP pointing to the beginning of a
     *         contiguous store of the values in \c X , which is
     *         <b>not</b> necessarily the beginning of the contiguous
     *         store of values in \c mv .
     *
     * \deprecated The other overloaded version of this function which
     * asks for a persisting store for the matrix values should be
     * favoured over this version.  This function uses the not-safe
     * `get1dViewNonConst'.  In any case, we eventually have to convert
     * the values into SuperLU's types, which makes a copy of the data
     * anyhow.
     */
    template <class MV>
    static
    Teuchos::ArrayRCP<typename TypeMap<Superlu,typename MV::scalar_type>::type>
    createMVDenseMatrix(
			const Teuchos::Ptr<MV>& mv,
			const Teuchos::Ptr<SLU::SuperMatrix>& X,
			Teuchos::Time& vecRedistTime
			)
    {
      typedef typename MV::scalar_type scalar_type;
      typedef typename TypeMap<Superlu,scalar_type>::type slu_type;
      SLU::Dtype_t dtype = TypeMap<Superlu,scalar_type>::dtype;

      int rows, cols, ldx;
      rows = Teuchos::as<int>(mv->getGlobalLength());
      cols = Teuchos::as<int>(mv->getGlobalNumVectors());
      // ldx  = Teuchos::as<int>(mv->getStride());
      ldx  = rows;

      Teuchos::ArrayRCP<scalar_type> vals_ptr;

      {
        Teuchos::TimeMonitor redistTimer( vecRedistTime );
        vals_ptr = mv->get1dViewNonConst();
      }
      typedef typename Teuchos::ArrayRCP<scalar_type>::size_type size_type;
      size_type vals_length = vals_ptr.size();

      typedef typename Teuchos::ArrayRCP<slu_type>::size_type slu_size_type;
      Teuchos::ArrayRCP<slu_type> slu_vals(Teuchos::as<slu_size_type>(vals_length));

      // Convert value types
      for ( size_type i = 0; i < vals_length; ++i ){
        slu_vals[i] = Teuchos::as<slu_type>(vals_ptr[i]);
      }

      FunctionMap<Superlu,scalar_type>::create_Dense_Matrix(
							    X.getRawPtr(), rows, cols, slu_vals.getRawPtr(), ldx,
							    SLU::SLU_DN, dtype, SLU::SLU_GE);

      return slu_vals;
    }
};                              // end struct MatrixHelper


} // end namespace Amesos

#endif  // end AMESOS2_SUPERLU_MATRIXHELPER_HPP
