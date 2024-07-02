// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_Superlumt_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Mon May 31 23:38:46 2010

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_SUPERLUMT_FUNCTIONMAP_HPP
#define AMESOS2_SUPERLUMT_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_Superlumt_TypeMap.hpp"


// External definitions of the SuperLU_MT functions
namespace SLUMT {
extern "C" {

typedef int int_t;

#include "slu_mt_util.h"
#include "pxgstrf_synch.h"	// preemptive inclusion

namespace S {
#include "pssp_defs.h"          // single-precision real definitions
}

namespace D {
#include "pdsp_defs.h"          // double-precision real definitions
}

#ifdef HAVE_TEUCHOS_COMPLEX
namespace C {
#include "pcsp_defs.h"          // single-precision complex definitions
}

namespace Z {
#include "pzsp_defs.h"          // double-precision complex definitions
}
#endif  // HAVE_TEUCHOS_COMPLEX

} // end extern "C"

} // end namespace SLUMT

namespace Amesos2 {

  template <class Matrix, class Vector> class Superlumt;

  /* ==================== Specializations ====================
   *
   * \cond SuperLU_MT_function_specializations
   */

  /*
   * Note that we don't need any generic declarations of the
   * SuperLU_MT functions that throw error in case the scalar type is
   * not supported.  This check is already performed in the factory
   * create method.  Just go straight for the specializations.
   */

  /**
   * \brief Pass function calls to SuperLU_MT based on data type.
   *
   * Helper class which passes on function calls to the appropriate
   * SuperLU_MT function based on the type of its scalar template
   * argument.
   *
   * SuperLU_MT has solver and matrix builder functions defined based on
   * data type.  One function for complex, one for double precision
   * complex, another for \c float , and yet another for \c double.  To
   * work elegantly with the Amesos2::SuperLU_MT interface we want to be
   * able to perform a single function call which is appropriate for the
   * scalar type of the Matrix and MultiVectors that we are working
   * with.  The \c FunctionMap class provides that capability.
   *
   * The class template is specialized for each data type that
   * SuperLU_MT supports, and errors are thrown for other data types.
   *
   * Please see the <a
   * href="http://crd.lbl.gov/~xiaoye/SuperLU/superlu_ug.pdf">Superlu Users'
   * Guide</a> for more information on the TPL functions.
   */
  template <>
  struct FunctionMap<Superlumt,float>
  {
    typedef TypeMap<Superlumt,float> type_map;

    /**
     * \brief Binds to the appropriate SuperLU_MT solver driver based on data type
     */
    static void gssvx(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, SLUMT::equed_t* equed, float* R, float* C,
		      SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U, void* work, int lwork,
		      SLUMT::SuperMatrix* B, SLUMT::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, SLUMT::superlu_memusage_t* mem_usage,
		      SLUMT::Gstat_t* stat, int* info)
    {
      options->etree = etree;
      options->perm_c = perm_c;
      options->perm_r = perm_r;

      options->work = work;
      options->lwork = lwork;
      
      SLUMT::S::psgssvx(options->nprocs, options, A, perm_c, perm_r,
			equed, R, C, L, U, B, X, recip_pivot_growth, rcond, ferr,
			berr, mem_usage, info);
    }

    /**
     * Solve the system A*X=B or A'*X=B using the L and U factors of A.
     */
    static void gstrs(SLUMT::trans_t trans, SLUMT::SuperMatrix* L,
		      SLUMT::SuperMatrix* U, int* perm_r, int* perm_c,
		      SLUMT::SuperMatrix* B, SLUMT::Gstat_t* Gstat, int* info)
    {
      SLUMT::S::sgstrs(trans, L, U, perm_r, perm_c, B, Gstat, info);
    }

    /**
     * \brief Computes an LU factorization of a general m-by-n matrix
     *
     * Uses a partial pivoting technique with row interchanges.  The
     * factorization has the form
     *
     * Pr * A = L * U
     *
     * where Pr is a row permutation matrix, L is lower triangular with unit
     * diagonal elements, and U is upper triangular.
     *
     * See Superlu documentation for a further description of function
     * arguments.
     *
     * \note The SuperLU factorization methods only accept SuperMatrix objects
     * in the SLUMT_NC format, so conversion must be done when necessary
     */
    static void gstrf(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_r, SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U,
		      SLUMT::Gstat_t* stat, int* info)
    {
      SLUMT::S::psgstrf(options, A, perm_r, L, U, stat, info);
    }

    /**
     * \brief Creates a SuperLU_MT CCS matrix using the appropriate function
     *
     * \internal Note: The following two functions are almost
     * identical to those found in FunctionMap<Superlu,D>.  We could
     * make use of these same functions, however, we would like to
     * keep the two packages as loosely coupled as possible, in order
     * to protect against unforseen incompatibilities arising as new
     * versions of each are released.
     * 
     * \throw std::runtime_error If there is no specialization of this type for
     *        the Scalar type
     */
    static void create_CompCol_Matrix(SLUMT::SuperMatrix* A, int m, int n, int nnz,
				      type_map::type* nzval, int* rowind, int* colptr,
				      SLUMT::Stype_t stype, SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::S::sCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				       stype, dtype, mtype);
    }

    /**
     * \brief Creates a SuperLU_MT Dense Matrix using the appropriate SuperLU_MT
     *         function.
     *
     * \param X SuperLU_MT SuperMatrix that is to be created
     * \param x vals in column major order
     * \param ldx leading dimension of x
     *
     * \throw std::runtime_error If there is no specialization of this type for
     *        the Scalar type
     */
    static void create_Dense_Matrix(SLUMT::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLUMT::Stype_t stype,
				    SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::S::sCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    /**
     * Equilibrates the matrix A.  Row scalings are placed in the array
     * \c r, and column scalings in the array \c c .  The estimated row
     * condition number is output in \c rowcnd , and the estimated
     * column condition number is output in \c colcnd .
     */
    static void gsequ(SLUMT::SuperMatrix* A,
		      type_map::magnitude_type* r,
		      type_map::magnitude_type* c,
		      type_map::magnitude_type* rowcnd,
		      type_map::magnitude_type* colcnd,
		      type_map::magnitude_type* amax,
		      int* info)
    {
      SLUMT::S::sgsequ(A, r, c, rowcnd, colcnd, amax, info);
    }

    /**
     * Apply equilibration to a matrix.  The parameters are expected to
     * be those as output from \c gsequ .  This function decides based
     * on \c rowcnd , \c colcnd , and \c amax whether it would be
     * worthwhile to apply the scalings, and outputs in \c equed what
     * type of equilibration was actually performed, whether \c ROW , \c
     * COL , or \c BOTH .
     */
    static void laqgs(SLUMT::SuperMatrix* A,
		      type_map::magnitude_type* r,
		      type_map::magnitude_type* c,
		      type_map::magnitude_type rowcnd,
		      type_map::magnitude_type colcnd,
		      type_map::magnitude_type amax,
		      SLUMT::equed_t* equed)
    {
      SLUMT::S::slaqgs(A, r, c, rowcnd, colcnd, amax, equed);
    }
  };


  template <>
  struct FunctionMap<Superlumt,double>
  {
    typedef TypeMap<Superlumt,double> type_map;

    static void gssvx(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, SLUMT::equed_t* equed, double* R, double* C,
		      SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U, void* work, int lwork,
		      SLUMT::SuperMatrix* B, SLUMT::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, SLUMT::superlu_memusage_t* mem_usage,
		      SLUMT::Gstat_t* stat, int* info)
    {
      options->etree = etree;
      options->perm_c = perm_c;
      options->perm_r = perm_r;

      options->work = work;
      options->lwork = lwork;
      
      SLUMT::D::pdgssvx(options->nprocs, options, A, perm_c, perm_r, 
			equed, R, C, L, U, B, X, recip_pivot_growth, rcond, ferr,
			berr, mem_usage, info);
    }

    static void gstrs(SLUMT::trans_t trans, SLUMT::SuperMatrix* L,
		      SLUMT::SuperMatrix* U, int* perm_r, int* perm_c,
		      SLUMT::SuperMatrix* B, SLUMT::Gstat_t* Gstat, int* info)
    {
      SLUMT::D::dgstrs(trans, L, U, perm_r, perm_c, B, Gstat, info);
    }

    static void gstrf(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_r, SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U,
		      SLUMT::Gstat_t* stat, int* info)
    {
      SLUMT::D::pdgstrf(options, A, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLUMT::SuperMatrix* A, int m, int n, int nnz,
				      type_map::type* nzval, int* rowind, int* colptr,
				      SLUMT::Stype_t stype, SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::D::dCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				       stype, dtype, mtype);
    }

    static void create_Dense_Matrix(SLUMT::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLUMT::Stype_t stype,
				    SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::D::dCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLUMT::SuperMatrix* A,
		      type_map::magnitude_type* r,
		      type_map::magnitude_type* c,
		      type_map::magnitude_type* rowcnd,
		      type_map::magnitude_type* colcnd,
		      type_map::magnitude_type* amax,
		      int* info)
    {
      SLUMT::D::dgsequ(A, r, c, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLUMT::SuperMatrix* A,
		      type_map::magnitude_type* r,
		      type_map::magnitude_type* c,
		      type_map::magnitude_type rowcnd,
		      type_map::magnitude_type colcnd,
		      type_map::magnitude_type amax,
		      SLUMT::equed_t* equed)
    {
      SLUMT::D::dlaqgs(A, r, c, rowcnd, colcnd, amax, equed);
    }
  };


#ifdef HAVE_TEUCHOS_COMPLEX
  /* The specializations for Teuchos::as<> for SLUMT::complex and
   * SLUMT::doublecomplex are provided in Amesos2_Superlumt_TypeMap.hpp
   */
  template <>
  struct FunctionMap<Superlumt,SLUMT::C::complex>
  {
    typedef TypeMap<Superlumt,SLUMT::C::complex> type_map;

    static void gssvx(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, SLUMT::equed_t* equed, float* R, float* C,
		      SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U, void* work, int lwork,
		      SLUMT::SuperMatrix* B, SLUMT::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, SLUMT::superlu_memusage_t* mem_usage,
		      SLUMT::Gstat_t* stat, int* info)
    {
      options->etree = etree;
      options->perm_c = perm_c;
      options->perm_r = perm_r;

      options->work = work;
      options->lwork = lwork;
      
      SLUMT::C::pcgssvx(options->nprocs, options, A, perm_c, perm_r, 
			equed, R, C, L, U, B, X, recip_pivot_growth, rcond, ferr,
			berr, mem_usage, info);
    }

    static void gstrs(SLUMT::trans_t trans, SLUMT::SuperMatrix* L,
		      SLUMT::SuperMatrix* U, int* perm_r, int* perm_c,
		      SLUMT::SuperMatrix* B, SLUMT::Gstat_t* Gstat, int* info)
    {
      SLUMT::C::cgstrs(trans, L, U, perm_r, perm_c, B, Gstat, info);
    }

    static void gstrf(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_r, SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U,
		      SLUMT::Gstat_t* stat, int* info)
    {
      SLUMT::C::pcgstrf(options, A, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLUMT::SuperMatrix* A, int m, int n, int nnz,
				      type_map::type* nzval, int* rowind, int* colptr,
				      SLUMT::Stype_t stype, SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::C::cCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				       stype, dtype, mtype);
    }

    static void create_Dense_Matrix(SLUMT::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLUMT::Stype_t stype,
				    SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::C::cCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLUMT::SuperMatrix* A, float* r, float* c,
		      float* rowcnd, float* colcnd, float* amax, int* info)
    {
      SLUMT::C::cgsequ(A, r, c, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLUMT::SuperMatrix* A, float* r, float* c, float rowcnd,
		      float colcnd, float amax, SLUMT::equed_t* equed)
    {
      SLUMT::C::claqgs(A, r, c, rowcnd, colcnd, amax, equed);
    }
  };


  template <>
  struct FunctionMap<Superlumt,SLUMT::Z::doublecomplex>
  {
    typedef TypeMap<Superlumt,SLUMT::Z::doublecomplex> type_map;

    static void gssvx(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, SLUMT::equed_t* equed, double* R, double* C,
		      SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U, void* work, int lwork,
		      SLUMT::SuperMatrix* B, SLUMT::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, SLUMT::superlu_memusage_t* mem_usage,
		      SLUMT::Gstat_t* stat, int* info)
    {
      options->etree = etree;
      options->perm_c = perm_c;
      options->perm_r = perm_r;

      options->work = work;
      options->lwork = lwork;
      
      SLUMT::Z::pzgssvx(options->nprocs, options, A, perm_c, perm_r, 
			equed, R, C, L, U, B, X, recip_pivot_growth, rcond, ferr,
			berr, mem_usage, info);
    }

    static void gstrs(SLUMT::trans_t trans, SLUMT::SuperMatrix* L,
		      SLUMT::SuperMatrix* U, int* perm_r, int* perm_c,
		      SLUMT::SuperMatrix* B, SLUMT::Gstat_t* Gstat, int* info)
    {
      SLUMT::Z::zgstrs(trans, L, U, perm_r, perm_c, B, Gstat, info);
    }

    static void gstrf(SLUMT::superlumt_options_t* options, SLUMT::SuperMatrix* A,
		      int* perm_r, SLUMT::SuperMatrix* L, SLUMT::SuperMatrix* U,
		      SLUMT::Gstat_t* stat, int* info)
    {
      SLUMT::Z::pzgstrf(options, A, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLUMT::SuperMatrix* A, int m, int n, int nnz,
				      type_map::type* nzval, int* rowind, int* colptr,
				      SLUMT::Stype_t stype, SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::Z::zCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				       stype, dtype, mtype);
    }

    static void create_Dense_Matrix(SLUMT::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLUMT::Stype_t stype,
				    SLUMT::Dtype_t dtype, SLUMT::Mtype_t mtype)
    {
      SLUMT::Z::zCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLUMT::SuperMatrix* A, double* r, double* c,
		      double* rowcnd, double* colcnd, double* amax, int* info)
    {
      SLUMT::Z::zgsequ(A, r, c, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLUMT::SuperMatrix* A, double* r, double* c, double rowcnd,
		      double colcnd, double amax, SLUMT::equed_t* equed)
    {
      SLUMT::Z::zlaqgs(A, r, c, rowcnd, colcnd, amax, equed);
    }
  };
#endif	// HAVE_TEUCHOS_COMPLEX

  /* \endcond SuperLU_MT_function_specializations */

} // end namespace Amesos2


#endif  // AMESOS2_SUPERLUMT_FUNCTIONMAP_HPP
