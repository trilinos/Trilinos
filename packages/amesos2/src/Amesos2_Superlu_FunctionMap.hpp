// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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

/**
   \file   Amesos2_Superlu_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Mon May 31 23:38:46 2010

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_SUPERLU_FUNCTIONMAP_HPP
#define AMESOS2_SUPERLU_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_Superlu_TypeMap.hpp"


/* External definitions of the Superlu functions
 *
 * Note that we do include the "slu_*defs.h" files provided for each
 * data-type.  This produces linker warnings, but keeps us from
 * including SuperLU code in our own code (even if only extern
 * declarations, which would eliminate linker warnings).  This is
 * because there are several declarations (as of SuperLU 4.1) across
 * these headers which conflict with each other in C linkage.  All of
 * the conflicting functions, on the other hand, we do not care about.
 */
namespace SLU {

  extern "C" {
    typedef int int_t;
#include "supermatrix.h"
#include "slu_util.h"

    namespace S {               // single-precision real definitions
      extern void
      sgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, float *, float *,
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      sgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, SLU::SuperLUStat_t*, int *);
      extern void
      sCreate_CompCol_Matrix(SLU::SuperMatrix *, int, int, int, float *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      sCreate_CompRow_Matrix(SLU::SuperMatrix *, int, int, int, float *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      sCreate_Dense_Matrix(SLU::SuperMatrix *, int, int, float *, int,
                           SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);

      extern void
      sgsequ (SLU::SuperMatrix *, float *, float *, float *,
	      float *, float *, int *);

      extern void 
      slaqgs (SLU::SuperMatrix *, float *, float *, float,
              float, float, char *);

//#include "slu_sdefs.h"
    }

    namespace D {               // double-precision real definitions
      extern void
      dgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, double *, double *,
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      dgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, SLU::SuperLUStat_t*, int *);
      extern void
      dCreate_CompCol_Matrix(SLU::SuperMatrix *, int, int, int, double *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      dCreate_CompRow_Matrix(SLU::SuperMatrix *, int, int, int, double *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      dCreate_Dense_Matrix(SLU::SuperMatrix *, int, int, double *, int,
                           SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);

      extern void
      dlaqgs (SLU::SuperMatrix *, double *, double *, double,
              double, double, char *);

      extern void
      dgsequ (SLU::SuperMatrix *, double *, double *, double *,
	      double *, double *, int *);

//#include "slu_ddefs.h"
    }

#ifdef HAVE_TEUCHOS_COMPLEX
    namespace C {              // single-precision complex definitions
      extern void
      cgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, float *, float *,
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      cgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, SLU::SuperLUStat_t*, int *);
      extern void
      cCreate_CompCol_Matrix(SLU::SuperMatrix *, int, int, int, complex *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      cCreate_CompRow_Matrix(SLU::SuperMatrix *, int, int, int, complex *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      cCreate_Dense_Matrix(SLU::SuperMatrix *, int, int, complex *, int,
                           SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);

       extern void
       cgsequ (SLU::SuperMatrix *, float *, float *, float *,
               float *, float *, int *);

       extern void
       claqgs (SLU::SuperMatrix *, float *, float *, float,
               float, float, char *);

//#include "slu_cdefs.h"
    }

    namespace Z {              // double-precision complex definitions
      extern void
      zgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, double *, double *,
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      zgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, SLU::SuperLUStat_t*, int *);
      extern void
      zCreate_CompCol_Matrix(SLU::SuperMatrix *, int, int, int, doublecomplex *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      zCreate_CompRow_Matrix(SLU::SuperMatrix *, int, int, int, doublecomplex *,
                             int *, int *, SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);
      extern void
      zCreate_Dense_Matrix(SLU::SuperMatrix *, int, int, doublecomplex *, int,
                           SLU::Stype_t, SLU::Dtype_t, SLU::Mtype_t);

      extern void
      zgsequ (SLU::SuperMatrix *, double *, double *, double *,
              double *, double *, int *);

      extern void
      zlaqgs (SLU::SuperMatrix *, double *, double *, double,
              double, double, char *);

//#include "slu_zdefs.h"
    }
#endif  // HAVE_TEUCHOS_COMPLEX

  } // end extern "C"

} // end namespace SLU


namespace Amesos2 {

  /* ==================== Specializations ====================
   *
   * \cond Superlu_function_specializations
   */

  /**
   * \brief Pass function calls to Superlu based on data type.
   *
   * Helper class which passes on function calls to the appropriate
   * Superlu function based on the type of its scalar template argument.
   *
   * Superlu has solver and matrix builder functions defined based on
   * data type.  One function for complex, one for double precision
   * complex, another for \c float , and yet another for \c double.  To
   * work elegantly with the Amesos2::Superlu interface we want to be
   * able to perform a single function call which is appropriate for the
   * scalar type of the Matrix and MultiVectors that we are working
   * with.  The \c FunctionMap class provides that capability.
   *
   * The class template is specialized for each data type that Superlu
   * supports.  The Amesos2::create function assures that an
   * unspecialized FunctionMap will never be called by the solver
   * interface.
   *
   * Please see the <a
   * href="http://crd.lbl.gov/~xiaoye/SuperLU/superlu_ug.pdf">Superlu Users'
   * Guide</a> for more information on the TPL functions.
   */
  template <>
  struct FunctionMap<Superlu,float>
  {
    typedef TypeMap<Superlu,float> type_map;

    /**
     * \brief Binds to the appropriate Superlu solver driver based on data type
     */
    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
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
     * The AC argument is given in the SuperLU \c NCPformat
     *
     * See Superlu documentation for a further description of function
     * arguments.
     *
     * \note The SuperLU factorization methods only accept SuperMatrix objects
     * in the SLU_NC format, so conversion must be done when necessary
     */
    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work,
		      int lwork, int* perm_c, int* perm_r, SLU::SuperMatrix* L,
		      SLU::SuperMatrix* U, SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    /**
     * \brief Creates a Superlu CCS matrix using the appropriate function
     */
    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, 
				      int nnz, type_map::type* nzval, int* rowind, int* colptr, 
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    /**
     * \brief Creates a Superlu CRS matrix using the appropriate function
     */
    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, 
				      int nnz, type_map::type* nzval, int* rowind, int* colptr,
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }


    /**
     * \brief Creates a Superlu Dense Matrix using the appropriate Superlu
     *         function.
     *
     * \param X Superlu SuperMatrix that is to be created
     * \param x vals in column major order
     * \param ldx leading dimension of x
     */
    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLU::Stype_t stype,
				    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    /**
     * \brief compute row and column scaling for the matrix A
     */
    static void gsequ(SLU::SuperMatrix* A, float* R, float* C,
		      float* rowcnd, float* colcnd, float* amax, int* info)
    {
      SLU::S::sgsequ(A, R, C, rowcnd, colcnd, amax, info);
    }

    /**
     * \brief Apply row and column scaling to the matrix A
     *
     * The row and column scaling in R and C are applied to A if its
     * determined that such scalings would improce the condition of the
     * matrix.
     *
     * On exit, equed says what type of equilibration were actually
     * applied:
     *  - 'N' no equilibration
     *  - 'R' row equilibration
     *  - 'C' column equilibration
     *  - 'B' both row and column equilibration
     */
    static void laqgs(SLU::SuperMatrix* A, float* R, float* C,
		      float rowcnd, float colcnd, float amax, char* equed)
    {
      SLU::S::slaqgs(A, R, C, rowcnd, colcnd, amax, equed);
    }
  };


  template <>
  struct FunctionMap<Superlu,double>
  {
    typedef TypeMap<Superlu,double> type_map;

    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, 
				      int nnz, type_map::type* nzval, int* rowind, int* colptr, 
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, 
				      int nnz, type_map::type* nzval, int* rowind, int* colptr, 
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, 
				    int n, type_map::type* x, int ldx, SLU::Stype_t stype, 
				    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLU::SuperMatrix* A, double* R, double* C,
		      double* rowcnd, double* colcnd, double* amax, int* info)
    {
      SLU::D::dgsequ(A, R, C, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLU::SuperMatrix* A, double* R, double* C,
		      double rowcnd, double colcnd, double amax, char* equed)
    {
      SLU::D::dlaqgs(A, R, C, rowcnd, colcnd, amax, equed);
    }

  };


#ifdef HAVE_TEUCHOS_COMPLEX

  /* The specializations for Teuchos::as<> for SLU::complex and
   * SLU::doublecomplex are provided in Amesos2_Superlu_Type.hpp
   */
  template <>
  struct FunctionMap<Superlu,SLU::C::complex>
  {
    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
				      SLU::C::complex* nzval, int* rowind, int* colptr,
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
				      SLU::C::complex* nzval, int* rowind, int* colptr,
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
				    SLU::C::complex* x, int ldx, SLU::Stype_t stype,
				    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLU::SuperMatrix* A, float* R, float* C,
		      float* rowcnd, float* colcnd, float* amax, int* info)
    {
      SLU::C::cgsequ(A, R, C, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLU::SuperMatrix* A, float* R, float* C,
		      float rowcnd, float colcnd, float amax, char* equed)
    {
      SLU::C::claqgs(A, R, C, rowcnd, colcnd, amax, equed);
    }
  };


  template <>
  struct FunctionMap<Superlu,SLU::Z::doublecomplex>
  {
    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
				      SLU::Z::doublecomplex* nzval, int* rowind, int* colptr,
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);

      TEUCHOS_TEST_FOR_EXCEPTION( A == NULL,
			  std::runtime_error,
			  "Supermatrix A not initialized properly!");
    }


    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
				      SLU::Z::doublecomplex* nzval, int* rowind, int* colptr,
				      SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);

      TEUCHOS_TEST_FOR_EXCEPTION( A == NULL,
			  std::runtime_error,
			  "Supermatrix A not initialized properly!");
    }

    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
				    SLU::Z::doublecomplex* x, int ldx, SLU::Stype_t stype,
				    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void gsequ(SLU::SuperMatrix* A, double* R, double* C,
		      double* rowcnd, double* colcnd, double* amax, int* info)
    {
      SLU::Z::zgsequ(A, R, C, rowcnd, colcnd, amax, info);
    }

    static void laqgs(SLU::SuperMatrix* A, double* R, double* C,
		      double rowcnd, double colcnd, double amax, char* equed)
    {
      SLU::Z::zlaqgs(A, R, C, rowcnd, colcnd, amax, equed);
    }
  };
#endif	// HAVE_TEUCHOS_COMPLEX

  /* \endcond Superlu_function_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_FUNCTIONMAP_HPP
