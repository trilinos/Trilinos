// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "superlu_enum_consts.h"

void
at_plus_a(
          const int n,      /* number of columns in matrix A. */
          const int nz,     /* number of nonzeros in matrix A */
          int *colptr,      /* column pointer of size n+1 for matrix A. */
          int *rowind,      /* row indices of size nz for matrix A. */
          int *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
          int **b_colptr,   /* out - size n+1 */
          int **b_rowind    /* out - size *bnz */
          );


    namespace S {               // single-precision real definitions

      extern float slangs (char *, SLU::SuperMatrix *);

      extern void sgscon (char *, SuperMatrix *, SuperMatrix *,
                          float, float *, SuperLUStat_t*, int *);

      extern void
      sCompRow_to_CompCol(int, int, int, float*, int*, int*,
             float **, int **, int **);
      extern void
      sgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, float *, float *,
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      sgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
      extern void
      sgsisx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      sgsitrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
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

      extern double dlangs (char *, SLU::SuperMatrix *);

      extern void dgscon (char *, SuperMatrix *, SuperMatrix *,
                          double, double *, SuperLUStat_t*, int *);

      extern void
      dCompRow_to_CompCol(int, int, int, double*, int*, int*,
             double **, int **, int **);
      extern void
      dgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, double *, double *,
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      dgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
      extern void
      dgsisx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      dgsitrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
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

      extern float clangs (char *, SLU::SuperMatrix *);

      extern void cgscon (char *, SuperMatrix *, SuperMatrix *,
                          float, float *, SuperLUStat_t*, int *);

      extern void
      cCompRow_to_CompCol(int, int, int, complex*, int*, int*,
            complex **, int **, int **);
      extern void
      cgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, float *, float *,
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      cgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
      extern void
      cgsisx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, float *, float *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             float *, float *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      cgsitrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
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

      extern double zlangs (char *, SLU::SuperMatrix *);

      extern void zgscon (char *, SuperMatrix *, SuperMatrix *,
                          double, double *, SuperLUStat_t*, int *);

      extern void
      zCompRow_to_CompCol(int, int, int, doublecomplex*, int*, int*,
            doublecomplex **, int **, int **);
      extern void
      zgssvx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, double *, double *,
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      zgstrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
      extern void
      zgsisx(SLU::superlu_options_t *, SLU::SuperMatrix *, int *, int *, int *,
             char *, double *, double *, SLU::SuperMatrix *, SLU::SuperMatrix *,
             void *, int, SLU::SuperMatrix *, SLU::SuperMatrix *,
             double *, double *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
             SLU::GlobalLU_t*,
#endif
             SLU::mem_usage_t *, SLU::SuperLUStat_t *, int *);
      extern void
      zgsitrf (SLU::superlu_options_t*, SLU::SuperMatrix*,
              int, int, int*, void *, int, int *, int *,
              SLU::SuperMatrix *, SLU::SuperMatrix *, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
              SLU::GlobalLU_t*,
#endif
              SLU::SuperLUStat_t*, int *);
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

    static float langs(char *norm, SLU::SuperMatrix *A)
    {
      return SLU::S::slangs(norm, A);
    }

    static void gscon (char *norm, SLU::SuperMatrix *L, SLU::SuperMatrix *U,
                       float anorm, float *rcond, SLU::SuperLUStat_t *stat, int *info)
    {
      SLU::S::sgscon (norm, L, U, anorm, rcond, stat, info);
    }

    /**
     * \brief Binds to the appropriate Superlu solver driver based on data type
     */
    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gsisx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgsisx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
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
		      SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         stat, info);
    }

    static void gsitrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work,
		      int lwork, int* perm_c, int* perm_r, SLU::SuperMatrix* L,
		      SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::S::sgsitrf(options, AC, relax, panel_size, etree,
		      work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          lu, 
#endif
          stat, info);
    }

    /**
     * \brief Creates a Superlu CCS matrix using the appropriate function
     */
    template<class view_t>
    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
          Teuchos::Array<float> & convert_nzval, view_t & nzval,
          int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      // conversion not necessay - pass view data directly
      SLU::S::sCreate_CompCol_Matrix(A, m, n, nnz, nzval.data(), rowind, colptr,
				     stype, dtype, mtype);
    }

    /**
     * \brief Creates a Superlu CRS matrix using the appropriate function
     */
    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, 
          int nnz, float* nzval, int* rowind, int* colptr,
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
    template<class view_t>
    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
          Teuchos::Array<float> & convert_x, view_t & x,
          int ldx, SLU::Stype_t stype,
          SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      // conversion not necessay - pass view data directly
      SLU::S::sCreate_Dense_Matrix(X, m, n, x.data(), ldx, stype, dtype, mtype);
    }

    template<class view_t>
    static void convert_back_Dense_Matrix(
          Teuchos::Array<float> & convert_x, view_t & x)
    {
      // conversion not necessay - pass view data directly
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

    static double langs(char *norm, SLU::SuperMatrix *A)
    {
      return SLU::D::dlangs(norm, A);
    }

    static void gscon (char *norm, SLU::SuperMatrix *L, SLU::SuperMatrix *U,
                       double anorm, double *rcond, SLU::SuperLUStat_t *stat, int *info)
    {
      SLU::D::dgscon (norm, L, U, anorm, rcond, stat, info);
    }

    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         stat, info);
    }

    static void gsisx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgsisx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gsitrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		       int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		       int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
#ifdef HAVE_AMESOS2_SUPERLU5_API
		       SLU::GlobalLU_t* lu, 
#endif
           SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::D::dgsitrf(options, AC, relax, panel_size, etree,
		      work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          lu, 
#endif
          stat, info);
    }

    template<class view_t>
    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
          Teuchos::Array<double> & convert_nzval, view_t & nzval,
          int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      // conversion not necessay - pass view data directly
      SLU::D::dCreate_CompCol_Matrix(A, m, n, nnz, nzval.data(), rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, 
          int nnz, double* nzval, int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    template<class view_t>
    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
          Teuchos::Array<double> & convert_x, view_t & x,
          int ldx, SLU::Stype_t stype,
          SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      // conversion not necessay - pass view data directly
      SLU::D::dCreate_Dense_Matrix(X, m, n, x.data(), ldx, stype, dtype, mtype);
    }

    template<class view_t>
    static void convert_back_Dense_Matrix(
          Teuchos::Array<double> & convert_x, view_t & x)
    {
      // conversion not necessay - pass view data directly
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

  template <>
  struct FunctionMap<Superlu, Kokkos::complex<float>>
  {

    static float langs(char *norm, SLU::SuperMatrix *A)
    {
      return SLU::C::clangs(norm, A);
    }

    static void gscon (char *norm, SLU::SuperMatrix *L, SLU::SuperMatrix *U,
                       float anorm, float *rcond, SLU::SuperLUStat_t *stat, int *info)
    {
      SLU::C::cgscon (norm, L, U, anorm, rcond, stat, info);
    }

    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, float* ferr, float* berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         stat, info);
    }

    static void gsisx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, float* recip_pivot_growth,
		      float* rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgsisx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gsitrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		       int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		       int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
           SLU::GlobalLU_t* lu, 
#endif
		       SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::C::cgsitrf(options, AC, relax, panel_size, etree,
		      work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          lu, 
#endif
          stat, info);
    }

    template<class view_t>
    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
          Teuchos::Array<SLU::C::complex> & convert_nzval, view_t & nzval,
          int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      convert_nzval.resize(nnz);
      for(int i = 0; i < nnz; ++i) {
        convert_nzval[i] = Teuchos::as<SLU::C::complex>(nzval(i));
      }
      SLU::C::cCreate_CompCol_Matrix(A, m, n, nnz, convert_nzval.data(), rowind, colptr,
				     stype, dtype, mtype);
    }

    static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
          SLU::C::complex* nzval, int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
				     stype, dtype, mtype);
    }

    template<class view_t>
    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
          Teuchos::Array<SLU::C::complex> & convert_x, view_t & x,
          int ldx, SLU::Stype_t stype,
          SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      convert_x.resize(m * n);
      int write_index = 0;
      for(int j = 0; j < n; ++j) {
        for(int i = 0; i < m; ++i) { // layout left
          convert_x[write_index++] = Teuchos::as<SLU::C::complex>(x(i,j));
        }
      }
      SLU::C::cCreate_Dense_Matrix(X, m, n, convert_x.data(), ldx, stype, dtype, mtype);
    }

    template<class view_t>
    static void convert_back_Dense_Matrix(
          Teuchos::Array<SLU::C::complex> & convert_x, view_t & x)
    {
      int read_index = 0;
      for(int j = 0; j < static_cast<int>(x.extent(1)); ++j) {
        for(int i = 0; i < static_cast<int>(x.extent(0)); ++i) { // layout left
          x(i,j) = Teuchos::as<Kokkos::complex<float>>(convert_x[read_index++]);
        }
      }
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
  struct FunctionMap<Superlu,Kokkos::complex<double>>
  {

    static double langs(char *norm, SLU::SuperMatrix *A)
    {
      return SLU::Z::zlangs(norm, A);
    }

    static void gscon (char *norm, SLU::SuperMatrix *L, SLU::SuperMatrix *U,
                       double anorm, double *rcond, SLU::SuperLUStat_t *stat, int *info)
    {
      SLU::Z::zgscon (norm, L, U, anorm, rcond, stat, info);
    }

    static void gssvx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, double* ferr, double* berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, ferr, berr, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		      int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		      int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgstrf(options, AC, relax, panel_size, etree,
		     work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         stat, info);
    }

    static void gsisx(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
		      int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C,
		      SLU::SuperMatrix* L, SLU::SuperMatrix* U, void* work, int lwork,
		      SLU::SuperMatrix* B, SLU::SuperMatrix* X, double* recip_pivot_growth,
		      double* rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          SLU::GlobalLU_t* lu, 
#endif
          SLU::mem_usage_t* mem_usage,
		      SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgsisx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
		     lwork, B, X, recip_pivot_growth, rcond, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
         lu, 
#endif
         mem_usage, stat, info);
    }

    static void gsitrf(SLU::superlu_options_t* options, SLU::SuperMatrix* AC,
		       int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
		       int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
           SLU::GlobalLU_t* lu, 
#endif
		       SLU::SuperLUStat_t* stat, int* info)
    {
      SLU::Z::zgsitrf(options, AC, relax, panel_size, etree,
		      work, lwork, perm_c, perm_r, L, U, 
#ifdef HAVE_AMESOS2_SUPERLU5_API
          lu, 
#endif
          stat, info);
    }

    template<class view_t>
    static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
          Teuchos::Array<SLU::Z::doublecomplex> & convert_nzval, view_t & nzval,
          int* rowind, int* colptr,
          SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      convert_nzval.resize(nnz);
      for(int i = 0; i < nnz; ++i) {
        convert_nzval[i] = Teuchos::as<SLU::Z::doublecomplex>(nzval(i));
      }
      SLU::Z::zCreate_CompCol_Matrix(A, m, n, nnz, convert_nzval.data(), rowind, colptr,
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

    template<class view_t>
    static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
          Teuchos::Array<SLU::Z::doublecomplex> & convert_x, view_t & x,
          int ldx, SLU::Stype_t stype,
          SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      convert_x.resize(m * n);
      int write_index = 0;
      for(int j = 0; j < n; ++j) {
        for(int i = 0; i < m; ++i) { // layout left
          convert_x[write_index++] = Teuchos::as<SLU::Z::doublecomplex>(x(i,j));
        }
      }
      SLU::Z::zCreate_Dense_Matrix(X, m, n, convert_x.data(), ldx, stype, dtype, mtype);
    }

    template<class view_t>
    static void convert_back_Dense_Matrix(
          Teuchos::Array<SLU::Z::doublecomplex> & convert_x, view_t & x)
    {
      int read_index = 0;
      for(int j = 0; j < static_cast<int>(x.extent(1)); ++j) {
        for(int i = 0; i < static_cast<int>(x.extent(0)); ++i) { // layout left
          x(i,j) = Teuchos::as<Kokkos::complex<double>>(convert_x[read_index++]);
        }
      }
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
