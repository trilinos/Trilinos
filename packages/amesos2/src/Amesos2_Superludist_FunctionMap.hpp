// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_Superludist_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Tue Jun 21 13:37:55 MDT 2011

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_SUPERLUDIST_FUNCTIONMAP_HPP
#define AMESOS2_SUPERLUDIST_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_config.h"
#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_Superludist_TypeMap.hpp"


// Declarations of SuperLU_DIST types and namespace are found in
// Superludist_TypeMap.hpp

#define AMESOS2_SLUD_GET_DIAG_SCALE(eq) (((eq)=='N') ? SLUD::NOEQUIL : ((eq)=='R') ? SLUD::ROW : ((eq)=='C') ? SLUD::COL : SLUD::BOTH)

#define AMESOS2_SLUD_GET_EQUED(ds) (((ds)==SLUD::NOEQUIL) ? 'N' : ((ds)==SLUD::ROW) ? 'R' : ((ds)==SLUD::COL) ? 'C' : 'B')

namespace Amesos2 {

  template <class Matrix, class Vector> class Superludist;

  SLUD::DiagScale_t get_diag_scale(char eq);
  char get_equed(SLUD::DiagScale_t ds);


  /* ==================== Specializations ====================
   *
   * \cond SuperLU_DIST_function_specializations 
   */

  /**
   * \brief Pass function calls to SuperLU_DIST based on data type.
   *
   * Helper class which passes on function calls to the appropriate
   * SuperLU_DIST function based on the type of its scalar template
   * argument.
   *
   * SuperLU_DIST has solver and matrix builder functions defined based
   * on data type.  One function for double precision complex, and one
   * for \c double.  To work elegantly with the Amesos2::SuperLU_DIST
   * interface we want to be able to perform a single function call
   * which is appropriate for the scalar type of the Matrix and
   * MultiVectors that we are working with.  The \c FunctionMap class
   * provides that capability.
   *
   * The class template is specialized for each data type that
   * SuperLU_DIST supports.
   *
   * Please see the <a
   * href="http://crd.lbl.gov/~xiaoye/SuperLU/superlu_ug.pdf">Superlu Users'
   * Guide</a> for more information on the TPL functions.
   */
  template <>
  struct FunctionMap<Superludist,double>
  {
    typedef TypeMap<Superludist,double> type_map;
  
    /**
     * \brief Perform numeric factorization in parallel.
     * 
     * The factorization has the form
     *
     * Pr * A = L * U
     *
     * where Pr is a row permutation matrix, L is lower triangular with unit
     * diagonal elements, and U is upper triangular.
     *
     * See Superlu documentation for a further description of function
     * arguments.
     */
    static void gstrf(SLUD::amesos2_superlu_dist_options_t* options, int m, int n, double anorm, 
		      type_map::LUstruct_t* LU, SLUD::gridinfo_t* grid, SLUD::SuperLUStat_t* stat, 
		      int* info)
    {
      SLUD::D::pdgstrf(options, m, n, anorm, LU, grid, stat, info);
    }

    /**
     * \brief Solve the system A*X=B or A'*X=B using the L and U factors
     * of A.
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gstrs(SLUD::amesos2_superlu_dist_options_t* options,
                      SLUD::int_t n, type_map::LUstruct_t* lu_struct, 
		      SLUD::D::dScalePermstruct_t* scale_perm_struct, SLUD::gridinfo_t* grid,
		      type_map::type* B, SLUD::int_t l_numrows, SLUD::int_t fst_global_row, 
		      SLUD::int_t ldb, int nrhs, type_map::SOLVEstruct_t* solve_struct, 
		      SLUD::SuperLUStat_t* stat, int* info)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void gstrs(SLUD::int_t n, type_map::LUstruct_t* lu_struct, 
		      SLUD::D::dScalePermstruct_t* scale_perm_struct, SLUD::gridinfo_t* grid,
		      type_map::type* B, SLUD::int_t l_numrows, SLUD::int_t fst_global_row, 
		      SLUD::int_t ldb, int nrhs, type_map::SOLVEstruct_t* solve_struct, 
		      SLUD::SuperLUStat_t* stat, int* info)
#else
    static void gstrs(SLUD::int_t n, type_map::LUstruct_t* lu_struct, 
		      SLUD::ScalePermstruct_t* scale_perm_struct, SLUD::gridinfo_t* grid,
		      type_map::type* B, SLUD::int_t l_numrows, SLUD::int_t fst_global_row, 
		      SLUD::int_t ldb, int nrhs, type_map::SOLVEstruct_t* solve_struct, 
		      SLUD::SuperLUStat_t* stat, int* info)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::D::pdgstrs(options, n, lu_struct, scale_perm_struct, grid, B, l_numrows,
		       fst_global_row, ldb, nrhs, solve_struct, stat, info);
#else
      SLUD::D::pdgstrs(n, lu_struct, scale_perm_struct, grid, B, l_numrows,
		       fst_global_row, ldb, nrhs, solve_struct, stat, info);
#endif
    }

    /**
     * A variation of gstrs which requires B to be globally replicated.
     *
     * This function is available in case the user requests it.  Most
     * useful if the matrix is small enough to fit in memory on a single
     * processor.
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gstrs_Bglobal(SLUD::amesos2_superlu_dist_options_t* options,
                              SLUD::int_t n, type_map::LUstruct_t* lu_struct,
			      SLUD::gridinfo_t* grid, type_map::type* B,
			      SLUD::int_t ldb, int nrhs,
			      SLUD::SuperLUStat_t* stat, int* info)

    {
      SLUD::D::pdgstrs_Bglobal(options, n, lu_struct, grid, B, ldb, nrhs, stat, info);
    }
#else
    static void gstrs_Bglobal(SLUD::int_t n, type_map::LUstruct_t* lu_struct,
			      SLUD::gridinfo_t* grid, type_map::type* B,
			      SLUD::int_t ldb, int nrhs,
			      SLUD::SuperLUStat_t* stat, int* info)

    {
      SLUD::D::pdgstrs_Bglobal(n, lu_struct, grid, B, ldb, nrhs, stat, info);
    }
#endif

    /**
     * \brief Use iterative refined to improve the solution.
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gsrfs(SLUD::amesos2_superlu_dist_options_t* options,
                      SLUD::int_t n, SLUD::SuperMatrix* A, double anorm, 
		      type_map::LUstruct_t* lu_struct,
		      SLUD::D::dScalePermstruct_t* scale_perm, 
		      SLUD::gridinfo_t* grid, type_map::type* B, SLUD::int_t ldb, 
		      type_map::type* X, SLUD::int_t ldx, int nrhs, 
		      type_map::SOLVEstruct_t* solve_struct, double* berr, 
		      SLUD::SuperLUStat_t* stat, int* info)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void gsrfs(SLUD::int_t n, SLUD::SuperMatrix* A, double anorm, 
		      type_map::LUstruct_t* lu_struct,
		      SLUD::D::dScalePermstruct_t* scale_perm, 
		      SLUD::gridinfo_t* grid, type_map::type* B, SLUD::int_t ldb, 
		      type_map::type* X, SLUD::int_t ldx, int nrhs, 
		      type_map::SOLVEstruct_t* solve_struct, double* berr, 
		      SLUD::SuperLUStat_t* stat, int* info)
#else
    static void gsrfs(SLUD::int_t n, SLUD::SuperMatrix* A, double anorm, 
		      type_map::LUstruct_t* lu_struct,
		      SLUD::ScalePermstruct_t* scale_perm, 
		      SLUD::gridinfo_t* grid, type_map::type* B, SLUD::int_t ldb, 
		      type_map::type* X, SLUD::int_t ldx, int nrhs, 
		      type_map::SOLVEstruct_t* solve_struct, double* berr, 
		      SLUD::SuperLUStat_t* stat, int* info)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::D::pdgsrfs(options, n, A, anorm, lu_struct, scale_perm, grid, B, ldb, 
		       X, ldx, nrhs, solve_struct, berr, stat, info);
#else
      SLUD::D::pdgsrfs(n, A, anorm, lu_struct, scale_perm, grid, B, ldb, 
		       X, ldx, nrhs, solve_struct, berr, stat, info);
#endif
    }

    /**
     * \brief Use iterative refined to improve the solution.
     *
     * This variation requires A, X, and B to be globally replicated on
     * all calling processors.  This method is available in case the
     * user requests such functionality.
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gsrfs_ABXglobal(SLUD::amesos2_superlu_dist_options_t* options,
                                SLUD::int_t n, SLUD::SuperMatrix* A, double anorm,
				type_map::LUstruct_t* lu_struct, SLUD::gridinfo_t* grid,
				type_map::type* B, SLUD::int_t ldb, type_map::type* X,
				SLUD::int_t ldx, int nrhs, double* berr,
				SLUD::SuperLUStat_t* stat, int* info)
    {
      SLUD::D::pdgsrfs_ABXglobal(options, n, A, anorm, lu_struct, grid, B, ldb, 
				 X, ldx, nrhs, berr, stat, info);
    }
#else
    static void gsrfs_ABXglobal(SLUD::int_t n, SLUD::SuperMatrix* A, double anorm,
				type_map::LUstruct_t* lu_struct, SLUD::gridinfo_t* grid,
				type_map::type* B, SLUD::int_t ldb, type_map::type* X,
				SLUD::int_t ldx, int nrhs, double* berr,
				SLUD::SuperLUStat_t* stat, int* info)
    {
      SLUD::D::pdgsrfs_ABXglobal(n, A, anorm, lu_struct, grid, B, ldb, 
				 X, ldx, nrhs, berr, stat, info);
    }
#endif
  
    /**
     * \brief Creates a SuperLU_DIST distributed CRS matrix using the
     * appropriate function
     */
    static void create_CompRowLoc_Matrix(SLUD::SuperMatrix* A, SLUD::int_t g_numrows,
					 SLUD::int_t g_numcols, SLUD::int_t l_nnz,
					 SLUD::int_t l_numrows, SLUD::int_t fst_global_row,
					 type_map::type* nzval, SLUD::int_t* colind,
					 SLUD::int_t* rowptr, SLUD::Stype_t storage_t,
					 SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)
    {
      SLUD::D::dCreate_CompRowLoc_Matrix_dist(A, g_numrows, g_numcols, l_nnz, 
					      l_numrows, fst_global_row,
					      nzval, colind, rowptr,
					      storage_t, data_t, mat_t);
    }

    /**
     * \brief create a compressed-column matrix in the SuperLU_DIST style.
     *
     * This function should be used for creating a global CCS matrix for
     * factoring and solving with a globally-replicated matrix A and RHS
     * vector B.
     */
    static void create_CompCol_Matrix(SLUD::SuperMatrix* A, SLUD::int_t numrows,
				      SLUD::int_t numcols, SLUD::int_t nnz,
				      type_map::type* nzval, SLUD::int_t* rowind,
				      SLUD::int_t* colptr, SLUD::Stype_t storage_t,
				      SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)
    {
      SLUD::D::dCreate_CompCol_Matrix_dist(A, numrows, numcols, nnz, 
					   nzval, rowind, colptr,
					   storage_t, data_t, mat_t);
    }

    /**
     * \brief Creates a SuperLU_DIST Dense Matrix using the appropriate SuperLU_DIST
     *         function.
     *
     * \param X SuperLU_DIST SuperMatrix that is to be created
     * \param x vals in column major order
     * \param ldx leading dimension of x
     */
    static void create_Dense_Matrix(SLUD::SuperMatrix* X, int m, int n,
				    type_map::type* x, int ldx, SLUD::Stype_t stype,
				    SLUD::Dtype_t dtype, SLUD::Mtype_t mtype)
    {
      SLUD::D::dCreate_Dense_Matrix_dist(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void permute_Dense_Matrix(SLUD::int_t fst_row, SLUD::int_t m_loc,
				     SLUD::int_t* row_to_proc, SLUD::int_t* perm,
				     type_map::type* X, int ldx, type_map::type* B,
				     int ldb, int nrhs, SLUD::gridinfo_t* grid)
    {
      SLUD::D::pdPermute_Dense_Matrix(fst_row, m_loc, row_to_proc, perm,
				      X, ldx, B, ldb, nrhs, grid);
    }
  
    /**
     * Equilibrates the matrix A.  Row scalings are placed in the array
     * \c r, and column scalings in the array \c c .  The estimated row
     * condition number is output in \c rowcnd , and the estimated
     * column condition number is output in \c colcnd .
     *
     * This form operates on a SuperMatrix having the NRformat_loc
     */
    static void gsequ_loc(SLUD::SuperMatrix* A, double* r, double* c, 
			  double* rowcnd, double* colcnd, double* amax, SLUD::int_t* info,
			  SLUD::gridinfo_t* grid)
    {
      SLUD::D::pdgsequ(A, r, c, rowcnd, colcnd, amax, info, grid);
    }

    /**
     * This form of gsequ operates on matrices with format NCformat,
     * suitable for a globally-replicated matrix.
     */
    static void gsequ(SLUD::SuperMatrix* A, double* r, double* c, 
		      double* rowcnd, double* colcnd, double* amax, SLUD::int_t* info)
    {
      SLUD::D::dgsequ_dist(A, r, c, rowcnd, colcnd, amax, info);
    }

    /**
     * This form suitable for working with matrices in NC format.
     */
    static void laqgs_loc(SLUD::SuperMatrix* A, double* r, double* c, 
			  double rowcnd, double colcnd, double amax,
			  SLUD::DiagScale_t* equed)
    {
      char eq = AMESOS2_SLUD_GET_EQUED(*equed);
      SLUD::D::pdlaqgs(A, r, c, rowcnd, colcnd, amax, &eq);
      *equed = AMESOS2_SLUD_GET_DIAG_SCALE(eq);
    }

    /**
     * Apply equilibration to a matrix.  The parameters are expected to
     * be those as output from \c gsequ .  This function decides based
     * on \c rowcnd , \c colcnd , and \c amax whether it would be
     * worthwhile to apply the scalings, and outputs in \c equed what
     * type of equilibration was actually performed, whether \c ROW , \c
     * COL , or \c BOTH .
     *
     * This form suitable for working with A in NR_loc format
     * 
     * \internal our dispatcher should handle the conversion from char*
     * to DiagScale_t
     */
    static void laqgs(SLUD::SuperMatrix* A, double* r, double* c, 
		      double rowcnd, double colcnd, double amax, SLUD::DiagScale_t* equed)
    {
      char eq = AMESOS2_SLUD_GET_EQUED(*equed);
      SLUD::D::dlaqgs_dist(A, r, c, rowcnd, colcnd, amax, &eq);
      *equed = AMESOS2_SLUD_GET_DIAG_SCALE(eq);
    }

    /*
     * This version suitable for A in NCPformat
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void distribute(SLUD::amesos2_superlu_dist_options_t* options,
                           SLUD::fact_t fact, SLUD::int_t n,
			   SLUD::SuperMatrix* A, SLUD::Glu_freeable_t* glu_freeable,
			   type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
    {
      SLUD::D::ddistribute(options, n, A, glu_freeable, lu, grid);
    }
#else
    static void distribute(SLUD::fact_t fact, SLUD::int_t n,
			   SLUD::SuperMatrix* A, SLUD::Glu_freeable_t* glu_freeable,
			   type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
    {
      SLUD::D::ddistribute(fact, n, A, glu_freeable, lu, grid);
    }
#endif

    /*
     * This version suitable for A in NR_loc format.
     *
     * This routine should be used in the case where fact ==
     * SamePattern_SameRowPerm, otherwise dist_psymbtonum should be
     * called.o
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void pdistribute(SLUD::amesos2_superlu_dist_options_t* options,
                            SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::D::dScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void pdistribute(SLUD::fact_t fact, SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::D::dScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#else
    static void pdistribute(SLUD::fact_t fact, SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::ScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::D::pddistribute(options, n, A, scale_perm, glu_freeable, lu, grid);
#else
      SLUD::D::pddistribute(fact, n, A, scale_perm, glu_freeable, lu, grid);
#endif
    }

    /*
     * Distributes the input matrix A onto the 2D process grid based on
     * the L/U graph data in pslu_freeable.  On exit the struct lu
     * contains the information necessary to perform a numeric
     * factorization using gstrf.
     *
     * This routine should always be called with fact == DOFACT
     */
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void dist_psymbtonum(SLUD::amesos2_superlu_dist_options_t* options,
                                SLUD::int_t n, SLUD::SuperMatrix* A,
				SLUD::D::dScalePermstruct_t* scale_perm,
				SLUD::Pslu_freeable_t* pslu_freeable,
				type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void dist_psymbtonum(SLUD::fact_t fact, SLUD::int_t n, SLUD::SuperMatrix* A,
				SLUD::D::dScalePermstruct_t* scale_perm,
				SLUD::Pslu_freeable_t* pslu_freeable,
				type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
#else
    static void dist_psymbtonum(SLUD::fact_t fact, SLUD::int_t n, SLUD::SuperMatrix* A,
		          	SLUD::ScalePermstruct_t* scale_perm,
				SLUD::Pslu_freeable_t* pslu_freeable,
				type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::D::ddist_psymbtonum(options, n, A, scale_perm, pslu_freeable, lu, grid);
#else
      SLUD::D::ddist_psymbtonum(fact, n, A, scale_perm, pslu_freeable, lu, grid);
#endif
    }

    /*
     * The parameter norm may be one of:
     *  - 'M' for the max absolute matrix entry value
     *  - '1' for the norm1(A)
     *  - 'I' for the normI(A)
     *  - 'F' for the Frobenius norm of A
     */
    static double plangs(char* norm, SLUD::SuperMatrix* A, SLUD::gridinfo_t* grid)
    {
      return SLUD::D::pdlangs(norm, A, grid);
    }

    static void SolveInit(SLUD::amesos2_superlu_dist_options_t* options, SLUD::SuperMatrix* A, 
			  SLUD::int_t* perm_r, SLUD::int_t* perm_c, SLUD::int_t nrhs, 
			  type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid, 
			  type_map::SOLVEstruct_t* solve_struct)
    {
      SLUD::D::dSolveInit(options, A, perm_r, perm_c, nrhs, lu, grid, solve_struct);
    }

    static void LUstructInit(SLUD::int_t m, SLUD::int_t n,
			     type_map::LUstruct_t* lu)
    {
      /// When we make sure that version 5 and higher is used
      /// we do not perform runtime check of the interface
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::D::dLUstructInit(n, lu);
#elif defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
      SLUD::D::LUstructInit(n, lu);
#else      
#ifdef HAVE_SUPERLUDIST_LUSTRUCTINIT_2ARG
      SLUD::D::LUstructInit(n, lu);
#else
      SLUD::D::LUstructInit(m, n, lu);
#endif
#endif
    }

    static void Destroy_LU(SLUD::int_t m, SLUD::gridinfo_t* grid,
			   type_map::LUstruct_t* lu)
    {
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::D::dDestroy_LU(m, grid, lu);
#else
      SLUD::D::Destroy_LU(m, grid, lu);
#endif
    }

    static void LUstructFree(type_map::LUstruct_t* lu)
    {
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::D::dLUstructFree(lu);
#else
      SLUD::D::LUstructFree(lu);
#endif
    }

    static void SolveFinalize(SLUD::amesos2_superlu_dist_options_t* options,
			      type_map::SOLVEstruct_t* solve_struct)
    {
      SLUD::D::dSolveFinalize(options, solve_struct);
    }

    static int ldperm_dist(int job, int n, SLUD::int_t nnz, SLUD::int_t colptr[], SLUD::int_t adjncy[],
            double nzval[], SLUD::int_t *perm, double u[], double v[])
    {
      return SLUD::D::dldperm_dist(job, n, nnz, colptr, adjncy, nzval, perm, u, v);
    }
  };


#if defined(HAVE_TEUCHOS_COMPLEX)  && !defined(__clang__)
  /* The specializations for Teuchos::as<> for SLUD::complex and
   * SLUD::doublecomplex are provided in Amesos2_Superlu_Type.hpp
   */
  template <>
  struct FunctionMap<Superludist,SLUD::Z::doublecomplex>
  {
    typedef TypeMap<Superludist,std::complex<double> > type_map;

    static void gstrf(SLUD::amesos2_superlu_dist_options_t* options, int m, int n, double anorm, 
		      type_map::LUstruct_t* LU, SLUD::gridinfo_t* grid,
		      SLUD::SuperLUStat_t* stat, int* info)
    {
      SLUD::Z::pzgstrf(options, m, n, anorm, LU, grid, stat, info);
    }

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gstrs(SLUD::amesos2_superlu_dist_options_t* options,
                      SLUD::int_t n, type_map::LUstruct_t* lu_struct,
		      SLUD::Z::zScalePermstruct_t* scale_perm_struct,
		      SLUD::gridinfo_t* grid, type_map::type* B,
		      SLUD::int_t l_numrows, SLUD::int_t fst_global_row,
		      SLUD::int_t ldb, int nrhs,
		      type_map::SOLVEstruct_t* solve_struct,
		      SLUD::SuperLUStat_t* stat, int* info)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void gstrs(SLUD::int_t n, type_map::LUstruct_t* lu_struct,
		      SLUD::Z::zScalePermstruct_t* scale_perm_struct,
		      SLUD::gridinfo_t* grid, type_map::type* B,
		      SLUD::int_t l_numrows, SLUD::int_t fst_global_row,
		      SLUD::int_t ldb, int nrhs,
		      type_map::SOLVEstruct_t* solve_struct,
		      SLUD::SuperLUStat_t* stat, int* info)
#else
    static void gstrs(SLUD::int_t n, type_map::LUstruct_t* lu_struct,
		      SLUD::ScalePermstruct_t* scale_perm_struct,
		      SLUD::gridinfo_t* grid, type_map::type* B,
		      SLUD::int_t l_numrows, SLUD::int_t fst_global_row,
		      SLUD::int_t ldb, int nrhs,
		      type_map::SOLVEstruct_t* solve_struct,
		      SLUD::SuperLUStat_t* stat, int* info)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::Z::pzgstrs(options, n, lu_struct, scale_perm_struct, grid, B, l_numrows,
		       fst_global_row, ldb, nrhs, solve_struct, stat, info);
#else
      SLUD::Z::pzgstrs(n, lu_struct, scale_perm_struct, grid, B, l_numrows,
		       fst_global_row, ldb, nrhs, solve_struct, stat, info);
#endif
    }

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void gstrs_Bglobal(SLUD::amesos2_superlu_dist_options_t* options,
                              SLUD::int_t n, type_map::LUstruct_t* lu_struct, 
			      SLUD::gridinfo_t* grid, type_map::type* B, 
			      SLUD::int_t ldb, int nrhs, SLUD::SuperLUStat_t* stat, int* info)
    {
      SLUD::Z::pzgstrs_Bglobal(options,n, lu_struct, grid, B, ldb, nrhs, stat, info);
    }
#else
    static void gstrs_Bglobal(SLUD::int_t n, type_map::LUstruct_t* lu_struct, 
			      SLUD::gridinfo_t* grid, type_map::type* B, 
			      SLUD::int_t ldb, int nrhs, SLUD::SuperLUStat_t* stat, int* info)
    {
      SLUD::Z::pzgstrs_Bglobal(n, lu_struct, grid, B, ldb, nrhs, stat, info);
    }
#endif
  
    static void create_CompRowLoc_Matrix(SLUD::SuperMatrix* A, SLUD::int_t g_numrows,
					 SLUD::int_t g_numcols, SLUD::int_t l_nnz,
					 SLUD::int_t l_numrows, SLUD::int_t fst_global_row,
					 type_map::type* nzval, SLUD::int_t* colind,
					 SLUD::int_t* rowptr, SLUD::Stype_t storage_t,
					 SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)
    {
      SLUD::Z::zCreate_CompRowLoc_Matrix_dist(A, g_numrows, g_numcols, l_nnz, 
					      l_numrows, fst_global_row,
					      nzval, colind, rowptr,
					      storage_t, data_t, mat_t);
    }

    static void create_CompCol_Matrix(SLUD::SuperMatrix* A, SLUD::int_t numrows,
				      SLUD::int_t numcols, SLUD::int_t nnz,
				      type_map::type* nzval, SLUD::int_t* rowind,
				      SLUD::int_t* colptr, SLUD::Stype_t storage_t,
				      SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)
    {
      SLUD::Z::zCreate_CompCol_Matrix_dist(A, numrows, numcols, nnz, 
					   nzval, rowind, colptr,
					   storage_t, data_t, mat_t);
    }

    static void create_Dense_Matrix(SLUD::SuperMatrix* X, int m, int n,
				    TypeMap<Superludist,std::complex<double> >::type* x, int ldx, 
				    SLUD::Stype_t stype, SLUD::Dtype_t dtype, SLUD::Mtype_t mtype)
    {
      SLUD::Z::zCreate_Dense_Matrix_dist(X, m, n, x, ldx, stype, dtype, mtype);
    }

    static void permute_Dense_Matrix(SLUD::int_t fst_row, SLUD::int_t m_loc,
				     SLUD::int_t* row_to_proc, SLUD::int_t* perm,
				     type_map::type* X, int ldx,
				     type_map::type* B, int ldb,
				     int nrhs, SLUD::gridinfo_t* grid)
    {
      SLUD::Z::pzPermute_Dense_Matrix(fst_row, m_loc, row_to_proc, perm,
				      X, ldx, B, ldb, nrhs, grid);
    }
  
    static void gsequ_loc(SLUD::SuperMatrix* A, double* r, double* c, 
			  double* rowcnd, double* colcnd, double* amax, SLUD::int_t* info, 
			  SLUD::gridinfo_t* grid)
    {
      SLUD::Z::pzgsequ(A, r, c, rowcnd, colcnd, amax, info, grid);
    }

    static void gsequ(SLUD::SuperMatrix* A, double* r, double* c, 
		      double* rowcnd, double* colcnd, double* amax, SLUD::int_t* info)
    {
      SLUD::Z::zgsequ_dist(A, r, c, rowcnd, colcnd, amax, info);
    }

    static void laqgs_loc(SLUD::SuperMatrix* A, double* r, double* c, 
			  double rowcnd, double colcnd, double amax, SLUD::DiagScale_t* equed)
    {
      char eq = AMESOS2_SLUD_GET_EQUED(*equed);
      SLUD::Z::pzlaqgs(A, r, c, rowcnd, colcnd, amax, &eq);
      *equed = AMESOS2_SLUD_GET_DIAG_SCALE(eq);
    }

    static void laqgs(SLUD::SuperMatrix* A, double* r, double* c, 
		      double rowcnd, double colcnd, double amax, SLUD::DiagScale_t* equed)
    {
      char eq = AMESOS2_SLUD_GET_EQUED(*equed);
      SLUD::Z::zlaqgs_dist(A, r, c, rowcnd, colcnd, amax, &eq);
      *equed = AMESOS2_SLUD_GET_DIAG_SCALE(eq);
    }

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void distribute(SLUD::amesos2_superlu_dist_options_t* options,
                           SLUD::int_t n,
			   SLUD::SuperMatrix* A, SLUD::Glu_freeable_t* glu_freeable,
			   type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
    {
      SLUD::Z::zdistribute(options, n, A, glu_freeable, lu, grid);
    }
#else
    static void distribute(SLUD::fact_t fact, SLUD::int_t n,
			   SLUD::SuperMatrix* A, SLUD::Glu_freeable_t* glu_freeable,
			   type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid)
    {
      SLUD::Z::zdistribute(fact, n, A, glu_freeable, lu, grid);
    }
#endif

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void pdistribute(SLUD::amesos2_superlu_dist_options_t* options,
                            SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::Z::zScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void pdistribute(SLUD::fact_t fact, SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::Z::zScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#else
    static void pdistribute(SLUD::fact_t fact, SLUD::int_t n, 
			    SLUD::SuperMatrix* A, SLUD::ScalePermstruct_t* scale_perm, 
			    SLUD::Glu_freeable_t* glu_freeable, type_map::LUstruct_t* lu,
			    SLUD::gridinfo_t* grid)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::Z::pzdistribute(options, n, A, scale_perm, glu_freeable, lu, grid);
#else
      SLUD::Z::pzdistribute(fact, n, A, scale_perm, glu_freeable, lu, grid);
#endif
    }

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
    static void dist_psymbtonum(SLUD::amesos2_superlu_dist_options_t* options,
                                SLUD::int_t n,
				SLUD::SuperMatrix* A, SLUD::Z::zScalePermstruct_t* scale_perm, 
				SLUD::Pslu_freeable_t* pslu_freeable, type_map::LUstruct_t* lu,
				SLUD::gridinfo_t* grid)
#elif SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
    static void dist_psymbtonum(SLUD::fact_t fact, SLUD::int_t n,
				SLUD::SuperMatrix* A, SLUD::Z::zScalePermstruct_t* scale_perm, 
				SLUD::Pslu_freeable_t* pslu_freeable, type_map::LUstruct_t* lu,
				SLUD::gridinfo_t* grid)
#else
    static void dist_psymbtonum(SLUD::fact_t fact, SLUD::int_t n,
		  	        SLUD::SuperMatrix* A, SLUD::ScalePermstruct_t* scale_perm, 
				SLUD::Pslu_freeable_t* pslu_freeable, type_map::LUstruct_t* lu,
				SLUD::gridinfo_t* grid)
#endif
    {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
      SLUD::Z::zdist_psymbtonum(options, n, A, scale_perm, pslu_freeable, lu, grid);
#else
      SLUD::Z::zdist_psymbtonum(fact, n, A, scale_perm, pslu_freeable, lu, grid);
#endif
    }

    static double plangs(char* norm, SLUD::SuperMatrix* A, SLUD::gridinfo_t* grid)
    {
      return SLUD::Z::pzlangs(norm, A, grid);
    }

    static void SolveInit(SLUD::amesos2_superlu_dist_options_t* options, SLUD::SuperMatrix* A,
			  SLUD::int_t* perm_r, SLUD::int_t* perm_c, SLUD::int_t nrhs,
			  type_map::LUstruct_t* lu, SLUD::gridinfo_t* grid, 
			  type_map::SOLVEstruct_t* solve_struct)
    {
      SLUD::Z::zSolveInit(options, A, perm_r, perm_c, nrhs, lu, grid, solve_struct);
    }

    static void LUstructInit(SLUD::int_t m, SLUD::int_t n, type_map::LUstruct_t* lu)
    {
      /// When we make sure that version 5 and higher is used
      /// we do not perform runtime check of the interface
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::Z::zLUstructInit(n, lu);
#elif defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
      SLUD::Z::LUstructInit(n, lu);
#else
#ifdef HAVE_SUPERLUDIST_LUSTRUCTINIT_2ARG
      SLUD::Z::LUstructInit(n, lu);
#else
      SLUD::Z::LUstructInit(m, n, lu);
#endif
#endif
    }

    static void Destroy_LU(SLUD::int_t m, SLUD::gridinfo_t* grid, type_map::LUstruct_t* lu)
    {
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::Z::zDestroy_LU(m, grid, lu);
#else
      SLUD::Z::Destroy_LU(m, grid, lu);
#endif
    }

    static void LUstructFree(type_map::LUstruct_t* lu)
    {
#if SUPERLU_DIST_MAJOR_VERSION > 6 || (SUPERLU_DIST_MAJOR_VERSION == 6 && SUPERLU_DIST_MINOR_VERSION > 2)
      SLUD::Z::zLUstructFree(lu);
#else
      SLUD::Z::LUstructFree(lu);
#endif
    }

    static void SolveFinalize(SLUD::amesos2_superlu_dist_options_t* options,
			      type_map::SOLVEstruct_t* solve_struct)
    {
      SLUD::Z::zSolveFinalize(options, solve_struct);
    }

    static int ldperm_dist(int job, int n, SLUD::int_t nnz, SLUD::int_t colptr[], SLUD::int_t adjncy[],
            SLUD::Z::doublecomplex nzval[], SLUD::int_t *perm, double u[], double v[])
    {
      return SLUD::Z::zldperm_dist(job, n, nnz, colptr, adjncy, nzval, perm, u, v);
    }
  };
#endif	// HAVE_TEUCHOS_COMPLEX

} // end namespace Amesos2


#endif  // AMESOS2_SUPERLUDIST_FUNCTIONMAP_HPP
