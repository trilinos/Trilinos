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

/**
   \file   Amesos2_Superlu_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Mon May 31 23:38:46 2010

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_SUPERLU_FUNCTIONMAP_HPP
#define AMESOS2_SUPERLU_FUNCTIONMAP_HPP

#include <complex>

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_Superlu_TypeMap.hpp"


// External definitions of the Superlu functions
namespace SLU {
extern "C" {
typedef int int_t;
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"

namespace S {
#include "slu_sdefs.h"          // single-precision real definitions
}

namespace D {
#include "slu_ddefs.h"          // double-precision real definitions
}

namespace C {
#include "slu_cdefs.h"          // single-precision complex definitions
}

namespace Z {
#include "slu_zdefs.h"          // double-precision complex definitions
}

} // end extern "C"


/**
 * \brief Gets the permutation matrix \c perm_c for \c A based on \c options .
 *
 * \param options A Superlu options structure.  Should specify what type of
 *                permutation to use.
 * \param A       A Superlu \c SuperMatrix , with SLU::NRformat Store, that is
 *                to be permuted.
 * \param perm_c  A permutation vector
 * \param etree   An elimination tree for \c A
 *
 * \tparam Scalar the scalar data type of entries in \c A and the returned
 * matrix.
 *
 * \return A Teuchos::RCP pointing to a SLU::SuperMatrix which is the permuted
 * form of \c A .  This new permuted matrix will have a SLU::SLU_NCP
 * (non-supernodal, column, permuted) format, which is appropriate for use in
 * SLU::gstrf().
 */
template <typename Scalar>
Teuchos::RCP<SLU::SuperMatrix> getPermMatrix(
  SLU::superlu_options_t* options,
  SLU::SuperMatrix* A,          // A is in NRFormat
  int* perm_c,
  int* etree);


} // end namespace SLU


namespace Amesos {



/**
 * \brief Pass function calls to Superlu based on data type.
 *
 * Helper class which passes on function calls to the appropriate Superlu
 * function based on the type of its scalar template argument.
 *
 * Superlu has solver and matrix builder functions defined based on data type.
 * One function for complex, one for double precision complex, another for \c
 * float , and yet another for \c double.  To work elegantly with the
 * Amesos::Superlu interface we want to be able to perform a single function
 * call which is appropriate for the scalar type of the Matrix and
 * MultiVectors that we are working with.  The \c FunctionMap class
 * provides that capability.
 *
 * The class template is specialized for each data type that Superlu supports,
 * and errors are thrown for other data types.
 *
 * Please see the <a
 * href="http://crd.lbl.gov/~xiaoye/SuperLU/superlu_ug.pdf">Superlu Users'
 * Guide</a> for more information on the TPL functions.
 */
template <typename Scalar>
struct FunctionMap<Superlu,Scalar>
{
  /**
   * \brief Binds to the appropriate Superlu solver driver based on data type
   *
   * \throw std::runtime_error If no specialization of this type exists for a
   *        particular scalar type
   */
  static void gssvx(
    SLU::superlu_options_t*                  options,
    SLU::SuperMatrix*                        A,
    int*                                     perm_c,
    int*                                     perm_r,
    int*                                     etree,
    char*                                    equed,
    typename TypeMap<Superlu,Scalar>::magnitude_type* R,
    typename TypeMap<Superlu,Scalar>::magnitude_type* C,
    SLU::SuperMatrix*                        L,
    SLU::SuperMatrix*                        U,
    void*                                    work,
    int                                      lwork,
    SLU::SuperMatrix*                        B,
    SLU::SuperMatrix*                        X,
    typename TypeMap<Superlu,Scalar>::magnitude_type* recip_pivot_growth,
    typename TypeMap<Superlu,Scalar>::magnitude_type* rcond,
    typename TypeMap<Superlu,Scalar>::magnitude_type* ferr,
    typename TypeMap<Superlu,Scalar>::magnitude_type* berr,
    SLU::mem_usage_t*                        mem_usage,
    SLU::SuperLUStat_t*                      stat,
    int*                                     info
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
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
   * in the SLU_NC format, so conversion must be done when necessary
   */
  static void gstrf(
    SLU::superlu_options_t*   options,
    SLU::SuperMatrix*         A,
    int                       relax,
    int                       panel_size,
    int*                      etree,
    void*                     work,
    int                       lwork,
    int*                      perm_c,
    int*                      perm_r,
    SLU::SuperMatrix*         L,
    SLU::SuperMatrix*         U,
    SLU::SuperLUStat_t*       stat,
    int*                      info
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }


  /**
   * \brief Creates a Superlu CCS matrix using the appropriate function
   *
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompCol_Matrix(
    SLU::SuperMatrix*                        A,
    int                                      numrows,
    int                                      numcols,
    int                                      nnz,
    typename TypeMap<Superlu,Scalar>::type*  nzval,
    int*                                     rowind,
    int*                                     colptr,
    SLU::Stype_t                             storage_t,
    SLU::Dtype_t                             data_t,
    SLU::Mtype_t                             mat_t
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }

  /**
   * \brief Creates a Superlu CRS matrix using the appropriate function
   *
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompRow_Matrix(
    SLU::SuperMatrix*                        A,
    int                                      numrows,
    int                                      numcols,
    int                                      nnz,
    typename TypeMap<Superlu,Scalar>::type*  nzval,
    int*                                     rowind,
    int*                                     colptr,
    SLU::Stype_t                             storage_t,
    SLU::Dtype_t                             data_t,
    SLU::Mtype_t                             mat_t
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }


  /**
   * \brief Creates a Superlu Dense Matrix using the appropriate Superlu
   *         function.
   *
   * \param X Superlu SuperMatrix that is to be created
   * \param x vals in column major order
   * \param ldx leading dimension of x
   *
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_Dense_Matrix(
    SLU::SuperMatrix*                        X,
    int                                      numrows,
    int                                      numcols,
    typename TypeMap<Superlu,Scalar>::type*  x,
    int                                      ldx,
    SLU::Stype_t                             storage_t,
    SLU::Dtype_t                             data_t,
    SLU::Mtype_t                             mat_t
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }
};


/* ==================== Specializations ====================
 *
 * \cond Superlu_function_specializations 
 */

template <>
struct FunctionMap<Superlu,float>
{
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

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {

      Teuchos::RCP<SLU::SuperMatrix> AC = SLU::getPermMatrix<float>(options,A,perm_c,etree);

      SLU::S::sgstrf(options, AC.getRawPtr(), relax, panel_size, etree,
        work, lwork, perm_c, perm_r, L, U, stat, info);

      SLU::Destroy_CompCol_Permuted(AC.getRawPtr());
    }

  static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,float>::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }

  static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,float>::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }


  static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
    TypeMap<Superlu,float>::type* x, int ldx, SLU::Stype_t stype,
    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::S::sCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }
};


template <>
struct FunctionMap<Superlu,double>
{
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

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = SLU::getPermMatrix<float>(options,A,perm_c,etree);

      SLU::D::dgstrf(options, AC.getRawPtr(), relax, panel_size, etree,
        work, lwork, perm_c, perm_r, L, U, stat, info);

      SLU::Destroy_CompCol_Permuted(AC.getRawPtr());
    }

  static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,double>::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }

  static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,double>::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }

  static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
    TypeMap<Superlu,double>::type* x, int ldx, SLU::Stype_t stype,
    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::D::dCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }
};


/* The specializations for Teuchos::as<> for SLU::complex and
 * SLU::doublecomplex are provided in Amesos2_Superlu_Type.hpp
 */
template <>
struct FunctionMap<Superlu,std::complex<float> >
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

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = SLU::getPermMatrix<float>(options,A,perm_c,etree);

      SLU::C::cgstrf(options, AC.getRawPtr(), relax, panel_size, etree,
        work, lwork, perm_c, perm_r, L, U, stat, info);

      SLU::Destroy_CompCol_Permuted(AC.getRawPtr());
    }

  static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,std::complex<float> >::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }

  static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,std::complex<float> >::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);
    }

  static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
    TypeMap<Superlu,std::complex<float> >::type* x, int ldx, SLU::Stype_t stype,
    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::C::cCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }
};


template <>
struct FunctionMap<Superlu,std::complex<double> >
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

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = SLU::getPermMatrix<float>(options,A,perm_c,etree);

      SLU::Z::zgstrf(options, AC.getRawPtr(), relax, panel_size, etree,
        work, lwork, perm_c, perm_r, L, U, stat, info);

      SLU::Destroy_CompCol_Permuted(AC.getRawPtr());
    }

  static void create_CompCol_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,std::complex<double> >::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);

      TEST_FOR_EXCEPTION( A == NULL,
        std::runtime_error,
        "Supermatrix A not initialized properly!");
    }


  static void create_CompRow_Matrix(SLU::SuperMatrix* A, int m, int n, int nnz,
    TypeMap<Superlu,std::complex<double> >::type* nzval, int* rowind, int* colptr,
    SLU::Stype_t stype, SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_CompRow_Matrix(A, m, n, nnz, nzval, rowind, colptr,
        stype, dtype, mtype);

      TEST_FOR_EXCEPTION( A == NULL,
        std::runtime_error,
        "Supermatrix A not initialized properly!");
    }

  static void create_Dense_Matrix(SLU::SuperMatrix* X, int m, int n,
    TypeMap<Superlu,std::complex<double> >::type* x, int ldx, SLU::Stype_t stype,
    SLU::Dtype_t dtype, SLU::Mtype_t mtype)
    {
      SLU::Z::zCreate_Dense_Matrix(X, m, n, x, ldx, stype, dtype, mtype);
    }
};

/* \endcond Superlu_function_specializations */


} // end namespace Amesos


namespace SLU {

template <typename Scalar>
Teuchos::RCP<SLU::SuperMatrix> getPermMatrix(
  superlu_options_t* options,
  SuperMatrix* A,          // A is in NRFormat
  int* perm_c,
  int* etree)
{
  SuperMatrix AA;      // A in SLU_NC format

  // AA in SLU_NCP format after preordering
  Teuchos::RCP<SuperMatrix> AC = Teuchos::rcp(new SuperMatrix());

  if( A->Stype == SLU_NR ){
    NRformat* Astore = (NRformat*)A->Store;
    // SLU::S::sCreate_CompCol_Matrix(&AA, A->ncol, A->nrow, Astore->nnz,
    //   (typename TypeMap<Amesos::Superlu,Scalar>::type*)Astore->nzval,
    //   Astore->colind, Astore->rowptr, SLU::SLU_NC, A->Dtype, A->Mtype);
    Amesos::FunctionMap<Amesos::Superlu,Scalar>::create_CompCol_Matrix(
      &AA, A->ncol, A->nrow, Astore->nnz,
      (typename Amesos::TypeMap<Amesos::Superlu,Scalar>::type*)Astore->nzval,
      Astore->colind, Astore->rowptr, SLU_NC, A->Dtype, A->Mtype);
  } else {
    AA = *A;
  }

  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = NATURAL:  natural ordering
   *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
   *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
   *   permc_spec = COLAMD:   approximate minimum degree column ordering
   *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
   */
  int permc_spec = options->ColPerm;
  if ( permc_spec != MY_PERMC ){
    get_perm_c(permc_spec, &AA, perm_c);
  }
  sp_preorder(options, &AA, perm_c, etree, AC.getRawPtr());

  // cleanup
  Destroy_SuperMatrix_Store(&AA);

  return(AC);
}

} // end namespace SLU


#endif  // AMESOS2_SUPERLU_FUNCTIONMAP_HPP
