/**
   \file   Amesos2_FunctionMap.hpp
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
} // end namespace SLU


namespace Amesos {

template <typename Scalar>
Teuchos::RCP<SLU::SuperMatrix> getPermMatrix(
  SLU::superlu_options_t* options,
  SLU::SuperMatrix* A,          // A is in NRFormat
  int* perm_c,
  int* etree)
{
  SLU::SuperMatrix AA;      // A in SLU_NC format

  // AA in SLU_NCP format after preordering
  Teuchos::RCP<SLU::SuperMatrix> AC = Teuchos::rcp(new SLU::SuperMatrix());

  if( A->Stype == SLU::SLU_NR ){
    SLU::NRformat* Astore = (SLU::NRformat*)A->Store;
    // SLU::S::sCreate_CompCol_Matrix(&AA, A->ncol, A->nrow, Astore->nnz,
    //   (typename TypeMap<Amesos::Superlu,Scalar>::type*)Astore->nzval,
    //   Astore->colind, Astore->rowptr, SLU::SLU_NC, A->Dtype, A->Mtype);
    FunctionMap<Superlu,Scalar>::create_CompCol_Matrix(&AA, A->ncol, A->nrow,
      Astore->nnz, (typename TypeMap<Amesos::Superlu,Scalar>::type*)Astore->nzval,
      Astore->colind, Astore->rowptr, SLU::SLU_NC, A->Dtype, A->Mtype);
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
  if ( permc_spec != SLU::MY_PERMC ){
    SLU::get_perm_c(permc_spec, &AA, perm_c);
  }
  SLU::sp_preorder(options, &AA, perm_c, etree, AC.getRawPtr());

  // cleanup
  SLU::Destroy_SuperMatrix_Store(&AA);

  return(AC);
}



/** 
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
 */
template <typename Scalar>
struct FunctionMap<Superlu,Scalar>
{
  /** \brief Binds to the appropriate Superlu solver driver based on data type
   * 
   * \throw std::runtime_error If no specialization of this type exists for a
   *        particular scalar type
   */
  static void gssvx(
    SLU::superlu_options_t*,    ///< options
    SLU::SuperMatrix*,          ///< A      
    int*,                       ///< perm_c 
    int*,                       ///< perm_r 
    int*,                       ///< etree  
    char*,                      ///< equed  
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< R      
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< C      
    SLU::SuperMatrix*,          ///< L      
    SLU::SuperMatrix*,          ///< U      
    void*,                      ///< work   
    int,                        ///< lwork
    SLU::SuperMatrix*,          ///< B
    SLU::SuperMatrix*,          ///< X
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< recip_pivot_growth
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< rcond
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< ferr
    typename TypeMap<Superlu,Scalar>::magnitude_type*, ///< berr
    SLU::mem_usage_t*,          ///< mem_usage
    SLU::SuperLUStat_t*,        ///< stat
    int*                        ///< info
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }

  /** \brief Computes an LU factorization of a general m-by-n matrix
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
    SLU::superlu_options_t*,    ///< options
    SLU::SuperMatrix*,          ///< A
    int,                        ///< relax
    int,                        ///< panel_size
    int*,                       ///< elimination tree
    void*,                      ///< work
    int,                        ///< lwork
    int*,                       ///< perm_c
    int*,                       ///< perm_r
    SLU::SuperMatrix*,          ///< L
    SLU::SuperMatrix*,          ///< U
    SLU::SuperLUStat_t*,        ///< stat
    int*                        ///< info
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }
  

  /** \brief Creates a Superlu CCS matrix using the appropriate function
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompCol_Matrix(
    SLU::SuperMatrix*,          ///< A
    int,                        ///< m => number of rows
    int,                        ///< n => number of cols
    int,                        ///< nnz
    typename TypeMap<Superlu,Scalar>::type*, ///< nzval
    int*,                       ///< rowind
    int*,                       ///< colptr
    SLU::Stype_t,               ///< SLU Storage type
    SLU::Dtype_t,               ///< SLU data type
    SLU::Mtype_t                ///< SLU matrix type
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }

  /** \brief Creates a Superlu CRS matrix using the appropriate function
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompRow_Matrix(
    SLU::SuperMatrix*,          ///< A
    int,                        ///< m => number of rows
    int,                        ///< n => number of cols
    int,                        ///< nnz
    typename TypeMap<Superlu,Scalar>::type*, ///< nzval
    int*,                       ///< rowind
    int*,                       ///< colptr
    SLU::Stype_t,               ///< SLU Storage type
    SLU::Dtype_t,               ///< SLU data type
    SLU::Mtype_t                ///< SLU matrix type
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }


  /** \brief Creates a Superlu Dense Matrix using the appropriate Superlu
   *         function.
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_Dense_Matrix(
    SLU::SuperMatrix*,          ///< X
    int,                        ///< m => number of rows
    int,                        ///< n => number of cols
    typename TypeMap<Superlu,Scalar>::type*, ///< x => vals in column major order
    int,                        ///< ldx => leading dimension of x
    SLU::Stype_t,               ///< SLU storage type
    SLU::Dtype_t,               ///< SLU data type
    SLU::Mtype_t                ///< SLU matrix type
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "Superlu does not support the data type");
    }
};


/* ==================== Specializations ==================== */

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
      // TEST_FOR_EXCEPTION( options == NULL,
      //   std::runtime_error,
      //   "options pointer is null!");
      // TEST_FOR_EXCEPTION( A == NULL,
      //   std::runtime_error,
      //   "Supermatrix pointer A is null!");
      // TEST_FOR_EXCEPTION( perm_c == NULL,
      //   std::runtime_error,
      //   "perm_c pointer is null!");
      // TEST_FOR_EXCEPTION( perm_r == NULL,
      //   std::runtime_error,
      //   "perm_r pointer is null!");
      // TEST_FOR_EXCEPTION( etree == NULL,
      //   std::runtime_error,
      //   "etree pointer is null!");
      // TEST_FOR_EXCEPTION( R == NULL,
      //   std::runtime_error,
      //   "R pointer is null!");
      // TEST_FOR_EXCEPTION( C == NULL,
      //   std::runtime_error,
      //   "C pointer is null!");
      // TEST_FOR_EXCEPTION( L == NULL,
      //   std::runtime_error,
      //   "Supermatrix L pointer is null!");
      // TEST_FOR_EXCEPTION( U == NULL,
      //   std::runtime_error,
      //   "Supermatrix U pointer is null!");
      // TEST_FOR_EXCEPTION( X == NULL,
      //   std::runtime_error,
      //   "Supermatrix X pointer is null!");
      // TEST_FOR_EXCEPTION( B == NULL,
      //   std::runtime_error,
      //   "Supermatrix B pointer is null!");
      // TEST_FOR_EXCEPTION( recip_pivot_growth == NULL,
      //   std::runtime_error,
      //   "recip_pivot_growth pointer is null!");
      // TEST_FOR_EXCEPTION( rcond == NULL,
      //   std::runtime_error,
      //   "rcond pointer is null!");
      // TEST_FOR_EXCEPTION( ferr == NULL,
      //   std::runtime_error,
      //   "ferr pointer is null!");
      // TEST_FOR_EXCEPTION( berr == NULL,
      //   std::runtime_error,
      //   "berr pointer is null!");
      
      SLU::S::sgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
        lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {

      Teuchos::RCP<SLU::SuperMatrix> AC = getPermMatrix<float>(options,A,perm_c,etree);

      SLU::S::sgstrf(options, AC.getRawPtr(), relax, panel_size, etree, work, lwork, perm_c,
        perm_r, L, U, stat, info);

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
      TEST_FOR_EXCEPTION( options == NULL,
        std::runtime_error,
        "options pointer is null!");
      TEST_FOR_EXCEPTION( A == NULL,
        std::runtime_error,
        "Supermatrix pointer A is null!");
      TEST_FOR_EXCEPTION( perm_c == NULL,
        std::runtime_error,
        "perm_c pointer is null!");
      TEST_FOR_EXCEPTION( perm_r == NULL,
        std::runtime_error,
        "perm_r pointer is null!");
      TEST_FOR_EXCEPTION( etree == NULL,
        std::runtime_error,
        "etree pointer is null!");
      TEST_FOR_EXCEPTION( R == NULL,
        std::runtime_error,
        "R pointer is null!");
      TEST_FOR_EXCEPTION( C == NULL,
        std::runtime_error,
        "C pointer is null!");
      TEST_FOR_EXCEPTION( L == NULL,
        std::runtime_error,
        "Supermatrix L pointer is null!");
      TEST_FOR_EXCEPTION( U == NULL,
        std::runtime_error,
        "Supermatrix U pointer is null!");
      TEST_FOR_EXCEPTION( X == NULL,
        std::runtime_error,
        "Supermatrix X pointer is null!");
      TEST_FOR_EXCEPTION( B == NULL,
        std::runtime_error,
        "Supermatrix B pointer is null!");
      TEST_FOR_EXCEPTION( recip_pivot_growth == NULL,
        std::runtime_error,
        "recip_pivot_growth pointer is null!");
      TEST_FOR_EXCEPTION( rcond == NULL,
        std::runtime_error,
        "rcond pointer is null!");
      TEST_FOR_EXCEPTION( ferr == NULL,
        std::runtime_error,
        "ferr pointer is null!");
      TEST_FOR_EXCEPTION( berr == NULL,
        std::runtime_error,
        "berr pointer is null!");
      
      SLU::D::dgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
        lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = getPermMatrix<float>(options,A,perm_c,etree);

      SLU::D::dgstrf(options, AC.getRawPtr(), relax, panel_size, etree, work, lwork, perm_c,
        perm_r, L, U, stat, info);

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
      TEST_FOR_EXCEPTION( options == NULL,
        std::runtime_error,
        "options pointer is null!");
      TEST_FOR_EXCEPTION( A == NULL,
        std::runtime_error,
        "Supermatrix pointer A is null!");
      TEST_FOR_EXCEPTION( perm_c == NULL,
        std::runtime_error,
        "perm_c pointer is null!");
      TEST_FOR_EXCEPTION( perm_r == NULL,
        std::runtime_error,
        "perm_r pointer is null!");
      TEST_FOR_EXCEPTION( etree == NULL,
        std::runtime_error,
        "etree pointer is null!");
      TEST_FOR_EXCEPTION( R == NULL,
        std::runtime_error,
        "R pointer is null!");
      TEST_FOR_EXCEPTION( C == NULL,
        std::runtime_error,
        "C pointer is null!");
      TEST_FOR_EXCEPTION( L == NULL,
        std::runtime_error,
        "Supermatrix L pointer is null!");
      TEST_FOR_EXCEPTION( U == NULL,
        std::runtime_error,
        "Supermatrix U pointer is null!");
      TEST_FOR_EXCEPTION( X == NULL,
        std::runtime_error,
        "Supermatrix X pointer is null!");
      TEST_FOR_EXCEPTION( B == NULL,
        std::runtime_error,
        "Supermatrix B pointer is null!");
      TEST_FOR_EXCEPTION( recip_pivot_growth == NULL,
        std::runtime_error,
        "recip_pivot_growth pointer is null!");
      TEST_FOR_EXCEPTION( rcond == NULL,
        std::runtime_error,
        "rcond pointer is null!");
      TEST_FOR_EXCEPTION( ferr == NULL,
        std::runtime_error,
        "ferr pointer is null!");
      TEST_FOR_EXCEPTION( berr == NULL,
        std::runtime_error,
        "berr pointer is null!");
      
      SLU::C::cgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
        lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = getPermMatrix<float>(options,A,perm_c,etree);
      
      SLU::C::cgstrf(options, AC.getRawPtr(), relax, panel_size, etree, work, lwork, perm_c,
        perm_r, L, U, stat, info);

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
      TEST_FOR_EXCEPTION( options == NULL,
        std::runtime_error,
        "options pointer is null!");
      TEST_FOR_EXCEPTION( A == NULL,
        std::runtime_error,
        "Supermatrix pointer A is null!");
      TEST_FOR_EXCEPTION( perm_c == NULL,
        std::runtime_error,
        "perm_c pointer is null!");
      TEST_FOR_EXCEPTION( perm_r == NULL,
        std::runtime_error,
        "perm_r pointer is null!");
      TEST_FOR_EXCEPTION( etree == NULL,
        std::runtime_error,
        "etree pointer is null!");
      TEST_FOR_EXCEPTION( R == NULL,
        std::runtime_error,
        "R pointer is null!");
      TEST_FOR_EXCEPTION( C == NULL,
        std::runtime_error,
        "C pointer is null!");
      TEST_FOR_EXCEPTION( L == NULL,
        std::runtime_error,
        "Supermatrix L pointer is null!");
      TEST_FOR_EXCEPTION( U == NULL,
        std::runtime_error,
        "Supermatrix U pointer is null!");
      TEST_FOR_EXCEPTION( X == NULL,
        std::runtime_error,
        "Supermatrix X pointer is null!");
      TEST_FOR_EXCEPTION( B == NULL,
        std::runtime_error,
        "Supermatrix B pointer is null!");
      TEST_FOR_EXCEPTION( recip_pivot_growth == NULL,
        std::runtime_error,
        "recip_pivot_growth pointer is null!");
      TEST_FOR_EXCEPTION( rcond == NULL,
        std::runtime_error,
        "rcond pointer is null!");
      TEST_FOR_EXCEPTION( ferr == NULL,
        std::runtime_error,
        "ferr pointer is null!");
      TEST_FOR_EXCEPTION( berr == NULL,
        std::runtime_error,
        "berr pointer is null!");
      
      SLU::Z::zgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work,
        lwork, B, X, recip_pivot_growth, rcond, ferr, berr, mem_usage, stat, info);
    }

  static void gstrf(SLU::superlu_options_t* options, SLU::SuperMatrix* A,
    int relax, int panel_size, int* etree, void* work, int lwork, int* perm_c,
    int* perm_r, SLU::SuperMatrix* L, SLU::SuperMatrix* U,
    SLU::SuperLUStat_t* stat, int* info)
    {
      Teuchos::RCP<SLU::SuperMatrix> AC = getPermMatrix<float>(options,A,perm_c,etree);
      
      SLU::Z::zgstrf(options, AC.getRawPtr(), relax, panel_size, etree, work, lwork, perm_c,
        perm_r, L, U, stat, info);

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




} // end namespace Amesos

#endif  // AMESOS2_SUPERLU_FUNCTIONMAP_HPP
