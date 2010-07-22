/**
   \file   NewSolver_MatrixHelper.hpp
   \author John Doe <jd@sandia.gov>
   \date   Fri Jul  9 13:44:44 2010
   
   \brief  Template file for defining helper methods for matrix creation.  This
           is just a skeleton.  There is not guarantee that the setup will be
           useful for a new solver package.
*/


#ifndef AMESOS2_NEWSOLVER_MATRIXHELPER_HPP
#define AMESOS2_NEWSOLVER_MATRIXHELPER_HPP

#include <Teuchos_ArrayRCP.hpp>

#include "Amesos2_MatrixHelper.hpp"
#include "Amesos2_NewSolver_FunctionMap.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace New_Solver {
extern "C" {
// #include files that define matrix types for this solver

// Forward declaration just for example.
class Matrix;

}
}

// Solver class must be forward-declared due to circular reference
template <class Matrix, class Vector> class NewSolver;

namespace Amesos {


template <>
struct MatrixHelper<NewSolver>
{

  /**
   * \brief Creates a NewSolver compressed-row Matrix from the given Matrix
   *
   * \tparam Matrix A matrix type conforming to the interface, in particular
   *         an Amesos::MatrixAdapter<>.
   *
   * \param [in]     mat    The matrix which will be converted to NewSolver format
   * \param [in,out] nzval  A user-provided persisting store for the nonzero
   *                        values of the matrix.  The NewSolver matrix will expect
   *                        the array to persist throughout its existence.
   * \param [in,out] colind A user-provided persisting store for the column
   *                        indices of the matrix.
   * \param [in,out] rowptr User-privded persisting store for row pointers
   * \param [out]    A      Pointer to the NewSolver SuperMatrix which is to be constructed
   */
  template <class Matrix>
  static void createCRSMatrix(
    const Teuchos::Ptr<Matrix>& mat,  
    const Teuchos::ArrayView<typename TypeMap<NewSolver,typename Matrix::scalar_type>::type>& nzval,
    const Teuchos::ArrayView<int>& colind, 
    const Teuchos::ArrayView<int>& rowptr, 
    const Teuchos::Ptr<New_Solver::Matrix>& A 
    )                         
    {
      typedef typename Matrix::scalar_type                          scalar_type;
      typedef typename Matrix::global_ordinal_type                      go_type;
      typedef typename Matrix::global_size_type                         gs_type;
      typedef typename TypeMap<Amesos::NewSolver,scalar_type>::type solver_type;

      // Extract the necessary information from mat and call NewSolver function
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
      
      Array<scalar_type> nzval_tmp(nzval.size());
      Array<go_type> colind_tmp(colind.size());
      Array<gs_type> rowptr_tmp(rowptr.size());
      size_t nnz_ret = 0;

      // Actually uses the compressed-row-store format, and tells NewSolver this
      // while creating a compressed-column store.
      mat->getCrs((ArrayView<scalar_type>)nzval_tmp,
        (ArrayView<go_type>)colind_tmp, (ArrayView<gs_type>)rowptr_tmp, nnz_ret);

      TEST_FOR_EXCEPTION( nnz_ret != Teuchos::as<size_t>(nnz),
        std::runtime_error,
        "Number of nonzeros returned by getCrs() different from getGlobalNNZ()");
      
      /* Convert types
       *
       * Note: We cannot simply convert when necessary.  That is, we cannot
       * check whether the matrix types are double, int, and int, and then do
       * a simple ArrayView::assign if they are the same, since the compiler
       * still has to do the static linking and cannot find an appropriate
       * function call in the case that the types are not double, int, and int
       */
      {
        typename Array<scalar_type>::size_type i, size;
        size = nzval_tmp.size();
        for ( i = 0; i < size; ++i ){
          nzval[i] = Teuchos::as<solver_type>(nzval_tmp[i]);
        }
      }
      {
        typename Array<go_type>::size_type i, size;
        size = colind_tmp.size();
        for ( i = 0; i < size; ++i ){
          colind[i] = Teuchos::as<int>(colind_tmp[i]);
        }
      }
      {
        typename Array<gs_type>::size_type i, size;
        size = rowptr_tmp.size();
        for ( i = 0; i < size; ++i ){
          rowptr[i] = Teuchos::as<int>(rowptr_tmp[i]);
        }
      }
      // end conversions

      // Replace with actual function call as defined in
      // NewSolver_FunctionMap.hpp
      FunctionMap<NewSolver,scalar_type>::create_CCS_Matrix( /* args */ );
    }


  /**
   * \brief Creates a NewSolver Dense Matrix from the given MultiVector
   *
   * \tparam MV A multivector type conforming to the interface, in particular
   *         an Amesos::MultiVecAdapter<>.
   *
   * \param [in]     mv   The MultiVector which will be converted to NewSolver format
   * \param [in,out] vals A user-provided persisting store for the values of
   *                      the multivector.  The NewSolver matrix will expect the array
   *                      to persist throughout its existence.
   * \param [out]    ldx  Leading dimension of \c vals
   * \param [out]    X    Pointer to the NewSolver Dense Matrix which is to be
   *                      constructed
   */
  template <class MV>
  static void createMVDenseMatrix(
    const Teuchos::Ptr<MV>& mv,
    const Teuchos::ArrayView<typename TypeMap<NewSolver,typename MV::scalar_type>::type>& vals,
    int& ldx,
    const Teuchos::Ptr<New_Solver::Matrix>& X
    )
    {
      typedef typename MV::scalar_type scalar_type;
      typedef typename TypeMap<NewSolver,scalar_type>::type solver_type;

      int rows, cols;
      rows = Teuchos::as<int>(mv->getGlobalLength());
      cols = Teuchos::as<int>(mv->getGlobalNumVectors());
      ldx  = Teuchos::as<int>(mv->getStride());

      if ( Util::is_same<scalar_type,solver_type>::value ){
        mv->get1dCopy(vals,ldx);
      } else {
        int vals_length = rows * cols;
        const Teuchos::Array<scalar_type> vals_tmp(vals_length);
        mv->get1dCopy(vals_tmp.view(0, vals_length));
        for ( int i = 0; i < vals_length; ++i ){
          vals[i] = Teuchos::as<solver_type>(vals_tmp[i]);
        }
      }

      // Again, replace with the actual function call as it is defined in
      // NewSolver_FunctionMap.hpp
      FunctionMap<NewSolver,scalar_type>::create_Dense_Matrix( /* args */ );
    }


  /**
   * \brief Creates a NewSolver Dense Matrix from the given MultiVector
   *
   * \tparam MV A multivector type conforming to the interface, in particular
   *         an Amesos::MultiVecAdapter<>.
   *
   * \param [in]     mv   The MultiVector which will be converted to NewSolver format
   * \param [out]    X    Pointer to the NewSolver Dense SuperMatrix which is to be
   *                      constructed
   *
   * \return A Teuchos::ArrayRCP pointing to the beginning of a contiguous
   * store of the values in \c X , which is <b>not</b> necessarily the beginning of
   * the contiguous store of values in \c mv .
   */
  template <class MV>
  static
  Teuchos::ArrayRCP<typename TypeMap<NewSolver,typename MV::scalar_type>::type>
  createMVDenseMatrix(
    const Teuchos::Ptr<MV>& mv,
    const Teuchos::Ptr<New_Solver::Matrix>& X
    )
    {
      typedef typename MV::scalar_type scalar_type;
      typedef typename TypeMap<NewSolver,scalar_type>::type solver_type;

      int rows, cols, ldx;
      rows = Teuchos::as<int>(mv->getGlobalLength());
      cols = Teuchos::as<int>(mv->getGlobalNumVectors());
      ldx  = Teuchos::as<int>(mv->getStride());

      Teuchos::ArrayRCP<scalar_type> vals_ptr;

      vals_ptr = mv->get1dViewNonConst();
      typedef typename Teuchos::ArrayRCP<scalar_type>::size_type size_type;
      size_type vals_length = vals_ptr.size();

      typedef typename Teuchos::ArrayRCP<solver_type>::size_type solver_size_type;
      Teuchos::ArrayRCP<solver_type> solver_vals(Teuchos::as<solver_size_type>(vals_length));

      // Convert value types
      for ( size_type i = 0; i < vals_length; ++i ){
        solver_vals[i] = Teuchos::as<solver_type>(vals_ptr[i]);
      }

      // Replace with actual function call, as defined
      FunctionMap<NewSolver,scalar_type>::create_Dense_Matrix( /* args */ );

      return solver_vals;
    }
};                              // end struct MatrixHelper


} // end namespace Amesos

#endif  // end AMESOS2_NEWSOLVER_MATRIXHELPER_HPP
