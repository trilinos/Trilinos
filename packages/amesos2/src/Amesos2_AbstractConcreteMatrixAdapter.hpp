// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_ABSTRACTCONCRETEMATRIXADAPTER_HPP
#define AMESOS2_ABSTRACTCONCRETEMATRIXADAPTER_HPP

namespace Amesos2 {

  /**
   * This class provides a statically polymorphic abstract base class
   * for concrete matrix adapters.
   *
   * If a matrix M inherits from another matrix class A, it's adapter
   * should inherit from::
   *
   * \code 
   *   AbstractConcreteMatrixAdapter<A,M>
   * \endcode
   *
   * and the common functionality for matrices inheritting from A
   * should reside in a specialization of
   * AbstractConcreteMatrixAdapter
   *
   * For example, Epetra_CrsMatrix and Epetra_VbrMatrix both inherit
   * from Epetra_RowMatrix.  There is much functionality which is
   * common, but some details, such as object construction is specific
   * to each, so the \c get_impl function, which must create new
   * instances should be adapted differently for each.
   *
   * \code 
   * template < typename M > 
   * AbstractConcreteMatrixAdapter<Epetra_RowMatrix, M>
   *   : MatrixAdapter<M> {
   *
   *   // < ...common functions... >
   *
   *   RCP<MatrixAdapter<M> > get_impl(){
   *     return static_cast<ConcreteMatrixAdapter<M>*>(this)->get_impl();
   *   }
   * }
   * \endcode
   *
   * So now the ConcreteMatrixAdapter specializations for
   * Epetra_CrsMatrix and Epetra_VbrMatrix must only provide the
   * \c get_impl function.
   */
  template <class Abstract, class Matrix>
  class AbstractConcreteMatrixAdapter {};

}

#endif	// AMESOS2_ABSTRACTCONCRETEMATRIXADAPTER_HPP
