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
