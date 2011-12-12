/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP
#define KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP

#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_BLAS_types.hpp>

/** \file Kokkos_DummySparseKernelClass.hpp
    \brief A file containing a stub for a new sparse kernel provider, as outlined in the \ref kokkos_crs_ops "Kokkos CRS API".
 */

namespace KokkosExamples {

  /** \class DummySparseKernel 
      \ingroup kokkos_crs_ops
      \brief A dummy-class illustrating the components necessary for a %Kokkos sparse operations provider.
   */
  template <class Node> 
  class DummySparseKernel {
  public:
    //@{ 
    //! @name Typedefs and structs

    //!
    typedef void  ScalarType;
    //!
    typedef void OrdinalType;
    //!
    typedef Node    NodeType;
    //! 
    typedef DummySparseKernel<Node> ThisType;

    /** \brief Rebind struct, for specifying type information for a different scalar.
        
        For a scalar type \c T, specify the class type that should be used for
        sparse kernels over sparse matrices of type \c T.
      */
    template <class T>
    struct rebind {
      typedef DummySparseKernel<Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! \brief Constructor accepting and retaining a node object.
    DummySparseKernel(const Teuchos::RCP<Node> & node) : node_(node) {}

    //! \brief Destructor.
    ~DummySparseKernel() {}

    //@} 
    //! @name Accessor routines.
    //@{ 

    /** \brief Node accessor.
        
        This returns the node type of type this::NodeType that was passed at object construction.
      */
    Teuchos::RCP<Node> getNode() const {return node_;}

    //@}
    //! @name Initialization of structure
    //@{

    /** \brief Initialize the structure of the sparse matrix.
        
        This is the mechanism by which the user specifies the structure for the sparse matrix. It always comes via 
        a Kokkos::CrsGraph<O,N,SO>, where
        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        The graph structure must be provided via initializeStructure() before the matrix values are specified to 
        initializeValues(). After calling initializeStructure(), the clear() method must be called before calling
        initializeStructure() again.

        In general, both initializeStructure() and initializeValues() must be called before calling multiply() or solve().
      */
    template <class Ordinal>
    void initializeStructure(const Kokkos::CrsGraph<Ordinal,Node,DummySparseKernel<Node> > & /* graph */) {}

    /** \brief Initialize the values of the sparse matrix.

        This is the mechanism by which the user specifies the values for the sparse matrix. It always comes via 
        a Kokkos::CrsMatrix<S,O,N,SO>, where
        - \c S is the scalar type this::ScalarType
        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        The graph structure must be provided via initializeStructure() before the matrix values are specified to 
        initializeValues(). The matrix values can be provided repeatadly via initializeValues() without the need
        to call clear().

        In general, both initializeStructure() and initializeValues() must be called before calling multiply() or solve().
      */
    template <class Scalar, class Ordinal>
    void initializeValues(const Kokkos::CrsMatrix<Scalar,Ordinal,Node,DummySparseKernel<Node> > & /* matrix */) {}

    /** \brief Clear the graph and matrix data.
        
        clear() must be called between successive calls to initializeStructure(). 

        In general, after calling clear(), no significant data (including persisting references) will be preserved, save for the pointer to the node. 

        In general, after calling clear(), multiply() and solve() cannot be called.
      */
    void clear() {}

    //@}
    //! @name Computational methods
    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const Kokkos::MultiVector<DomainScalar,Node> &X, Kokkos::MultiVector<RangeScalar,Node> &Y) const {}

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const Kokkos::MultiVector<DomainScalar,Node> &X, RangeScalar beta, Kokkos::MultiVector<RangeScalar,Node> &Y) const {}

    //! Solves the matrix for a given set of right-hand-sides.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const Kokkos::MultiVector<DomainScalar,Node> &Y, Kokkos::MultiVector<RangeScalar,Node> &X) const {}

    //@}
  protected:
    //! Copy constructor (protected and unimplemented)
    DummySparseKernel(const DummySparseKernel& source);

    Teuchos::RCP<Node> node_;
  };

  /** \example DummySparseKernelDriver.cpp 
    * This is an example that demonstrates the basic use case for a sparse kernel provider. It also verifies that the stub builds correctly.
    */

} // end of namespace KokkosExamples

#endif
