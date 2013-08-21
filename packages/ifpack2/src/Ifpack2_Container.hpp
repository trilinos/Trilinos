/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CONTAINER_HPP
#define IFPACK2_CONTAINER_HPP

/// \file Ifpack2_Container.hpp
/// \brief Ifpack2::Container class declaration

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Describable.hpp>
#include <iostream>

namespace Ifpack2 {

/// \class Container
/// \brief Interface for creating and solving a local linear problem.
///
/// This class is mainly useful for BlockRelaxation and other
/// operators that need to solve linear systems with diagonal blocks
/// of a sparse matrix.  Users of BlockRelaxation do not normally need
/// to interact with the Container interface.  However, they do need
/// to specify a specific Container subclass to use.  Implementations
/// of Container specify
/// - the kind of data structure used to store diagonal blocks, and
/// - how linear systems with those diagonal blocks are solved.
///
/// For example, the SparseContainer subclass uses a Tpetra sparse
/// matrix to store each diagonal block, and can use any given Ifpack2
/// Preconditioner subclass to solve linear systems.
///
/// A Container can create, populate, and solve a local linear
/// system. The local linear system matrix, B, is a submatrix of the
/// local components of a distributed matrix, A. The idea of Container
/// is to specify the rows of A that are contained in B, then extract
/// B from A, and compute all it is necessary to solve a linear system
/// in B.  Then, set the initial guess (if necessary) and right-hand
/// side for B, and solve the linear system in B.
///
/// If you are writing a class (comparable to BlockRelaxation) that
/// uses Container, you should use it in the following way:
/// <ol>
/// <li> Create an container object, specifying the number of rows of B.</li>
/// <li> Optionally, set linear solve parameters using setParameters().</li>
/// <li> Initialize the container by calling initialize().</li>
/// <li> Specify the ID of the local rows of A that are contained in B,
///      using ID().</li>
/// <li> Prepare the linear system solver using compute().</li>
/// <li> Solve the linear system using apply().</li>
/// </ol>
///
/// \c MatrixType and \c InverseType may store values of different
/// types, and may have different template parameters (e.g., local or
/// global ordinal types).  You may mix and match so long as implicit
/// conversions are available.  The most obvious use case for this
/// are:
/// - <tt>MatrixType::global_ordinal_type=long long</tt> and
///   <tt>InverseType::global_ordinal_type=short</tt>
/// - <tt>MatrixType::scalar_type=float</tt> and
///   <tt>InverseType::scalar_type=double</tt>
template<class MatrixType, class InverseType>
class Container : public Teuchos::Describable {
public:
  typedef typename MatrixType::scalar_type          MatrixScalar;
  typedef typename MatrixType::local_ordinal_type   MatrixLocalOrdinal;
  typedef typename MatrixType::global_ordinal_type  MatrixGlobalOrdinal;
  typedef typename MatrixType::node_type            MatrixNode;

  typedef typename InverseType::scalar_type         InverseScalar;
  typedef typename InverseType::local_ordinal_type  InverseLocalOrdinal;
  typedef typename InverseType::global_ordinal_type InverseGlobalOrdinal;
  typedef typename InverseType::node_type           InverseNode;

  typedef typename MatrixType::scalar_type          scalar_type;
  typedef typename MatrixType::local_ordinal_type   local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type  global_ordinal_type;
  typedef typename MatrixType::node_type            node_type;

  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  //! Destructor.
  virtual ~Container() {};

  //! The number of rows in the diagonal block.
  virtual size_t getNumRows() const = 0;

  /// \brief A nonconst reference to the local index associated with local row i.
  ///
  /// The set of (local) rows assigned to this container is defined
  /// by calling ID(i) = j, where i (from 0 to NumRows()) indicates
  /// the container-row, and j indicates the local row in the calling
  /// process.
  ///
  /// This method is usually used to recorder the local row index (on
  /// the calling process) of the i-th row in the container.
  virtual MatrixLocalOrdinal & ID (const size_t i) = 0;

  //! Initialize the container, by performing all operations that only require matrix structure.
  virtual void initialize () = 0;

  //! Finalize the linear system matrix and prepare for the application of the inverse.
  virtual void compute (const Teuchos::RCP<const row_matrix_type>& Matrix) = 0;

  //! Set parameters.
  virtual void setParameters (const Teuchos::ParameterList& List) = 0;

  //! Return \c true if the container has been successfully initialized.
  virtual bool isInitialized () const = 0;

  //! Return \c true if the container has been successfully computed.
  virtual bool isComputed () const = 0;

  /// \brief Compute Y = alpha * M^{-1} X + beta*Y
  ///
  /// X and Y have as many global rows as the dimension of the global
  /// problem from which the container was originally extracted.
  virtual void
  apply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
         Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         MatrixScalar alpha = Teuchos::ScalarTraits<MatrixScalar>::one (),
         MatrixScalar beta = Teuchos::ScalarTraits<MatrixScalar>::zero ()) = 0;

  /// \brief Compute Y = alpha * diag(D) * M^{-1} (diag(D) X) + beta*Y
  ///
  /// D, X, and Y have as many global rows as the dimension of the
  /// global problem from which the container was originally
  /// extracted.  The method is designed to support techniques with
  /// overlap, such as Schwarz methods.
  virtual void
  weightedApply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
                 Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
                 const Tpetra::Vector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& D,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 MatrixScalar alpha = Teuchos::ScalarTraits<MatrixScalar>::one (),
                 MatrixScalar beta = Teuchos::ScalarTraits<MatrixScalar>::zero ()) = 0;

  //! Print basic information about the container to \c os.
  virtual std::ostream& print (std::ostream& os) const = 0;
};

template <class MatrixType, class InverseType>
inline std::ostream&
operator<< (std::ostream& os, const Ifpack2::Container<MatrixType,InverseType>& obj)
{
  return obj.print (os);
}

}

#endif // IFPACK2_CONTAINER_HPP
