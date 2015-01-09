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
#ifndef IFPACK2_TRIDICONTAINER_DECL_HPP
#define IFPACK2_TRIDICONTAINER_DECL_HPP

/// \file Ifpack2_TriDiContainer_decl.hpp
/// \brief Ifpack2::TriDiContainer class declaration

#include "Ifpack2_Container.hpp"
#include "Ifpack2_Details_MultiVectorLocalGatherScatter.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialTriDiMatrix.hpp"

namespace Ifpack2 {

/// \class TriDiContainer
/// \brief Store and solve a local TriDi linear problem.
///
/// Please refer to the documentation of the Container
/// interface. Currently, Containers are used by BlockRelaxation.
/// Block relaxations need to be able to do two things:
/// <ol>
/// <li> Store the diagonal blocks </li>
/// <li> Solve linear systems with each diagonal block </li>
/// </ol>
/// TriDiContainer stores the diagonal blocks as TriDi matrices, and
/// solves them using either LAPACK (for the four Scalar types that it
/// supports) or a custom LU factorization (for Scalar types not
/// supported by LAPACK).
///
/// As with Ifpack2::Container, <tt>MatrixType</tt> must be a
/// specialization of Tpetra::RowMatrix or of its subclass
/// Tpetra::CrsMatrix.  Using a TriDi matrix for each block is a good
/// idea when the blocks are small.  For large and / or sparse blocks,
/// it would probably be better to use an implementation of Container
/// that stores the blocks sparsely, in particular SparseContainer.
///
/// This class may store the TriDi local matrix using values of a
/// different type (\c LocalScalarType) than those in \c MatrixType.
/// You may mix and match so long as implicit conversions are
/// available between \c LocalScalarType and
/// <tt>MatrixType::scalar_type</tt>.
///
/// This class currently assumes the following about the column and
/// row Maps of the input matrix:
/// <ol>
/// <li> On all processes, the column and row Maps begin with the same
///      set of on-process entries, in the same order.  That is,
///      on-process row and column indices are the same.</li>
/// <li> On all processes, all off-process indices in the column Map
///      of the input matrix occur after that initial set.</li>
/// </ol>
/// These assumptions may be violated if \c MatrixType is a
/// Tpetra::CrsMatrix specialization and was constructed with a
/// user-provided column Map.  The assumptions are not mathematically
/// necessary and could be relaxed at any time.  Implementers who wish
/// to do so will need to modify the extract() method, so that it
/// translates explicitly between local row and column indices,
/// instead of just assuming that they are the same.
template<typename MatrixType, typename LocalScalarType>
class TriDiContainer : public Container<MatrixType> {
public:
  //! \name Public typedefs
  //@{

  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization or a
  /// Tpetra::CrsMatrix specialization.  It may have entirely
  /// different template parameters (e.g., \c scalar_type) than
  /// <tt>InverseType</tt>.
  typedef MatrixType matrix_type;
  //! The second template parameter of this class.
  typedef LocalScalarType local_scalar_type;

  //! The type of entries in the input (global) matrix.
  typedef typename MatrixType::scalar_type scalar_type;
  //! The type of local indices in the input (global) matrix.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  //! The type of global indices in the input (global) matrix.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  //! The Node type of the input (global) matrix.
  typedef typename MatrixType::node_type node_type;

  /// \brief The (base class) type of the input matrix.
  ///
  /// The input matrix to the constructor may be either a
  /// Tpetra::RowMatrix specialization or a Tpetra::CrsMatrix
  /// specialization.  However, we want to make the constructor as
  /// general as possible, so we always accept the matrix as a
  /// Tpetra::RowMatrix.  This typedef is the appropriate
  /// specialization of Tpetra::RowMatrix.
  typedef typename Container<MatrixType>::row_matrix_type row_matrix_type;

  //@}
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor.
  ///
  /// \brief matrix [in] The original input matrix.  This Container
  ///   will construct a local diagonal block from the rows given by
  ///   <tt>localRows</tt>.
  ///
  /// \param localRows [in] The set of (local) rows assigned to this
  ///   container.  <tt>localRows[i] == j</tt>, where i (from 0 to
  ///   <tt>getNumRows() - 1</tt>) indicates the SparseContainer's
  ///   row, and j indicates the local row in the calling process.
  ///   <tt>localRows.size()</tt> gives the number of rows in the
  ///   local matrix on each process.  This may be different on
  ///   different processes.
  TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                  const Teuchos::ArrayView<const local_ordinal_type>& localRows);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~TriDiContainer ();

  //@}
  //! \name Get and set methods
  //@{

  /// \brief The number of rows in the local matrix on the calling process.
  ///
  /// Local matrices must be square.  Each process has exactly one
  /// matrix.  Those matrices may vary in dimensions.
  virtual size_t getNumRows() const;

  //! Whether the container has been successfully initialized.
  virtual bool isInitialized() const;

  //! Whether the container has been successfully computed.
  virtual bool isComputed() const;

  //! Set all necessary parameters.
  virtual void setParameters(const Teuchos::ParameterList& List);

  //@}
  //! \name Mathematical functions
  //@{

  //! Do all set-up operations that only require matrix structure.
  virtual void initialize ();

  //! Extract the local diagonal block and prepare the solver.
  virtual void compute ();

  //! Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  virtual void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode=Teuchos::NO_TRANS,
         scalar_type alpha=Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta=Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                 Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                 const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& D,
                 Teuchos::ETransp mode=Teuchos::NO_TRANS,
                 scalar_type alpha=Teuchos::ScalarTraits<scalar_type>::one(),
                 scalar_type beta=Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //@}
  //! \name Miscellaneous methods
  //@{

  /// \brief Print information about this object to the given output stream.
  ///
  /// operator<< uses this method.
  virtual std::ostream& print (std::ostream& os) const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  virtual std::string description () const;

  //! Print the object with some verbosity level to the given FancyOStream.
  virtual void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

  //@}
private:
  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  TriDiContainer (const TriDiContainer<MatrixType, LocalScalarType>& rhs);

  //! Extract the submatrix identified by the local indices set by the constructor.
  void extract (const Teuchos::RCP<const row_matrix_type>& globalMatrix);

  /// \brief Factor the extracted submatrix.
  ///
  /// Call this after calling extract().
  void factor ();

  typedef Tpetra::MultiVector<local_scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> local_mv_type;

  /// \brief Post-permutation, post-view version of apply().
  ///
  /// apply() first does any necessary subset permutation and view
  /// creation (or copying data), then calls this method to solve the
  /// linear system with the diagonal block.
  ///
  /// \param X [in] Subset permutation of the input X of apply().
  /// \param Y [in] Subset permutation of the input/output Y of apply().
  void
  applyImpl (const local_mv_type& X,
             local_mv_type& Y,
             Teuchos::ETransp mode,
             const local_scalar_type alpha,
             const local_scalar_type beta) const;

  //! Number of rows in the local matrix.
  size_t numRows_;

  //! The local diagonal block, which compute() extracts.
  Teuchos::SerialTriDiMatrix<int, local_scalar_type> diagBlock_;

  //! Permutation array from LAPACK (GETRF).
  Teuchos::Array<int> ipiv_;

  //! Map of the input and output Tpetra::MultiVector arguments of applyImpl().
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> > localMap_;

  //! Solution vector.
  mutable Teuchos::RCP<local_mv_type> Y_;

  //! Input vector for local problems
  mutable Teuchos::RCP<local_mv_type> X_;

  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;

  //! If \c true, the container has been successfully computed.
  bool IsComputed_;
};

}// namespace Ifpack2

#endif // IFPACK2_TRIDICONTAINER_DECL_HPP
