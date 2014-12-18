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

#ifndef IFPACK2_SPARSECONTAINER_DECL_HPP
#define IFPACK2_SPARSECONTAINER_DECL_HPP

/// \file Ifpack2_SparseContainer_decl.hpp
/// \brief Ifpack2::SparseContainer class declaration

#include "Ifpack2_Container.hpp"
#include "Ifpack2_Details_MultiVectorLocalGatherScatter.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Ifpack2 {

/// \class SparseContainer
/// \brief Store and solve a local sparse linear problem.
///
/// Please refer to the documentation of the Container
/// interface. Currently, Containers are used by BlockRelaxation.
/// Block relaxations need to be able to do two things:
/// <ol>
/// <li> Store the diagonal blocks </li>
/// <li> Solve linear systems with each diagonal block </li>
/// </ol>
/// These correspond to the two template parameters:
/// <ol>
/// <li> \c MatrixType, which stores a sparse matrix </li>
/// <li> \c InverseType, which solves linear systems with that matrix </li>
/// </ol>
/// This class stores each block as a sparse matrix.  Thus,
/// <tt>MatrixType</tt> must be a specialization of Tpetra::RowMatrix
/// or of its subclass Tpetra::CrsMatrix.  Using a sparse matrix for
/// each block is a good idea when the blocks are large and sparse.
/// For small and / or dense blocks, it would probably be better to
/// use an implementation of Container that stores the blocks densely.
///
/// The \c InverseType template parameter represents the class to use
/// for solving linear systems with a block.  In SparseContainer, this
/// template parameter must be a specialization of Preconditioner.
/// Specifically, \c InverseType must implement the following methods:
/// <ul>
/// <li> A constructor that takes an <tt>RCP<const MatrixType></tt> </li>
/// <li> <tt>setParameters(Teuchos::ParameterList&)</tt> </li>
/// <li> <tt>initialize()</tt> </li>
/// <li> <tt>compute()</tt> </li>
/// <li> <tt>apply (const MV& X, MV& Y, ...)</tt>, where <tt>MV</tt>
///      is the appropriate specialization of Tpetra::MultiVector </li>
/// </ul>
/// We also assume that \c InverseType has the following typedefs:
/// <ul>
/// <li> \c scalar_type </li>
/// <li> \c local_ordinal_type </li>
/// <li> \c global_ordinal_type </li>
/// <li> \c node_type </li>
/// </ul>
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
///
/// SparseContainer currently assumes the following about the column
/// and row Maps of the input matrix:
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
template<typename MatrixType, typename InverseType>
class SparseContainer : public Container<MatrixType> {
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
  /// \brief The second template parameter of this class.
  ///
  /// This must be a specialization of Ifpack2::Preconditioner or one
  /// of its subclasses.  It may have entirely different template
  /// parameters (e.g., \c scalar_type) than \c MatrixType.
  typedef InverseType inverse_type;

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
  SparseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::ArrayView<const local_ordinal_type>& localRows);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~SparseContainer();

  //@}
  //! \name Get and set methods
  //@{

  /// \brief The number of rows in the local matrix on the calling process.
  ///
  /// Local matrices must be square.  Each process has exactly one matrix.
  /// Those matrices may vary in dimensions.
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
  virtual void initialize();

  //! Extract the local diagonal block and prepare the solver.
  virtual void compute ();

  //! Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  virtual void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode=Teuchos::NO_TRANS,
         scalar_type alpha=Teuchos::ScalarTraits< scalar_type >::one(),
         scalar_type beta=Teuchos::ScalarTraits< scalar_type >::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                 Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                 const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& D,
                 Teuchos::ETransp mode=Teuchos::NO_TRANS,
                 scalar_type alpha=Teuchos::ScalarTraits< scalar_type >::one(),
                 scalar_type beta=Teuchos::ScalarTraits< scalar_type >::zero()) const;

  //@}
  //! \name Miscellaneous methods
  //@{

  /// \brief Print information about this object to the given output stream.
  ///
  /// operator<< uses this method.
  virtual std::ostream& print(std::ostream& os) const;

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
  typedef typename InverseType::scalar_type         InverseScalar;
  typedef typename InverseType::local_ordinal_type  InverseLocalOrdinal;
  typedef typename InverseType::global_ordinal_type InverseGlobalOrdinal;
  typedef typename InverseType::node_type           InverseNode;

  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  SparseContainer (const SparseContainer<MatrixType,InverseType>& rhs);

  //! Extract the submatrices identified by the local indices set by the constructor.
  void extract (const Teuchos::RCP<const row_matrix_type>& globalMatrix);

  /// \brief Post-permutation, post-view version of apply().
  ///
  /// apply() first does any necessary subset permutation and view
  /// creation (or copying data), then calls this method to solve the
  /// linear system with the diagonal block.
  ///
  /// \param X [in] Subset permutation of the input X of apply(),
  ///   suitable for the first argument of Inverse_->apply().
  ///
  /// \param Y [in] Subset permutation of the input/output Y of apply(),
  ///   suitable for the second argument of Inverse_->apply().
  void
  applyImpl (const Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>& X,
             Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>& Y,
             Teuchos::ETransp mode,
             InverseScalar alpha,
             InverseScalar beta) const;

  //! Number of rows in the local matrix.
  size_t numRows_;
  //! Row Map of the local diagonal block.
  Teuchos::RCP<Tpetra::Map<InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > localMap_;
  //! The local diagonal block, which compute() extracts.
  Teuchos::RCP<Tpetra::CrsMatrix<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > diagBlock_;
  //! Solution vector.
  mutable Teuchos::RCP<Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > Y_;
  //! Input vector for local problems
  mutable Teuchos::RCP<Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > X_;
  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;
  //! Serial communicator (containing only MPI_COMM_SELF if MPI is used).
  Teuchos::RCP<Teuchos::Comm<int> > localComm_;

  /// \brief Local operator.
  ///
  /// InverseType must be a specialization of Ifpack2::Preconditioner,
  /// with the same template parameters (in the same order) as those
  /// of \c diagBlock_ above.  Its apply() method defines the action
  /// of the inverse of the local matrix.  See the class documentation
  /// for more details.
  Teuchos::RCP<InverseType> Inverse_;

  //! Parameters for the InverseType linear solve operator.
  Teuchos::ParameterList List_;
};

}// namespace Ifpack2



#endif // IFPACK2_SPARSECONTAINER_HPP
