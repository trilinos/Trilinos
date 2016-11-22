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
#ifndef IFPACK2_BANDEDCONTAINER_DECL_HPP
#define IFPACK2_BANDEDCONTAINER_DECL_HPP

/// \file Ifpack2_BandedContainer_decl.hpp
/// \brief Ifpack2::BandedContainer class declaration

#include "Ifpack2_Container.hpp"
#include "Ifpack2_Details_MultiVectorLocalGatherScatter.hpp"
#include "Ifpack2_Details_LapackSupportsScalar.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialBandDenseMatrix.hpp"

namespace Ifpack2 {

/// \class BandedContainer
/// \brief Store and solve a local Banded linear problem.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// Please refer to the documentation of the Container
/// interface. Currently, Containers are used by BlockRelaxation.
/// Block relaxations need to be able to do two things:
/// <ol>
/// <li> Store the diagonal blocks </li>
/// <li> Solve linear systems with each diagonal block </li>
/// </ol>
/// BandedContainer stores the diagonal blocks as Banded matrices, and
/// solves them using either LAPACK (for the four Scalar types that it
/// supports) or a custom LU factorization (for Scalar types not
/// supported by LAPACK).
///
/// As with Ifpack2::Container, <tt>MatrixType</tt> must be a
/// specialization of Tpetra::RowMatrix.  Using a Banded matrix for
/// each block is a good idea when the blocks are small.  For large
/// and / or sparse blocks, it would probably be better to use an
/// implementation of Container that stores the blocks sparsely, in
/// particular SparseContainer.
///
/// This class may store the Banded local matrix using values of a
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
/// These assumptions may be violated if the input matrix is a
/// Tpetra::CrsMatrix that was constructed with a user-provided column
/// Map.  The assumptions are not mathematically necessary and could
/// be relaxed at any time.  Implementers who wish to do so will need
/// to modify the extract() method, so that it translates explicitly
/// between local row and column indices, instead of just assuming
/// that they are the same.
template<class MatrixType,
         class LocalScalarType,
         bool supportsLapackScalar = ::Ifpack2::Details::LapackSupportsScalar<LocalScalarType>::value>
class BandedContainer;

template<class MatrixType, class LocalScalarType>
class BandedContainer<MatrixType, LocalScalarType, true> :
    public Container<MatrixType> {
  //! @name Internal typedefs (private)
  //@{
private:
  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization or a
  /// Tpetra::CrsMatrix specialization.  It may have entirely
  /// different template parameters (e.g., \c scalar_type) than
  /// <tt>InverseType</tt>.
  typedef MatrixType matrix_type;
  //! The second template parameter of this class.
  typedef LocalScalarType local_scalar_type;
  //! The internal representation of LocalScalarType in Kokkos::View
  typedef typename Kokkos::Details::ArithTraits<local_scalar_type>::val_type local_impl_scalar_type;

  //! The type of entries in the input (global) matrix.
  typedef typename Container<MatrixType>::scalar_type scalar_type;
  //! The type of local indices in the input (global) matrix.
  typedef typename Container<MatrixType>::local_ordinal_type local_ordinal_type;
  //! The type of global indices in the input (global) matrix.
  typedef typename Container<MatrixType>::global_ordinal_type global_ordinal_type;
  //! The Node type of the input (global) matrix.
  typedef typename Container<MatrixType>::node_type node_type;

  typedef typename Container<MatrixType>::mv_type mv_type;
  typedef typename Container<MatrixType>::map_type map_type;
  typedef Tpetra::MultiVector<local_scalar_type, local_ordinal_type, global_ordinal_type, node_type> local_mv_type;
  typedef typename Container<MatrixType>::vector_type vector_type;
  typedef typename Container<MatrixType>::partitioner_type partitioner_type;
  typedef typename Container<MatrixType>::import_type import_type;

  typedef typename Container<MatrixType>::HostView HostView;
  typedef typename local_mv_type::dual_view_type::t_host HostViewLocal;

  static_assert(std::is_same<MatrixType,
                  Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> >::value,
                "Ifpack2::BandedContainer: Please use MatrixType = Tpetra::RowMatrix.");

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
public:
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
  ///   <tt>number of subdiagonals
  ///   <tt>number of superdiagonals. Note: Internally, we store a Teuchos::SerialBandedMatrix
  ///       with kl+ku superdiagonals, as we need the addtional storage for the LU decomposition.
  BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                   const Teuchos::RCP<const import_type>& importer,
                   int OverlapLevel,
                   scalar_type DampingFactor);

  BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::Array<local_ordinal_type>& localRows);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~BandedContainer ();

  //@}
  //! \name Get and set methods
  //@{

  //! Whether the container has been successfully initialized.
  virtual bool isInitialized () const {
    return IsInitialized_;
  }

  //! Whether the container has been successfully computed.
  virtual bool isComputed () const {
    return IsComputed_;
  }

  //! Set all necessary parameters.
  virtual void setParameters (const Teuchos::ParameterList& List);

  //@}
  //! \name Mathematical functions
  //@{

  //! Do all set-up operations that only require matrix structure.
  virtual void initialize ();

  //! Initialize and compute each block.
  virtual void compute ();

  void clearBlocks();

  //! Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  virtual void
  apply (HostView& X,
         HostView& Y,
         int blockIndex,
         int stride,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply (HostView& X,
                 HostView& Y,
                 HostView& D,
                 int blockIndex,
                 int stride,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                 scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

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

  /// \brief Get the name of this container type for Details::constructContainer()
  static std::string getName();

private:
  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  BandedContainer (const BandedContainer<MatrixType, LocalScalarType>& rhs);

  //! Extract the submatrix identified by the local indices set by the constructor.
  void extract ();

  /// \brief Factor the extracted submatrix.
  ///
  /// Call this after calling extract().
  void factor ();

  /// \brief Post-permutation, post-view version of apply().
  ///
  /// apply() first does any necessary subset permutation and view
  /// creation (or copying data), then calls this method to solve the
  /// linear system with the diagonal block.
  ///
  /// \param X [in] Subset permutation of the input X of apply().
  /// \param Y [in] Subset permutation of the input/output Y of apply().
  void
  applyImpl (HostViewLocal& X,
             HostViewLocal& Y,
             int blockIndex,
             int stride,
             Teuchos::ETransp mode,
             const local_scalar_type alpha,
             const local_scalar_type beta) const;

  //! The local diagonal block, which compute() extracts.
  std::vector<Teuchos::SerialBandDenseMatrix<int, local_scalar_type> > diagBlocks_;

  //! Temporary X vector used in apply().
  mutable std::vector<HostViewLocal> X_local;

  //! Temporary Y vector used in apply().
  mutable std::vector<HostViewLocal> Y_local;

  //! Permutation array from LAPACK (GETRF).
  Teuchos::Array<int> ipiv_;

  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;

  //! If \c true, the container has been successfully computed.
  bool IsComputed_;

  Teuchos::Array<local_ordinal_type> kl_; //< number of subdiagonals
  Teuchos::Array<local_ordinal_type> ku_; //< number of superdiagonals

  //! Scalar data for all blocks
  local_scalar_type* scalars_;

  //! Offsets in scalars_ array for all blocks
  Teuchos::Array<local_ordinal_type> scalarOffsets_;
};

template<class MatrixType, class LocalScalarType>
class BandedContainer<MatrixType, LocalScalarType, false> :
    public Container<MatrixType> {
  //! @name Internal typedefs (private)
  //@{
private:
  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization or a
  /// Tpetra::CrsMatrix specialization.  It may have entirely
  /// different template parameters (e.g., \c scalar_type) than
  /// <tt>InverseType</tt>.
  typedef MatrixType matrix_type;
  //! The second template parameter of this class.
  typedef LocalScalarType local_scalar_type;
  //! The internal representation of LocalScalarType in Kokkos::View
  typedef typename Kokkos::Details::ArithTraits<local_scalar_type>::val_type local_impl_scalar_type;

  //! The type of entries in the input (global) matrix.
  typedef typename Container<MatrixType>::scalar_type scalar_type;
  //! The type of local indices in the input (global) matrix.
  typedef typename Container<MatrixType>::local_ordinal_type local_ordinal_type;
  //! The type of global indices in the input (global) matrix.
  typedef typename Container<MatrixType>::global_ordinal_type global_ordinal_type;
  //! The Node type of the input (global) matrix.
  typedef typename Container<MatrixType>::node_type node_type;

  typedef typename Container<MatrixType>::mv_type mv_type;
  typedef typename Container<MatrixType>::map_type map_type;
  typedef Tpetra::MultiVector<local_scalar_type, local_ordinal_type, global_ordinal_type, node_type> local_mv_type;
  typedef typename Container<MatrixType>::vector_type vector_type;
  typedef typename Container<MatrixType>::partitioner_type partitioner_type;
  typedef typename Container<MatrixType>::import_type import_type;

  typedef typename Container<MatrixType>::HostView HostView;
  typedef typename local_mv_type::dual_view_type::t_host HostViewLocal;

  static_assert(std::is_same<MatrixType,
                  Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> >::value,
                "Ifpack2::BandedContainer: Please use MatrixType = Tpetra::RowMatrix.");

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
public:
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
  ///   <tt>number of subdiagonals
  ///   <tt>number of superdiagonals. Note: Internally, we store a Teuchos::SerialBandedMatrix
  ///       with kl+ku superdiagonals, as we need the addtional storage for the LU decomposition.
  BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                   const Teuchos::RCP<const import_type>& importer,
                   int OverlapLevel,
                   scalar_type DampingFactor);

  BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::Array<local_ordinal_type>& localRows);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~BandedContainer ();

  //@}
  //! \name Get and set methods
  //@{

  //! Whether the container has been successfully initialized.
  virtual bool isInitialized () const {
    return IsInitialized_;
  }

  //! Whether the container has been successfully computed.
  virtual bool isComputed () const {
    return IsComputed_;
  }

  //! Set all necessary parameters.
  virtual void setParameters (const Teuchos::ParameterList& List);

  //@}
  //! \name Mathematical functions
  //@{

  //! Do all set-up operations that only require matrix structure.
  virtual void initialize ();

  //! Initialize and compute each block.
  virtual void compute ();

  void clearBlocks();

  //! Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  virtual void
  apply (HostView& X,
         HostView& Y,
         int blockIndex,
         int stride,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply (HostView& X,
                 HostView& Y,
                 HostView& D,
                 int blockIndex,
                 int stride,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                 scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

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

  /// \brief Get the name of this container type for Details::constructContainer()
  static std::string getName();

private:
  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  BandedContainer (const BandedContainer<MatrixType, LocalScalarType>& rhs);

  //! Extract the submatrix identified by the local indices set by the constructor.
  void extract ();

  /// \brief Factor the extracted submatrix.
  ///
  /// Call this after calling extract().
  void factor ();

  /// \brief Post-permutation, post-view version of apply().
  ///
  /// apply() first does any necessary subset permutation and view
  /// creation (or copying data), then calls this method to solve the
  /// linear system with the diagonal block.
  ///
  /// \param X [in] Subset permutation of the input X of apply().
  /// \param Y [in] Subset permutation of the input/output Y of apply().
  void
  applyImpl (HostViewLocal& X,
             HostViewLocal& Y,
             int blockIndex,
             int stride,
             Teuchos::ETransp mode,
             const local_scalar_type alpha,
             const local_scalar_type beta) const;

  //! The local diagonal block, which compute() extracts.
  std::vector<Teuchos::SerialBandDenseMatrix<int, local_scalar_type> > diagBlocks_;

  //! Temporary X vector used in apply().
  mutable std::vector<HostViewLocal> X_local;

  //! Temporary Y vector used in apply().
  mutable std::vector<HostViewLocal> Y_local;

  //! Permutation array from LAPACK (GETRF).
  Teuchos::Array<int> ipiv_;

  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;

  //! If \c true, the container has been successfully computed.
  bool IsComputed_;

  Teuchos::Array<local_ordinal_type> kl_; //< number of subdiagonals
  Teuchos::Array<local_ordinal_type> ku_; //< number of superdiagonals

  //! Scalar data for all blocks
  local_scalar_type* scalars_;

  //! Offsets in scalars_ array for all blocks
  Teuchos::Array<local_ordinal_type> scalarOffsets_;
};

}// namespace Ifpack2

#endif // IFPACK2_BANDEDCONTAINER_DECL_HPP
