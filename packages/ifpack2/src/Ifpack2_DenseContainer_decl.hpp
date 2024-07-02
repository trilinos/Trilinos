// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DENSECONTAINER_DECL_HPP
#define IFPACK2_DENSECONTAINER_DECL_HPP

/// \file Ifpack2_DenseContainer_decl.hpp
/// \brief Ifpack2::DenseContainer class declaration

#include "Ifpack2_Container.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_decl.hpp"

namespace Ifpack2 {

/// \class DenseContainer
/// \brief Store and solve a local dense linear problem.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// Please refer to the documentation of the Container
/// interface. Currently, Containers are used by BlockRelaxation.
/// Block relaxations need to be able to do two things:
/// <ol>
/// <li> Store the diagonal blocks </li>
/// <li> Solve linear systems with each diagonal block </li>
/// </ol>
/// DenseContainer stores the diagonal blocks as dense matrices, and
/// solves them using LAPACK.
///
/// As with Ifpack2::Container, \c MatrixType must be a specialization
/// of Tpetra::RowMatrix.  Using a dense matrix for each block is a
/// good idea when the blocks are small.  For large and / or sparse
/// blocks, it would probably be better to use an implementation of
/// Container that stores the blocks sparsely, in particular
/// SparseContainer.
///
/// This class may store the dense local matrix using values of a
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

template<class MatrixType, class LocalScalarType>
class DenseContainer
: public ContainerImpl<MatrixType, LocalScalarType>
{
  //! @name Internal typedefs (private)
  //@{
private:
  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization or a
  /// Tpetra::CrsMatrix specialization.  It may have entirely
  /// different template parameters (e.g., \c scalar_type) than
  /// <tt>InverseType</tt>.
  using matrix_type = MatrixType;
  //! The second template parameter of this class.
  using LSC = LocalScalarType;

  //! The type of entries in the input (global) matrix.
  using typename Container<MatrixType>::SC;
  //! The type of local indices in the input matrix.
  using typename Container<MatrixType>::LO;
  //! The type of global indices in the input matrix.
  using typename Container<MatrixType>::GO;
  //! The Node type of the input matrix.
  using typename Container<MatrixType>::NO;

  using typename Container<MatrixType>::mv_type;
  using local_mv_type = Tpetra::MultiVector<LSC, LO, GO, NO>;
  using typename Container<MatrixType>::map_type;
  using typename Container<MatrixType>::vector_type;
  using typename Container<MatrixType>::import_type;

  using typename Container<MatrixType>::ISC;
  using typename ContainerImpl<MatrixType, LocalScalarType>::LISC;
  using typename Container<MatrixType>::HostView;
  using typename ContainerImpl<MatrixType, LocalScalarType>::HostViewLocal;
  using typename ContainerImpl<MatrixType, LocalScalarType>::HostSubviewLocal;
  using typename ContainerImpl<MatrixType, LocalScalarType>::ConstHostSubviewLocal;

  static_assert(std::is_same<MatrixType, Tpetra::RowMatrix<SC, LO, GO, NO>>::value,
                "Ifpack2::DenseContainer: Please use MatrixType = Tpetra::RowMatrix.");

  /// \brief The (base class) type of the input matrix.
  ///
  /// The input matrix to the constructor may be either a
  /// Tpetra::RowMatrix specialization or a Tpetra::CrsMatrix
  /// specialization.  However, we want to make the constructor as
  /// general as possible, so we always accept the matrix as a
  /// Tpetra::RowMatrix.  This typedef is the appropriate
  /// specialization of Tpetra::RowMatrix.
  using typename Container<MatrixType>::row_matrix_type;

  using block_crs_matrix_type = Tpetra::BlockCrsMatrix<SC, LO, GO, NO>;
  //@}

public:
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor.
  ///
  /// \param matrix [in] The original input matrix.  This Container
  ///   will construct local diagonal blocks from its rows according to
  ///   <tt>partitions</tt>.
  /// \param partitioner [in] The Partitioner object that assigns
  ///   local rows of the input matrix to blocks.
  /// \param pointIndexed [in] If the input matrix is a \c Tpetra::BlockCrsMatrix,
  ///    whether elements of \c partitions[k] identify rows within blocks (true) or
  ///    whole blocks (false).
  DenseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                  const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                  const Teuchos::RCP<const import_type>& importer,
                  bool pointIndexed);

  DenseContainer (const DenseContainer<MatrixType, LocalScalarType>& rhs) = delete;

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~DenseContainer ();

  //@}
  //! \name Get and set methods
  //@{

  //@}
  //! \name Mathematical functions
  //@{

  //! Do all set-up operations that only require matrix structure.
  virtual void initialize ();

  //! Extract the local diagonal block and prepare the solver.
  virtual void compute ();

  void clearBlocks();

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
  //! Populate the diagonal blocks
  void extract();

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
  /// \param Y [in/out] Subset permutation of the input/output Y of apply().
  void
  solveBlock(ConstHostSubviewLocal X,
             HostSubviewLocal Y,
             int blockIndex,
             Teuchos::ETransp mode,
             const LSC alpha,
             const LSC beta) const;

  //! The local diagonal blocks, which compute() extracts.
  std::vector<Teuchos::SerialDenseMatrix<int, LSC>> diagBlocks_;

  //! Permutation array from LAPACK (GETRF).
  mutable Teuchos::Array<int> ipiv_;

  //! True if the input matrix type is a <tt>BlockCrsMatrix</tt>, false otherwise.
  bool hasBlockCrsMatrix_;

  //! Scalar array for all of the <tt>diagBlocks_</tt>.
  Teuchos::Array<LSC> scalars_;

  //! The indices in \c scalars_ where each block matrix starts.
  Teuchos::Array<GO> scalarOffsets_;
};

} // namespace Ifpack2

#endif // IFPACK2_DENSECONTAINER_DECL_HPP
