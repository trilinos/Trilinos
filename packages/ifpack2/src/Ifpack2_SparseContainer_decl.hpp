// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "Ifpack2_ILUT_decl.hpp"
#include <vector>
#ifdef HAVE_IFPACK2_AMESOS2
#include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif

namespace Ifpack2 {

/// \class SparseContainer
/// \brief Store and solve a local sparse linear problem.
/// \tparam A specialization of Tpetra::RowMatrix.
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
/// This class stores each block as a sparse matrix.  Using a sparse
/// matrix for each block is a good idea when the blocks are large and
/// sparse.  For small and / or dense blocks, it would probably be
/// better to use an implementation of Container that stores the
/// blocks densely, like DenseContainer.  You may also want to
/// consider BandedContainer.
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
/// <li> <tt>apply (const mv_type& X, mv_type& Y, ...)</tt>, where <tt>mv_type</tt>
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
/// These assumptions may be violated if the input matrix is a
/// Tpetra::CrsMatrix that was constructed with a user-provided column
/// Map.  The assumptions are not mathematically necessary and could
/// be relaxed at any time.  Implementers who wish to do so will need
/// to modify the extract() method, so that it translates explicitly
/// between local row and column indices, instead of just assuming
/// that they are the same.
template<typename MatrixType, typename InverseType>
class SparseContainer
: public ContainerImpl<MatrixType, typename InverseType::scalar_type>
{

  //! @name Internal type aliases (private)
  //@{
private:
  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization or a
  /// Tpetra::CrsMatrix specialization.  It may have entirely
  /// different template parameters (e.g., \c scalar_type) than
  /// <tt>InverseType</tt>.
  using matrix_type = MatrixType;
  /// \brief The second template parameter of this class.
  ///
  /// This must be a specialization of Ifpack2::Preconditioner or one
  /// of its subclasses.  It may have entirely different template
  /// parameters (e.g., \c scalar_type) than \c MatrixType.
  using inverse_type = InverseType;

  using typename Container<MatrixType>::SC;
  using typename Container<MatrixType>::LO;
  using typename Container<MatrixType>::GO;
  using typename Container<MatrixType>::NO;

  using typename Container<MatrixType>::mv_type;
  using typename Container<MatrixType>::map_type;
  using typename Container<MatrixType>::vector_type;
  using typename Container<MatrixType>::import_type;

  using InverseScalar = typename InverseType::scalar_type;
  using InverseLocalOrdinal = typename InverseType::local_ordinal_type;
  using InverseGlobalOrdinal = typename InverseType::global_ordinal_type;
  using InverseNode = typename InverseType::node_type;

  using typename ContainerImpl<MatrixType, InverseScalar>::block_crs_matrix_type;

  using inverse_mv_type = Tpetra::MultiVector<InverseScalar, InverseLocalOrdinal, InverseGlobalOrdinal, InverseNode>;
  using InverseCrs = Tpetra::CrsMatrix<InverseScalar, InverseLocalOrdinal, InverseGlobalOrdinal, InverseNode>;
  using InverseMap = typename Tpetra::Map<InverseLocalOrdinal, InverseGlobalOrdinal, InverseNode>;
  using InverseGraph = typename InverseCrs::crs_graph_type;
  using typename Container<MatrixType>::HostView;
  using typename Container<MatrixType>::ConstHostView;
  using HostViewInverse = typename inverse_mv_type::dual_view_type::t_host;

  static_assert(std::is_same<MatrixType,
                  Tpetra::RowMatrix<SC, LO, GO, NO>>::value, "Ifpack2::SparseContainer: Please use MatrixType = Tpetra::RowMatrix.");

  /// \brief The (base class) type of the input matrix.
  ///
  /// The input matrix to the constructor may be either a
  /// Tpetra::RowMatrix specialization or a Tpetra::CrsMatrix
  /// specialization.  However, we want to make the constructor as
  /// general as possible, so we always accept the matrix as a
  /// Tpetra::RowMatrix.  This type is the appropriate
  /// specialization of Tpetra::RowMatrix.
  using typename Container<MatrixType>::row_matrix_type;
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
  SparseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                   const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                   const Teuchos::RCP<const import_type>& importer,
                   bool pointIndexed);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~SparseContainer();

  //@}
  //! \name Get and set methods
  //@{

  //! Set all necessary parameters.
  virtual void setParameters(const Teuchos::ParameterList& List);

  //@}
  //! \name Mathematical functions
  //@{

  //! Do all set-up operations that only require matrix structure.
  virtual void initialize();

  //! Initialize and compute all blocks.
  virtual void compute ();

  //! Free all per-block resources: <tt>Inverses_</tt>, and <tt>diagBlocks_</tt>. 
  //! Called by \c BlockRelaxation when the input matrix is changed. Also calls
  //! \c Container::clearBlocks()
  void clearBlocks ();

  //! Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  virtual void
  apply (ConstHostView X,
         HostView Y,
         int blockIndex,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         SC alpha = Teuchos::ScalarTraits<SC>::one(),
         SC beta = Teuchos::ScalarTraits<SC>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply (ConstHostView X,
                 HostView Y,
                 ConstHostView W,
                 int blockIndex,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 SC alpha = Teuchos::ScalarTraits<SC>::one(),
                 SC beta = Teuchos::ScalarTraits<SC>::zero()) const;

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

  /// \brief Get the name of this container type for Details::constructContainer()
  static std::string getName();

private:

  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  SparseContainer (const SparseContainer<MatrixType,InverseType>& rhs);

  //! Extract the submatrices identified by the local indices set by the constructor.
  void extract ();
  void extractGraph ();
  void extractValues ();

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
  solveBlockMV(const inverse_mv_type& X,
               inverse_mv_type& Y,
               int blockIndex,
               Teuchos::ETransp mode,
               InverseScalar alpha,
               InverseScalar beta) const;

  //! The local diagonal block, which compute() extracts.
  std::vector<Teuchos::RCP<InverseCrs>> diagBlocks_;

  //! Scratch copy of X, used in solveBlock, # of rows is size of corresponding block
  mutable std::vector<inverse_mv_type> invX;
  //! Scratch copy of Y, used in solveBlock, # of rows is size of corresponding block
  mutable std::vector<inverse_mv_type> invY;

  /// \brief Local operators.
  ///
  /// InverseType must be a specialization of Ifpack2::Preconditioner,
  /// with the same template parameters (in the same order) as those
  /// of \c diagBlocks_ above.  Its apply() method defines the action
  /// of the inverse of the local matrix.  See the class documentation
  /// for more details.
  mutable Teuchos::Array<Teuchos::RCP<InverseType>> Inverses_;
  //! Serial communicator (containing only MPI_COMM_SELF if MPI is used).
  Teuchos::RCP<Teuchos::Comm<int>> localComm_;


  //! Parameters for the InverseType linear solve operator.
  Teuchos::ParameterList List_;
};

}// namespace Ifpack2

#endif // IFPACK2_SPARSECONTAINER_HPP
