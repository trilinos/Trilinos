// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKRELAXATION_DECL_HPP
#define IFPACK2_BLOCKRELAXATION_DECL_HPP

/// \file Ifpack2_BlockRelaxation_decl.hpp
/// \brief Ifpack2::BlockRelaxation class declaration

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_ContainerFactory.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include <type_traits>

namespace Ifpack2 {

/// \class BlockRelaxation
/// \brief Block relaxation preconditioners (or smoothers) for
///   Tpetra::RowMatrix and Tpetra::CrsMatrix sparse matrices.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
/// \tparam ContainerType DO NOT SPECIFY THIS EXPLICITLY.
///   This exists only for backwards compatibility.
///
/// This class implements the construction and application of block
/// relaxation preconditioners and smoothers, for sparse matrices
/// represented as Tpetra::RowMatrix or Tpetra::CrsMatrix.  This class
/// implements Tpetra::Operator, and its apply() method applies the
/// block relaxation.
///
/// BlockRelaxation implements block variants of the following
/// relaxations:
/// <ul>
/// <li> Damped Jacobi </li>
/// <li> Damped Gauss-Seidel, i.e., SOR </li>
/// <li> Damped symmetric Gauss-Seidel, i.e., symmetric SOR </li>
/// </ul>
///
/// For a list of supported parameters, please refer to the
/// documentation of setParameters().
template<class MatrixType, class ContainerType = Container<MatrixType> >
class BlockRelaxation :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
public:
  //! @name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  //! Node type of the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization corresponding to \c MatrixType.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::BlockRelaxation: Please use MatrixType = Tpetra::RowMatrix.");

  static_assert (std::is_same<ContainerType, Container<row_matrix_type> >::value,
                 "Ifpack2::BlockRelaxation: Do NOT specify the (second) "
                 "ContainerType template parameter explicitly.  The default "
                 "value is fine.  Please instead specify the container type to "
                 "use by setting the \"relaxation: container\" parameter.");

  //! Tpetra::Importer specialization for use with \c MatrixType and compatible MultiVectors.
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;

private:

  /// \brief Variant of setParameters() that takes a nonconst Teuchos::ParameterList.
  ///
  /// This variant fills in default values for any valid parameters
  /// that are not in the input list.
  void setParametersImpl(Teuchos::ParameterList& params);

  void computeImporter() const;

  //! \name Internal typedefs (handy for brevity and code clarity)
  //@{
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vector_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                                 global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::BlockCrsMatrix<scalar_type, local_ordinal_type,
                                 global_ordinal_type, node_type> block_crs_matrix_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
public:
  //@}
  //@}
  // \name Constructors and Destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param Matrix [in] The matrix for which to make the constructor.
  ///   Tpetra::RowMatrix is the base class of Tpetra::CrsMatrix, so
  ///   you may give either a Tpetra::RowMatrix or a Tpetra::CrsMatrix
  ///   here.
  ///
  /// The results of apply() are undefined if you change the sparse
  /// matrix after invoking this constructor, without first calling
  /// initialize() and compute() (in that order) to reinitialize the
  /// preconditioner.
  ///
  /// The "explicit" keyword just means that you must invoke the
  /// Relaxation constructor explicitly; you aren't allowed to use it
  /// as an implicit conversion ("cast").  For example, you may do
  /// this (namespaces and Tpetra template parameters omitted for
  /// brevity):
  /// \code
  /// RCP<const CrsMatrix<...> > A = ...;
  /// BlockRelaxation<RowMatrix<...> > R (A);
  /// \endcode
  /// but you may not do this:
  /// \code
  /// // Declaration of some user-defined function.
  /// void foo (const BlockRelaxation<RowMatrix<...> >& R);
  ///
  /// RCP<const CrsMatrix<...> > A = ...;
  /// foo (A);
  /// \endcode
  explicit BlockRelaxation (const Teuchos::RCP<const row_matrix_type>& Matrix);

  //! Destructor.
  virtual ~BlockRelaxation ();

  //@}
  //! \name Preconditioner computation methods
  //@{

  //! Sets all the parameters for the preconditioner
  /**
     Valid parameters are:
     <ul>
      <li> "relaxation: type"<br>
        Valid values (string):<br>
        <ul>
         <li> "Jacobi"
         <li> "Gauss-Seidel"
         <li> "Symmetric Gauss-Seidel"
        </ul>
      <li> "relaxation: sweeps" (int)
      <li> "relaxation: damping factor" (scalar)
      <li> "relaxation: zero starting solution" (bool)
      <li> "relaxation: backward mode" (bool)
      <li> "partitioner: type" <br>
        Valid values (string):<br>
        <ul>
         <li> "linear"
         <li> "line"
         <li> "user"
        </ul>
      <li> "partitioner: local parts" (local ordinal)
      <li> "partitioner: overlap" (int)
     </ul>

     \see Ifpack2::Details::UserPartitioner.
  */
  void setParameters(const Teuchos::ParameterList& params);

  bool supportsZeroStartingSolution() { return true; }

  void setZeroStartingSolution (bool zeroStartingSolution) { ZeroStartingSolution_ = zeroStartingSolution; };

  //! Return a list of all the parameters that this class accepts.
  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters () const;

  //! Initialize
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! compute the preconditioner for the specified matrix, diagonal perturbation thresholds and relaxation parameters.
  void compute();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}
  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param A [in] The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method with a matrix different than the current
  /// matrix resets the preconditioner's state.  After calling this
  /// method with a nonnull input, you must first call initialize()
  /// and compute() (in that order) before you may call apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}
  //! @name Methods implementing the Tpetra::Operator interface.
  //@{

  //! Applies the preconditioner to X, returns the result in Y.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO compliant.
  */
  void apply(const MV& X,
             MV& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
             scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const;

  bool hasTransposeApply() const;

  //! Applies the matrix to a Tpetra::MultiVector.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void applyMat(const MV& X,
                MV& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Attribute accessor methods
  //@{

  //! The communicator over which the input matrix is distributed.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The input matrix of this preconditioner's constructor.
  Teuchos::RCP<const row_matrix_type> getMatrix() const;

  //! Returns the number of calls to initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in initialize().
  double getInitializeTime() const;

  //! Returns the time spent in compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;  

  //@}
  //! @name Implementation of the Teuchos::Describable interface
  //@{

  //! A one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

  //@}

  //! For diagnostic purposes
  Teuchos::RCP<Ifpack2::Partitioner<Tpetra::RowGraph<local_ordinal_type,global_ordinal_type,node_type> > > getPartitioner(){return Partitioner_;}

private:

  //! Copy constructor; do not use (declared but unimplemented)
  BlockRelaxation (const BlockRelaxation<MatrixType, ContainerType> & RHS);

  //! Assignment operator; do not use (declared but unimplemented)
  BlockRelaxation<MatrixType,ContainerType>&
  operator= (const BlockRelaxation<MatrixType, ContainerType>& RHS);

  virtual void ApplyInverseJacobi (const MV& X, MV& Y) const;

  virtual void ApplyInverseGS (const MV& X, MV& Y) const;

  virtual void ApplyInverseSGS (const MV& X, MV& Y) const;

  //! Initialize structural information in \c Container_ using
  //! <tt>Partitioner</tt>.
  void ExtractSubmatricesStructure();

  //@}
  //! \name Internal data and parameters
  //@{

  //! The sparse matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> A_;

  //! Contains the (block) diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Container<row_matrix_type> > Container_;

  // FIXME (mfh 06 Oct 2014) This doesn't comply with the naming
  // convention for instance members of a class.  Furthermore, the
  // class should keep the Vector, not the ArrayRCP to the data _in_
  // the Vector.
  // FIXED! (amk 10 Nov 2015)
  mutable Teuchos::RCP<vector_type> DiagRCP_;

  //! Contains information about non-overlapping partitions.
  Teuchos::RCP<Ifpack2::Partitioner<Tpetra::RowGraph<local_ordinal_type,global_ordinal_type,node_type> > > Partitioner_;

  //! Which partitioner class to use; is \c linear by default
  //! but can be specified as \c line or \c user in the parameter list.
  std::string PartitionerType_;

  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;

  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;

  //! Number of local blocks
  local_ordinal_type NumLocalBlocks_;

  //! How to solve each block; the "container type"
  std::string containerType_;

  //! Which type of point relaxation approach to use
  Details::RelaxationType PrecType_;

  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;

  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;

  //! True when the input matrix is a Tpetra::BlockCrsMatrix, which
  //! means that multiple DOFs can be counted as a single "row"
  bool hasBlockCrsMatrix_;

  //! Backward-Mode Gauss Seidel
  bool DoBackwardGS_;

  //! Number of rows of overlap between adjacent blocks.
  int OverlapLevel_;

  //! Combine overlapping soltuions in nonsymmetric way
  //only valid with block Jacobi relaxation and overlapping blocks (e.g., defined by "partitioner: type" "user" and "parts: " or "global ID parts:"). Average solutions in overlapped regions (i.e., after summing different solutions divide by number of blocks contain this dof). When false (the default) symmetric averaging performed (i.e., average residuals and solutions).
  bool nonsymCombine_;

  //! Corresponds to "schwarz: combine mode" to distinguish between RAS and ADD when computing weights
  std::string schwarzCombineMode_;

  //! Damping factor.
  scalar_type DampingFactor_;

  //! Whether to decouple DOFs
  bool decouple_;

  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;

  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;

  //! Contains the number of successful calls to initialize().
  int NumInitialize_;

  //! Contains the number of successful call to compute().
  int NumCompute_;

  //! Whether apply should register and use a timer
  bool TimerForApply_;

  //! Contains the number of successful call to apply().
  mutable int NumApply_;

  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;

  //! Contains the time for all successful calls to compute().
  double ComputeTime_;

  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;

  //! The number of input matrix rows on the local process.
  local_ordinal_type NumLocalRows_;

  //! The total number of rows in the input matrix.
  global_ordinal_type NumGlobalRows_;

  //! The total number of nonzero entries in the input matrix.
  global_ordinal_type NumGlobalNonzeros_;

  //! Weight array used in \c Container::weightedApply(). If it is required,
  //! \c W_ is created with the row map of <tt>A_</tt>.
  Teuchos::RCP<vector_type> W_;

  mutable Teuchos::RCP<const Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>> Importer_;

  //@}
}; //class BlockRelaxation

}//namespace Ifpack2

#endif // IFPACK2_BLOCKRELAXATION_DECL_HPP

