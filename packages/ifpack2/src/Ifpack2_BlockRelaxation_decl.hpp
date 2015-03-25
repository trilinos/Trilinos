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

#ifndef IFPACK2_BLOCKRELAXATION_DECL_HPP
#define IFPACK2_BLOCKRELAXATION_DECL_HPP

/// \file Ifpack2_BlockRelaxation_decl.hpp
/// \brief Ifpack2::BlockRelaxation class declaration

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Partitioner.hpp>
#include <Ifpack2_Details_CanChangeMatrix.hpp>
#include <Teuchos_Time.hpp>
#include <string>
#include <iostream>
#include <sstream>


namespace Ifpack2 {

/// \class BlockRelaxation
/// \brief Block relaxation preconditioners (or smoothers) for
///   Tpetra::RowMatrix and Tpetra::CrsMatrix sparse matrices.
/// \tparam MatrixType A specialization of Tpetra::CrsMatrix (better)
///   or Tpetra::RowMatrix (acceptable).
/// \tparam ContainerType A specialization or subclass of Container; a
///   type that knows how to solve linear systems with diagonal blocks
///   of MatrixType.  Those blocks may be either sparse or dense; the
///   subclass of Container controls the representation.
///
/// This class implements the construction and application of block
/// relaxation preconditioners and smoothers, for sparse matrices
/// represented as Tpetra::RowMatrix or Tpetra::CrsMatrix.  This class
/// implements Tpetra::Operator, and its apply() method applies the
/// block relaxation.
///
/// BlockRelaxation implements block variants of the following relaxations:
/// - (Damped) Jacobi;
/// - (Damped) Gauss-Seidel, i.e., SOR
/// - (Damped) symmetric Gauss-Seidel, i.e., symmetric SOR
///
/// For a list of supported parameters, please refer to the
/// documentation of setParameters().
template<class MatrixType, class ContainerType>
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
  /// BlockRelaxation<CrsMatrix<...> > R (A);
  /// \endcode
  /// but you may not do this:
  /// \code
  /// // Declaration of some user-defined function.
  /// void foo (const BlockRelaxation<CrsMatrix<...> >& R);
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
     Valid parameters include the following:
     <ul>
      <li> "relaxation: type"<br>
        Valid values (string):<br>
        <ul>
         <li> "Jacobi"
         <li> "Gauss-Seidel"
         <li> "Symmetric Gauss-Seidel"
        </ul>
      <li> "relaxation: sweeps" (int)
      <li> "relaxation: damping factor" (floating-point)
      <li> "relaxation: min diagonal value" (floating-point)
      <li> "relaxation: zero starting solution" (bool)
      <li> "relaxation: backward mode" (bool)
     </ul>
  */
  void setParameters(const Teuchos::ParameterList& params);

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
  void apply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
             Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
             scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getRangeMap() const;

  bool hasTransposeApply() const;

  //! Applies the matrix to a Tpetra::MultiVector.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void applyMat(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Attribute accessor methods
  //@{

  //! The communicator over which the input matrix is distributed.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The input matrix of this preconditioner's constructor.
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;

  //! Returns the number of flops in the computation phase.
  double getComputeFlops() const;

  //! Returns the number of flops for the application of the preconditioner.
  double getApplyFlops() const;

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

private:
  //! \name Internal typedefs (handy for brevity and code clarity)
  //@{
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  //@}

  //! Copy constructor; do not use (declared but unimplemented)
  BlockRelaxation (const BlockRelaxation<MatrixType, ContainerType> & RHS);

  //! Assignment operator; do not use (declared but unimplemented)
  BlockRelaxation<MatrixType,ContainerType>&
  operator= (const BlockRelaxation<MatrixType, ContainerType>& RHS);

  virtual void ApplyInverseJacobi (const MV& X, MV& Y) const;

  virtual void DoJacobi (const MV& X, MV& Y) const;

  virtual void ApplyInverseGS (const MV& X, MV& Y) const;

  virtual void DoGaussSeidel (MV& X, MV& Y) const;

  virtual void ApplyInverseSGS (const MV& X, MV& Y) const;

  virtual void DoSGS (MV& X, MV& Y) const;

  void ExtractSubmatrices ();

  //@}
  //! \name Internal data and parameters
  //@{

  //! The sparse matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> A_;

  //! Timer
  Teuchos::RCP<Teuchos::Time> Time_;

  //! Import object for parallel GS and SGS
  Teuchos::RCP<const Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> > Importer_;

  //! Weights for overlapping overlapped Jacobi only.
  Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > W_;

  // Level of overlap among blocks (for overlapped Jacobi only).
  int OverlapLevel_;

  //! Contains the (block) diagonal elements of \c Matrix.
  mutable std::vector<Teuchos::RCP<ContainerType> > Containers_;

  //  mutable Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Diagonal_;

  // FIXME (mfh 06 Oct 2014) This doesn't comply with the naming
  // convention for instance members of a class.  Furthermore, the
  // class should keep the Vector, not the ArrayRCP to the data _in_
  // the Vector.
  Teuchos::ArrayRCP< const scalar_type > DiagRCP;

  //! Contains information about non-overlapping partitions.
  Teuchos::RCP<Ifpack2::Partitioner<Tpetra::RowGraph<local_ordinal_type,global_ordinal_type,node_type> > > Partitioner_;

  std::string PartitionerType_;

  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;

  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;

  //! Number of local blocks
  local_ordinal_type NumLocalBlocks_;

  //! Which type of point relaxation approach to use
  Details::RelaxationType PrecType_;

  //! Damping factor.
  scalar_type DampingFactor_;

  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;

  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;

  //! Backward-Mode Gauss Seidel
  bool DoBackwardGS_;

  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;

  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;

  //! Contains the number of successful calls to initialize().
  int NumInitialize_;

  //! Contains the number of successful call to compute().
  int NumCompute_;

  //! Contains the number of successful call to apply().
  mutable int NumApply_;

  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;

  //! Contains the time for all successful calls to compute().
  double ComputeTime_;

  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;

  //! Contains the number of flops for compute().
  double ComputeFlops_;

  //! Contain sthe number of flops for apply().
  mutable double ApplyFlops_;

  //! Number of local rows.
  size_t NumMyRows_;

  //! Number of global rows.
  global_size_t NumGlobalRows_;

  //! Number of global nonzeros.
  global_size_t NumGlobalNonzeros_;
  //@}
}; //class BlockRelaxation

}//namespace Ifpack2

#endif // IFPACK2_BLOCKRELAXATION_DECL_HPP

