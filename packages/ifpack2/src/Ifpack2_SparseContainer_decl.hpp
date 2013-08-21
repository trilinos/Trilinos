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
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

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
/// This class stores each block as a sparse matrix.  Thus, \c
/// MatrixType must be a specialization of Tpetra::RowMatrix or of its
/// subclass Tpetra::CrsMatrix.  Using a sparse matrix for each block
/// is a good idea when the blocks are large and sparse.  For small
/// and / or dense blocks, it would probably be better to use an
/// implementation of Container that stores the blocks densely.
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
/// \warning Please don't rely too much on this interface, because the
///   interface needs to be reworked to make it more rational.
template<typename MatrixType, typename InverseType >
class SparseContainer : public Container<MatrixType,InverseType> {
public:
  typedef typename MatrixType::scalar_type          MatrixScalar;
  typedef typename MatrixType::local_ordinal_type   MatrixLocalOrdinal;
  typedef typename MatrixType::global_ordinal_type  MatrixGlobalOrdinal;
  typedef typename MatrixType::node_type            MatrixNode;

  typedef typename InverseType::scalar_type         InverseScalar;
  typedef typename InverseType::local_ordinal_type  InverseLocalOrdinal;
  typedef typename InverseType::global_ordinal_type InverseGlobalOrdinal;
  typedef typename InverseType::node_type           InverseNode;

  //! \name Constructor and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param NumRows [in] Number of rows in the local matrix on each
  ///   process.  This may be different on different processes.
  ///
  /// \param NumVectors [in] Hint for the number of columns to expect
  ///   in the \c X multivector input argument of apply().
  ///
  /// \warning FIXME (mfh 20 Aug 2013) It should not be necessary to
  ///   ask for the number of rows in the local matrix on each
  ///   process, if the local matrix is a Tpetra::RowMatrix
  ///   specialization.  We do need to revise this interface.
  SparseContainer (const size_t NumRows, const size_t NumVectors = 1);

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

  //! Returns the ID associated to local row i.
  /*!
   * The set of (local) rows assigned to this container is defined
   * by calling ID(i) = j, where i (from 0 to NumRows()) indicates
   * the container-row, and j indicates the local row in the calling
   * process.
   *
   * This is usually used to recorder the local row ID (on calling process)
   * of the i-th row in the container.
   */
  virtual MatrixLocalOrdinal & ID(const size_t i);

  //! Whether the container has been successfully initialized.
  virtual bool isInitialized() const;

  //! Whether the container has been successfully computed.
  virtual bool isComputed() const;

  //! Set all necessary parameters.
  virtual void setParameters(const Teuchos::ParameterList& List);

  //@}
  //! \name Mathematical functions
  //@{

  /*!
   * \brief Initializes the container, by completing all the operations based
   * on matrix structure.
   *
   * \note After a call to Initialize(), no new matrix entries can be
   * added.
   */
  virtual void initialize();

  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual void
  compute (const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix);

  /// \brief Computes Y = alpha * M^{-1} X + beta*Y
  ///
  /// Here X and Y are the size of the global problem the container
  /// was extracted from to begin with.
  virtual void
  apply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
         Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
         Teuchos::ETransp mode=Teuchos::NO_TRANS,
         MatrixScalar alpha=Teuchos::ScalarTraits< MatrixScalar >::one(),
         MatrixScalar beta=Teuchos::ScalarTraits< MatrixScalar >::zero());

  /// \brief Computes Y = alpha * diag(D) * M^{-1} (diag(D) X) + beta*Y
  ///
  /// Here D, X and Y are the size of the global problem the container
  /// was extracted from to begin with.  D has to be
  /// a single Vector, while X and Y can be MultiVectors.  This
  /// function is designed to support techniques with overlap,
  /// such as Schwarz methods.
  virtual void
  weightedApply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
                 Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
                 const Tpetra::Vector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& D,
                 Teuchos::ETransp mode=Teuchos::NO_TRANS,
                 MatrixScalar alpha=Teuchos::ScalarTraits< MatrixScalar >::one(),
                 MatrixScalar beta=Teuchos::ScalarTraits< MatrixScalar >::zero());

  //@}
  //! \name Miscellaneous methods
  //@{

  //! Destroy all data.
  virtual void destroy();

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& print(std::ostream& os) const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  virtual std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  //! Copy constructor: Declared but not implemented, to forbid copy construction.
  SparseContainer(const SparseContainer<MatrixType,InverseType>& rhs);

  //! Extract the submatrices identified by the ID set int ID().
  virtual void extract(const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix);

  //! Number of rows in the local matrix.
  size_t numRows_;
  //! Linear map on which the local matrix is based.
  Teuchos::RCP<Tpetra::Map<InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > Map_;
  //! Pointer to the local matrix.
  Teuchos::RCP<Tpetra::CrsMatrix<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > Matrix_;
  //! Solution vector.
  mutable Teuchos::RCP<Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > Y_;
  //! Input vector for local problems
  mutable Teuchos::RCP<Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> > X_;
  //! Contains the subrows/subcols of A that will be inserted in Matrix_.
  std::vector<MatrixLocalOrdinal> GID_;
  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;
  //! Serial communicator (containing only MPI_COMM_SELF if MPI is used).
  Teuchos::RCP<Teuchos::Comm<int> > LocalComm_;
  //! Pointer to an Ifpack2_Preconditioner object whose ApplyInverse() defined the action of the inverse of the local matrix.
  Teuchos::RCP<InverseType> Inverse_;

  Teuchos::ParameterList List_;

  //! Whether apply() and weightedApply() need to permute the input X and output Y.
  bool needPermutation_;
};

}// namespace Ifpack2



#endif // IFPACK2_SPARSECONTAINER_HPP
