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

#ifndef IFPACK2_DIAGONAL_DECL_HPP
#define IFPACK2_DIAGONAL_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"

namespace Ifpack2 {

//! A class for diagonal preconditioning.
/**
Ifpack2::Diagonal wraps a vector as a diagonal preconditioner.
The preconditioner is simply defined by
\f[
z_i = D_i r_i,
\f]
where \f$r,z\f$ are the vector to be preconditioned and the preconditioned vector, and \f$D_i\f$ is the i-th element of the scaling vector.

When Ifpack2::Diagonal is constructed with a matrix, \f$D\f$ contains the inverted diagonal elements from the matrix.

When Ifpack2::Diagonal is constructed with a vector, \f$D\f$ is the caller-supplied vector.

\author Michael Heroux, ETHZ/D-INFK

\date Ifpack2 conversion (from Ifpack code) 31-Mar-2010
 */
template<class MatrixType>
class Diagonal :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type> {
public:
  typedef TEUCHOS_DEPRECATED typename MatrixType::scalar_type Scalar;
  typedef TEUCHOS_DEPRECATED typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef TEUCHOS_DEPRECATED typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef TEUCHOS_DEPRECATED typename MatrixType::node_type Node;
  typedef TEUCHOS_DEPRECATED typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType magnitudeType;

  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  //! Tpetra::CrsMatrix specialization used by this class.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

  /// \brief Constructor that takes a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  Diagonal (const Teuchos::RCP<const row_matrix_type>& A);

  /// \brief Constructor that takes a Tpetra::CrsMatrix.
  ///
  /// This constructor exists to avoid "ambiguous constructor"
  /// warnings.  It does the same thing as the constructor that takes
  /// a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  Diagonal (const Teuchos::RCP<const crs_matrix_type>& A_in);

  /// \brief Constructor that accepts a Tpetra::Vector of diagonal entries.
  ///
  /// If your compiler complains about this constructor being ambigous
  /// with the other constructor overload, instead call the
  /// free-standing function Ifpack2::createDiagonalPreconditioner
  /// which is located at the bottom of this header file.  (This issue
  /// may arise if this constructor is called with a
  /// <tt>Teuchos::RCP<Tpetra::Vector></tt> that isn't const-qualified
  /// exactly as declared here.)
  Diagonal (const Teuchos::RCP<const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& diag);

  //! Destructor
  virtual ~Diagonal();

  //! Sets parameters on this object.
  /**
    Currently doesn't need any parameters.
  */
  void setParameters (const Teuchos::ParameterList& params);

  //! Initialize
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(isInitialized_);
  }

  //! compute the preconditioner
  void compute();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return(isComputed_);
  }

  //! @name Methods implementing a Tpetra::Operator interface.
  //@{

  /// \brief Apply the preconditioner to X, putting the result in Y.
  ///
  /// If the result of applying this preconditioner to a vector X is
  /// \f$F \cdot X$, then this method computes \f$\beta Y + \alpha F \cdot X\f$.
  /// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! The Tpetra::Map representing this operator's domain.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getDomainMap () const {
    return domainMap_;
  }

  //! The Tpetra::Map representing this operator's range.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getRangeMap () const {
    return rangeMap_;
  }

  //@}
  //! \name Mathematical functions.
  //@{

  /// \brief Compute the condition number estimate and return its value.
  ///
  /// \warning This method is DEPRECATED.  It was inherited from
  ///   Ifpack, and Ifpack never clearly stated what this method
  ///   computes.  Furthermore, Ifpack's method just estimates the
  ///   condition number of the matrix A, and ignores the
  ///   preconditioner -- which is probably not what users thought it
  ///   did.  If there is sufficient interest, we might reintroduce
  ///   this method with a different meaning and a better algorithm.
  magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  magnitude_type TEUCHOS_DEPRECATED getCondEst() const {
    return condEst_;
  }

  //! Return the communicator associated with this matrix operator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Return a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >
  getMatrix () const {
    return matrix_;
  }

  //! Return the number of flops in the computation phase.
  double getComputeFlops() const;

  //! Return the number of flops for the application of the preconditioner.
  double getApplyFlops() const;

  //! Return the number of calls to initialize().
  int getNumInitialize() const;

  //! Return the number of calls to compute().
  int getNumCompute() const;

  //! Return the number of calls to apply().
  int getNumApply() const;

  //! Return the time spent in initialize().
  double getInitializeTime() const;

  //! Return the time spent in compute().
  double getComputeTime() const;

  //! Return the time spent in apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vector_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  bool isInitialized_;
  bool isComputed_;
  Teuchos::RCP<const map_type> domainMap_;
  Teuchos::RCP<const map_type> rangeMap_;
  Teuchos::RCP<const row_matrix_type> matrix_;
  Teuchos::RCP<const vector_type> inversediag_;
  Teuchos::ArrayRCP<size_t> offsets_;

  mutable int numInitialize_;
  mutable int numCompute_;
  mutable int numApply_;

  double initializeTime_;
  double computeTime_;
  double applyTime_;

  magnitude_type condEst_;
};

/** Function to construct a Diagonal preconditioner with vector input.
* The input vector is assumed to contain the equivalent of the inverted
* diagonal of a matrix.
*
* Example usage:<br>
* typedef Tpetra::CrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> TCrsMatrix;<br>
* typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> TVector;<br>
* typedef Tpetra::Preconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type> TPrec;
*
* Teuchos::RCP<TVector> myvec = ...
*
* Teuchos::RCP<TPrec> prec = Ifpack2::createDiagonalPreconditioner<TCrsMatrix,TVector>(myvec);
*
*/
template<class MatrixType,class VectorType>
Teuchos::RCP<Ifpack2::Diagonal<MatrixType> >
createDiagonalPreconditioner(const Teuchos::RCP<const VectorType>& invdiag)
{
  return Teuchos::rcp(new Ifpack2::Diagonal<MatrixType>(invdiag));
}

}//namespace Ifpack2

#endif
