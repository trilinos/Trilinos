/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CHEBYSHEV_DECL_HPP
#define IFPACK2_CHEBYSHEV_DECL_HPP


#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <iostream>
#include <string>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {

/// \class Chebyshev 
/// \brief Diagonally scaled Chebyshev iteration for Tpetra sparse matrices.
/// \tparam MatrixType A specialization of Tpetra::CrsMatrix.
///
/// \section Ifpack_Chebyshev_Summary Summary
///
/// This class implements a Chebyshev polynomial preconditioner or
/// smoother for a Tpetra::RowMatrix or Tpetra::CrsMatrix.  Chebyshev
/// is derived from Preconditioner, which itself is derived from
/// Tpetra::Operator.  Therefore, this object can be used as
/// preconditioner everywhere an apply() method is required in the
/// preconditioning step.
///
/// \warning Our implementation currently <i>only</i> works with a
///   real symmetric positive definite (SPD) matrix.  Results for
///   matrices that are not SPD, or for complex-valued Scalar types,
///   are not defined.
///
/// \section Ifpack_Chebyshev_Algorithm Algorithm
///
/// This algorithm requires that the matrix be symmetric positive
/// definite.  As a result, all of its eigenvalues must lie in a
/// positive interval on the real line, \f$[\lambda_{min},
/// \lambda_{max}]\f$.  We require that users give us at least (an
/// estimate of) the maximum eigenvalue \f$\lambda_{max}\f$.  They may
/// optionally also give us the (estimated) ratio \f$\eta =
/// \lambda_{max} / \lambda_{min}\, or (an estimate of) the minimum
/// eigenvalue \f$\lambda_{min}\f$.
///
/// When using Chebyshev iteration to solve linear systems directly,
/// it is important to have good estimates of both the minimum and
/// maximum eigenvalues.  However, when using a small number of
/// Chebyshev iterations as a smoother in multigrid, the maximum
/// eigenvalue estimate is more important.  (The point of a smoother
/// is to smooth out the high-frequency components of the error.)
/// This is why we use a ratio \f$\eta = \lambda_{max} /
/// \lambda_{min}$, rather than requiring a guess for $\lambda_{min}$.
/// In fact, we only use \f$\lambda_{min}\f$ for error checking, not
/// when determining the Chebyshev coefficients.  Often, if users give
/// us \f$\lambda_{max}\f$, our default value of \f$\eta\f$ suffices.
///
/// Some Chebyshev implementations attempt to estimate the eigenvalue
/// interval automatically.  Steve Ashby's CHEBYCODE is the original
/// example.  We do not attempt to do this.  Users who want estimates
/// of the smallest and largest eigenvalues should run a few
/// iterations of Lanczos or CG.  Since the largest eigenvalue is more
/// important for smoother applications, a few iterations of the power
/// method may be enough.
///
/// Call the setParameters() method to give this instance your
/// estimates of \f$\lambda_{max}\f$ and \f$\eta\f$, as well as to set
/// other options controlling the behavior of Chebyshev iteration.
///
/// \section Ifpack_Chebyshev_Performance Performance
///
/// Chebyshev should spend most of its time in Tpetra's native sparse
/// matrix-vector multiply kernel.  This should give good performance,
/// since we have spent a lot of effort tuning that kernel.  Depending
/// on the Kokkos Node type of your Tpetra matrix, the kernel may also
/// exploit threads for additional parallelism within each MPI process
/// ("hybrid parallelism" a.k.a. "MPI + X").  If your application
/// depends on hybrid parallelism for performance, you should favor
/// Chebyshev smoothers whenever possible over "serial within a
/// process" smoothers like Gauss-Seidel or SOR (Symmetric
/// Over-Relaxation).
///
/// \section Ifpack_Chebyshev_History History
///
/// The original implementation of this class was an adaptation of
/// ML's ML_Cheby routine.  The original author was Ulrich Hetmaniuk,
/// a Sandia employee in what was then (2006) Org 1416.
template<class MatrixType>
class Chebyshev : 
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
					   typename MatrixType::local_ordinal_type,
					   typename MatrixType::global_ordinal_type,
					   typename MatrixType::node_type> 
{
public:
  //! \name Typedefs
  //@{ 
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //@}
  // \name Constructors and destructors
  //@{

  //! Chebyshev constructor with given Tpetra::RowMatrix input.
  explicit Chebyshev (const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);

  //! Destructor.
  virtual ~Chebyshev ();

  //@}
  //! \name Preconditioner computation methods
  //@{ 

  /// \brief Set parameters for the preconditioner.
  ///
  /// The following parameters control the Chebyshev coefficients, the
  /// number of iterations, various other properties of the algorithm,
  /// and error checking.
  /// - "chebyshev: max eigenvalue" (\c Scalar): (An estimate of) the
  ///   largest eigenvalue \f$\lambda_{max}\f$ of the matrix.  You
  ///   should always provide this value, since otherwise Chebyshev
  ///   will not be an effective smoother.
  /// - "chebyshev: ratio eigenvalue" (\c Scalar): (Estimated) ratio
  ///   \f$\eta\f$ between the largest and smallest eigenvalue.  The
  ///   default value is 30.0.  We use this to define the interval on
  ///   the real line that determines the Chebyshev coefficients.
  /// - "chebyshev: min eigenvalue" (\c Scalar): (An estimate of) the
  ///   smallest eigenvalue \f$\lambda_{min}\f$ of the matrix.  This
  ///   parameter is optional and is only used to check whether the
  ///   input matrix is the identity matrix.  We do <i>not</i> use
  ///   this to define the interval on the real line that determines
  ///   the Chebyshev coefficients.
  /// - "chebyshev: degree" (\c int): The polynomial degree; the
  ///   number of times that apply() invokes sparse matrix-vector
  ///   multiply.
  /// - "chebyshev: min diagonal value" (\c Scalar): Threshold for
  ///   diagonal values.  Values smaller than this are not inverted.
  /// - "chebyshev: zero starting solution" (\c bool): If true, the
  ///   input vector(s) is/are filled with zeros on entry to the
  ///   apply() method.
  /// - "chebyshev: operator inv diagonal" (<tt>Tpetra::Vector</tt>).
  ///   A (raw) pointer to the inverse of the diagonal entries of the
  ///   matrix.  If not provided, we compute this ourselves.  If
  ///   provided, we make a deep copy.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Initialize the preconditioner.
  ///
  /// You must call this method before you may call compute() or apply().
  void initialize();

  //! Whether the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  /// \brief Compute the preconditioner.
  ///
  /// You must call this method after calling initialize(), but before
  /// calling apply().  If any of the diagonal entries of the matrix
  /// have changed, you must call this method again before calling
  /// apply().
  void compute();

  //! Whether the preconditioner has been successfully computed.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}
  //! @name Implementation of Tpetra::Operator
  //@{ 

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// \param X [in] A (multi)vector to which to apply the preconditioner.
  /// \param Y [in/out] A (multi)vector containing the result of
  ///   applying the preconditioner to X.
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  bool hasTransposeApply() const;

  //! Applies the matrix to a Tpetra::MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void applyMat(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                      Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}

  //@{
  //! \name Mathematical functions.

  //! Computes the estimated condition number and returns the value.
  magnitudeType computeCondEst(CondestType CT = Cheap, 
                               LocalOrdinal MaxIters = 1550,
                               magnitudeType Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

  //@}

  //@{ 
  //! \name Attribute accessor methods

  //! The estimated condition number, or -1.0 if it has not yet been computed.
  magnitudeType getCondEst() const;

  //! The communicator over which the matrix is distributed.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const;

  //! The total number of flops taken by all calls to compute().
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

  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

  // @{ \name Utility methods

  //! Simple power method to compute lambda_max.
  static void PowerMethod(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator,
                         const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal,
                         const int MaximumIterations, 
                         Scalar& LambdaMax);

  //! Not currently implemented: Use CG to estimate lambda_min and lambda_max.
  static void CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
                const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
                const int MaximumIterations, 
                Scalar& lambda_min, Scalar& lambda_max);

  //@}

private:
  
  //! Copy constructor (should never be used)
  Chebyshev(const Chebyshev<MatrixType>& src);

  //! operator= (should never be used)
  Chebyshev<MatrixType>& operator=(const Chebyshev<MatrixType>& src);

  // @{ Internal data and parameters

  //! reference to the matrix to be preconditioned
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Reference to the communicator object
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! Contains the inverse of diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InvDiagonal_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Contains the degree of Chebyshev polynomial.
  int PolyDegree_;
  //! The ratio such that [LambdaMax_ / EigRatio_, LambdaMax_], the interval of interest for the Chebyshev polynomial.
  magnitudeType EigRatio_;
  //! An approximation to the smallest eigenvalue.
  Scalar LambdaMin_;
  //! An approximation to the largest eigenvalue.
  Scalar LambdaMax_;
  //! The minimum value on the diagonal.
  Scalar MinDiagonalValue_;
  //! The estimated condition number
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  magnitudeType Condest_;
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
  //! Contains the number of flops for apply().
  mutable double ApplyFlops_;
  //! Number of local rows.
  size_t NumMyRows_;
  //! Number of global rows.
  global_size_t NumGlobalRows_;
  //! Number of global nonzeros.
  global_size_t NumGlobalNonzeros_;

  //@}

}; // class Chebyshev

}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DECL_HPP

