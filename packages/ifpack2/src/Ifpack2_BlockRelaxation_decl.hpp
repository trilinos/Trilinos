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

#ifndef IFPACK2_BLOCKRELAXATION_DECL_HPP
#define IFPACK2_BLOCKRELAXATION_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_Partitioner.hpp"

#include <Tpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <string>
#include <iostream>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {
enum RelaxationType {
  JACOBI,
  GS,
  SGS
};

//! Ifpack2::BlockRelaxation: defines relaxation preconditioners for Tpetra::RowMatrix objects.

/*! 
  The Ifpack2::BlockRelaxation class enables the construction of relaxation
  preconditioners for Tpetra::RowMatrix. Ifpack2::Relaxation 
  is derived from Ifpack2::Preconditioner, which is itself derived
  from Tpetra::Operator.
  Therefore this object can be used as preconditioner everywhere an
  apply() method is required in the preconditioning step.
 
This class enables the construction of the following simple preconditioners:
- Jacobi;
- Gauss-Seidel;
- Symmetric Gauss-Seidel.

<P>We now briefly describe the main features of the above preconditioners.
Consider a linear system of type
\f[
A x = b,
\f]
where \f$A\f$ is a square, real matrix, and \f$x, b\f$ are two real
vectors. We begin with the decomposition
\f[
A = D - E - F
\f]
where \f$D\f$ is the diagonal of A, \f$-E\f$ is the strict lower part, and
\f$-F\f$ is the strict upper part. It is assumed that the diagonal entries
of \f$A\f$ are different from zero.

<P>Given an starting solution \f$x_0\f$, an iteration of the (damped) Jacobi
method can be written in matrix form as follows:
\f[
x_{k+1} = \omega D^{-1}(E + F) x_k + D_{-1}b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter.

Using Ifpack2::Jacobi, the user can apply the specified number of sweeps
(\f$k_{max}\f$), and the damping parameter. If only one sweep is used, then
the class simply applies the inverse of the diagonal of A to the input
vector.

<P>Given a starting solution \f$x_0\f$, an iteration of the (damped) GaussSeidel
method can be written in matrix form as follows:
\f[
(D - E) x_{k+1} = \omega F x_k + b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
the Gauss-Seidel preconditioner can be defined as
\f[
P_{GS}^{-1} = (D - E)^{-1}.
\f]
Clearly, the role of E and F can be interchanged. This is what the
"relaxation: backward mode" option is for.

<P>For a list of supported parameters, please refer to the BlockRelaxation::setParameters method.

    \author Michael Heroux, SNL 9214.

    \date Last modified on 22-Jan-05.
*/
template<class MatrixType,class ContainerType>
class BlockRelaxation : virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {

public:
  typedef typename MatrixType::scalar_type         Scalar;
  typedef typename MatrixType::local_ordinal_type  LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type           Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  // \name Constructors and Destructors
  //@{

  //! BlockRelaxation constructor with given Tpetra::RowMatrix input.
  explicit BlockRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix);

  //! BlockRelaxation Destructor.
  virtual ~BlockRelaxation();

  //@}

  //@{ \name Preconditioner computation methods

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

  //! @name Methods implementing a Tpetra::Operator interface.
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

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  magnitudeType getCondEst() const;

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const;

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

  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:

  // @{ Internal methods

  //! Copy constructor (should never be used)
  BlockRelaxation(const BlockRelaxation<MatrixType,ContainerType> & RHS);

  //! operator= (should never be used)
  BlockRelaxation<MatrixType,ContainerType>& operator=(const BlockRelaxation<MatrixType,ContainerType>& RHS);

  virtual void ApplyInverseJacobi(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
				  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual void DoJacobi(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual void ApplyInverseGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual void DoGaussSeidel(Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual void ApplyInverseSGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual void DoSGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xtmp,
		     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;
 
  void ExtractSubmatrices();

  //@}

  // @{ Internal data and parameters

  //! reference to the matrix to be preconditioned
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Importer for parallel GS and SGS
  Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > Importer_;

  //! Weights for overlapping overlapped Jacobi only.
  Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > W_;
  // Level of overlap among blocks (for overlapped Jacobi only).
  int OverlapLevel_;
  //! Contains the (block) diagonal elements of \c Matrix.
  mutable std::vector<Teuchos::RCP<ContainerType> > Containers_;
  //! Contains information about non-overlapping partitions.
  Teuchos::RCP<Ifpack2::Partitioner<Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> > > Partitioner_;
  std::string PartitionerType_;

  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;

  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;
  //! Number of local blocks
  size_t NumLocalBlocks_;
  //! Which type of point relaxation approach to use
  int PrecType_;
  //! Minimum diagonal value
  Scalar MinDiagonalValue_;
  //! Damping factor.
  Scalar DampingFactor_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! Backward-Mode Gauss Seidel 
  bool DoBackwardGS_;
  //! Condition number estimate
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

