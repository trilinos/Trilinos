#ifndef IFPACK_POINTRELAXATION_H
#define IFPACK_POINTRELAXATION_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Preconditioner.h"
namespace Teuchos {
  class ParameterList;
}
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Time;
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_Import;

//! Ifpack_PointRelaxation: a class to define point relaxation preconditioners of Ifpack_Preconditioner's.

/*! 
  The Ifpack_PointRelaxation class enables the construction of point
  relaxation
  preconditioners of an Epetra_RowMatrix. Ifpack_PointRelaxation 
  is derived from 
  the Ifpack_Preconditioner class, which is derived from Epetra_Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.
 
This class enables the construction of the following simple preconditioners:
- Jacobi;
- Gauss-Seidel;
- symmetric Gauss-Seidel;
- SOR;
- SSOR.

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

<P>Given an starting solution \f$x_0\f$, an iteration of the Jacobi
method can be written in matrix form as follows:
\f[
x_{k+1} = \omega D^{-1}(E + F) x_k + D_{-1}b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter.

Using Ifpack_Jacobi, the user can apply the specified number of sweeps
(\f$k_{max}\f$), and the damping parameter. If only one sweep is used, then
the class simply applies the inverse of the diagonal of A to the input
vector.

<P>Given an starting solution \f$x_0\f$, an iteration of the GaussSeidel
method can be written in matrix form as follows:
\f[
(D - E) x_{k+1} = \omega F x_k + b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
the Gauss-Seidel preconditioner can be defined as
\f[
P_{GS}^{-1} = (D - E)^{-1}.
\f]
Clearly, the role of E and F can be interchanged. However,
Ifpack_GaussSeidel does not consider backward Gauss-Seidel methods.

<P>For a list of supported parameters, please refer to page \ref ifp_params.

\note The ApplyInverse() implementation of this class is \e not AztecOO
complaint, as it does assume that the two input vectors X and Y actually
refer to two different memory location. In fact, this case is handled
by class Ifpack_AdditiveSchwarz, which takes care of calling methods ApplyInverse()
of Ifpack_PointRelaxation with two vectors pointing to different memory
locations.

\author Marzio Sala, SNL 9214.

\date Last modified: Nov-04.
  
*/
class Ifpack_PointRelaxation : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_PointRelaxation constructor with given Epetra_RowMatrix.
  /*! Creates an instance of Ifpack_PointRelaxation class.
   *
   * \param In
   * Matrix - Pointer to matrix to precondition.
   */
  Ifpack_PointRelaxation(const Epetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Ifpack_PointRelaxation();

  //@}

  /*! This flag can be used to apply the preconditioner to the transpose of
   * the input operator. 
   * 
   * \return Integer error code, set to 0 if successful.  
   * Set to -1 if this implementation does not support transpose.
    */
  virtual inline int SetUseTranspose(bool UseTranspose)
  {
    UseTranspose_ = UseTranspose;
    return(0);
  }

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Epetra_MultiVector.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual inline int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO complaint.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

  virtual const char * Label() const
  {
    return(Label_.c_str());
  }

  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(UseTranspose_);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;

  virtual int Initialize();
  
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  virtual int Compute();

  //@}
 
  //@{ \name Miscellaneous

  virtual const Epetra_RowMatrix& Matrix() const 
  {
    return(*Matrix_);
  }

  //! Returns the condition number estimate.
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Epetra_RowMatrix* Matrix = 0);

  virtual double Condest() const
  {
    return(Condest_);
  }

  //! Sets all the parameters for the preconditioner
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Print object to an output stream
  //! Print method
  virtual ostream& Print(ostream & os) const;

  //@}

  //@{ \name Timing and flop count

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the computation phase.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the number of flops for the application of the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  // @}

private:
 
  // @{ Application of the preconditioner
  
  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseJacobi(const Epetra_MultiVector& X, 
                                 Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseGS(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSGS(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

#ifdef FIXME
  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSOR(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSSOR(const Epetra_MultiVector& X, 
                               Epetra_MultiVector& Y) const;
#endif
  //@}

private:
  
  //! Sets the label.
  virtual void SetLabel();

  //! Copy constructor (PRIVATE, should not be used)
  Ifpack_PointRelaxation(const Ifpack_PointRelaxation& rhs)
  {}
  
  //! operator = (PRIVATE, should not be used)
  Ifpack_PointRelaxation& operator=(const Ifpack_PointRelaxation& rhs)
  {
    return(*this);
  }

  // @{ Initializations, timing and flops
  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  // @}

  // @{ Settings
  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;
  //! Damping factor.
  double DampingFactor_;
  //! If true, use the tranpose of \c Matrix_.
  bool UseTranspose_;
  //! Contains the estimated condition number
  double Condest_;
  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
  //! Contains the label of this object.
  string Label_;
  int PrecType_;
  double MinDiagonalValue_;
  // @}

  // @{ Other data
  //! Number of local rows.
  int NumMyRows_;
  //! Number of local nonzeros.
  int NumMyNonzeros_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global nonzeros.
  int NumGlobalNonzeros_;
  //! Pointers to the matrix to be preconditioned.
  const Epetra_RowMatrix* Matrix_;
  //! Importer for parallel GS and SGS
  Epetra_Import* Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Epetra_Vector* Diagonal_;
  //! Time object to track timing.
  Epetra_Time* Time_;
  bool IsParallel_;
  bool ZeroStartingSolution_;
  bool UseWithAztecOO_;
  // @}

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_POINTRELAXATION_H
