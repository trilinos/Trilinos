#ifndef IFPACK_POINTRELAXATION_H
#define IFPACK_POINTRELAXATION_H

#include "Ifpack_ConfigDefs.h"
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

//! Ifpack_PointRelaxation: a class to define point preconditioners of Ifpack_Preconditioner's.

/*! 
  The Ifpack_PointRelaxation class enables the construction of point
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

\note
Ifpack_PointRelaxation is \e not supposed to be used as stand-alone
preconditioner, but rather as a template for the construction
of Ifpack_AdditiveSchwarz objects. 

<P>The following parameters are recognized by Ifpack_PointRelaxation:
- \c "point: sweeps" [int, default = 1]: number of sweeps;
- \c "point: damping factor" [double, default = 1.0]: dampig parameter;
- \c "point: print frequency" [int, default = 0]: if different from zero, every specified
  number of sweeps the class computed the true residual, and prints on
  cout. As the computation of the residual is an expensive operation,
  this parameter should be considered for debugging purposed only;
- \c "point: min diagonal value" [double, default = 1e-9]: replace
  diagonal values below this value with this value.
\author Marzio Sala, SNL 9214.

\date Last modified: Oct-04.
  
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

  //! Copy constructor
  Ifpack_PointRelaxation(const Ifpack_PointRelaxation& rhs);
  
  //! Destructor.
  virtual ~Ifpack_PointRelaxation();

  // operator =
  Ifpack_PointRelaxation& operator=(const Ifpack_PointRelaxation& rhs);

  //@}

  //@{ \name Atribute set methods.

  //! If set true, applies the preconditioner to the transpose of the input operator.
  /*! This flag can be used to apply the preconditioner to the transpose of
   * the input operator. 
   * 
   * \return Integer error code, set to 0 if successful.  
   * Set to -1 if this implementation does not support transpose.
    */
  virtual int SetUseTranspose(bool UseTranspose)
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
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

  virtual char * Label() const
  {
    return((char*)Label_);
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
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  virtual int Compute();

  //! Returns a pointer to the matrix.
  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //@}

  //@{ \name Miscellaneous

  //! Returns the condition number estimate.
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Epetra_RowMatrix* Matrix = 0);

  virtual double Condest() const
  {
    return(Condest_);
  }

#ifdef HAVE_IFPACK_TEUCHOS
  //! Sets all the parameters for the preconditioner
  virtual int SetParameters(Teuchos::ParameterList& List);
#endif

  //! Sets the label.
  virtual void SetLabel();

  //! Print object to an output stream
  //! Print method
  virtual ostream& Print(ostream & os) const;

  //@}

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

  virtual long int ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual long int ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  virtual const Epetra_Vector* Diagonal() const
  {
    return(Diagonal_);
  }
  
  virtual bool ZeroStartingSolution() const
  {
    return(ZeroStartingSolution_);
  }

  virtual double MinDiagonalValue() const
  {
    return(MinDiagonalValue_);
  }

  virtual int PrecType() const
  {
    return(PrecType_);
  }

protected:
 
  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseJacobi(const Epetra_MultiVector& X, 
                                 Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseGS(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSGS(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  //! Utility function for SGS
  virtual int ApplyInverseSGS2(Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSOR(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSSOR(const Epetra_MultiVector& X, 
                               Epetra_MultiVector& Y) const;

  //! Utility function for SSOR
  virtual int ApplyInverseSSOR2(Epetra_MultiVector& Y) const;

  //@{ \name Setting functions

  //! Sets the number of sweeps.
  inline int SetNumSweeps(const int NumSweeps)
  {
    NumSweeps_ = NumSweeps;
    return(0);
  }

  //! Gets the number of sweeps.
  inline int NumSweeps() const
  {
    return(NumSweeps_);
  }
 
  //! Sets the damping factor
  inline int SetDampingFactor(const double DampingFactor)
  {
    DampingFactor_ = DampingFactor;
    return(0);
  }

  //! Gets the damping factor.
  inline double DampingFactor() const
  {
    return(DampingFactor_);
  }

  inline int SetPrintFrequency(const int PrintFrequency)
  {
    PrintFrequency_ = PrintFrequency;
    return(0);
  }

  inline int PrintFrequency() const
  {
    return(PrintFrequency_);
  }

  //! Sets the label.


  //@}

  //! Contains the diagonal elements of \c Matrix.
  mutable Epetra_Vector* Diagonal_;
  //! Contains the label of this object.
  char Label_[80];
  //! If true, use zero vector as starting solution.
  bool ZeroStartingSolution_;
  
private:

  //! Extracts a copy of the diagonal, stores the elements in \c Diagonal_.
  int ExtractDiagonal();

  //! Returns the i-th element stored in Diagonal_.
  double Diagonal(const int i);
  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;
  //! Damping factor.
  double DampingFactor_;
  //! Number of local rows.
  int NumMyRows_;
  //! Pointers to the matrix to be preconditioned.
  const Epetra_RowMatrix* Matrix_;
  //! If true, use the tranpose of \c Matrix_.
  bool UseTranspose_;
  //! Toggles the frequency, 0 means no output.
  int PrintFrequency_;
  //! Contains the estimated condition number
  double Condest_;
  //! If true, Compute() also computed the condition number estimate.
  bool ComputeCondest_;

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
  Epetra_Time* Time_;

  //! Contains the number of flops for Compute().
  long int ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  long int ApplyInverseFlops_;
  int PrecType_;
  double MinDiagonalValue_;
};

#endif // IFPACK_POINTRELAXATION_H
