#ifndef IFPACK_POINTPRECONDITIONER_H
#define IFPACK_POINTPRECONDITIONER_H

class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#endif

//! Ifpack_PointPreconditioner: a class to define point preconditioners of Ifpack_Preconditioner's.

/*! 
  The Ifpack_PointPreconditioner class enables the construction of point
  preconditioners of an Epetra_RowMatrix. Ifpack_PointPreconditioner 
  is derived from 
  the Ifpack_Preconditioner class, which is derived from Epetra_Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.
 
  This class is semi-virtual, as it does not implement ApplyInverse().
  Ifpack_PointPreconditioner is intended to furnish the basic utilities
  for all preconditioners like Jacobi, Gauss-Seidel, symmetric
  Gauss-Seidel, SOR, SSOR, etc.

 \date Sep-04
  
*/
class Ifpack_PointPreconditioner : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_PointPreconditioner constructor with given Epetra_RowMatrix.
  /*! Creates an instance of Ifpack_PointPreconditioner class.
   *
   * \param In
   * Matrix - Pointer to matrix to precondition.
   */
  Ifpack_PointPreconditioner(const Epetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Ifpack_PointPreconditioner();

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
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

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
  virtual double Condest() const
  {
    return(Condest_);
  }

#ifdef HAVE_IFPACK_TEUCHOS
  //! Sets all the parameters for the preconditioner
  virtual int SetParameters(Teuchos::ParameterList& List);
#endif

  //! Sets integer parameters
  virtual int SetParameter(const string Name, const int value)
  {
    if (Name == "point: sweeps")
      SetNumSweeps(value);
    else if (Name == "point: print frequency")
      SetPrintFrequency(value);
    else if (Name == "point: zero starting solution") {
      if (value)
	ZeroStartingSolution_ = true;
      else
	ZeroStartingSolution_ = false;
    }
    SetLabel();

    return(0);
  }

  //! Sets double parameters
  virtual int SetParameter(const string Name, const double value)
  {
    if (Name == "point: damping factor")
      SetDampingFactor(value);

    SetLabel();

    return(0);
  }

  //! Sets the label.
  virtual void SetLabel() = 0;

  //! Print object to an output stream
  //! Print method
  virtual ostream& Print(ostream & os) const;

  //@}

protected:
 
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
  //! If \c true, the Jacobi preconditioner has been computed successfully.
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

};

#endif // IFPACK_POINTPRECONDITIONER_H
