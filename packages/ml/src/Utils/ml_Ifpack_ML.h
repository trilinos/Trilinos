/*!
 * \file ml_Ifpack_ML.h
 *
 * \class Ifpack_ML
 *
 * \brief Wrapper for Ifpack_Preconditioner
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last update to Doxygen: 08-Mar-05.
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_IFPACK_ML_H
#define ML_IFPACK_ML_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)

#include "ml_epetra.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Amesos.h"
#include "ml_MultiLevelPreconditioner.h"

namespace ML_Epetra {

/*!
 * \class Ifpack_ML
 *
 * \brief Wraps an ML preconditioner as an Ifpack_Preconditioner
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 14-Mar-05.
 *
 */

class Ifpack_ML : public Ifpack_Preconditioner {

public:

  //! Constructor.
  Ifpack_ML(Epetra_RowMatrix* A) :
    A_(A),
    MLPrec_(0)
  {};

  //! Destructor.
  virtual ~Ifpack_ML()
  {
    if (MLPrec_)
      delete MLPrec_;
  }

  //! Sets all the parameters for the preconditioner from the list.
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    std::string listName = List.get("ML sublist name","ML list");
    try{MLList_ = List.sublist(listName,true);}
    catch(...) {
      if (A_->Comm().MyPID()==0)
        std::cout << "Did not find sublist \"" << listName
             << "\" for ML subdomain solver.  Setting \"SA\" defaults." << std::endl;
      SetDefaults("SA",MLList_);
    };
    return(0);
  }

  //! Initialize the preconditioner.
  virtual int Initialize()
  {
    return(0);
  };

  //! Returns true if the  preconditioner has been successfully initialized, false otherwise
  virtual bool IsInitialized() const
  {
    return(true);
  }

  //! Computes all it is necessary to apply the preconditioner.
  virtual int Compute()
  {
    if (MLPrec_)
      delete MLPrec_;

    MLPrec_ = new ML_Epetra::MultiLevelPreconditioner(*A_, MLList_);
    if (MLPrec_->IsPreconditionerComputed() == false) {
      ML_CHK_ERR(-1);
    }
    else
      return(0);
  }

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool IsComputed() const
  {
    return(MLPrec_->IsPreconditionerComputed());
  }

  //! Computes the condition number estimate, returns its value.
  virtual double Condest(const Ifpack_CondestType /* CT */ = Ifpack_Cheap,
                         const int /* MaxIters */ = 1550,
                         const double /* Tol */ = 1e-9,
                         Epetra_RowMatrix* /* matrix */ = 0)
  {
    return(-1.0);
  }

  //! Returns the computed condition number estimate, or -1.0 if not computed.
  virtual double Condest() const
  {
    return(-1.0);
  }


  //! Applies the preconditioner to vector X, returns the result in Y.
  virtual int ApplyInverse(const Epetra_MultiVector& X,
                           Epetra_MultiVector& Y) const
  {
    int i = MLPrec_->ApplyInverse(X, Y);
    ML_RETURN(i);
  }

  //! Returns a pointer to the matrix to be preconditioned.
  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*A_);
  }

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(-1);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(-1);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(-1);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(0.0);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(0.0);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(0.0);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  //! Returns the number of flops in the computation phase.
  virtual double ComputeFlops() const
  {
    return(0.0);
  }

  //! Returns the number of flops in the application of the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(0.0);
  }

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& Print(std::ostream& os) const
  {
    return(os);
  }

  //! Sets the use of transpose (NOT SUPPORTED)
  int SetUseTranspose(bool useTranspose)
  {
    if (useTranspose) {
      ML_CHK_ERR(-1);
    }
    return 0;
  }

  //! Applies the matrix to a vector (NOT SUPPORTED)
  int Apply(const Epetra_MultiVector&, Epetra_MultiVector&) const
  {
    ML_CHK_ERR(-1);
  }

  //! Returns the norm inf (NOT SUPPORTED)
  double NormInf() const
  {
    return(-1.0);
  }

  //! Returns the label of \c this object.
  const char* Label() const
  {
    return("Ifpack_ML");
  }

  //! Returns \c true if the transpose is used.
  bool UseTranspose() const
  {
    ML_CHK_ERR(-1);
  }

  //! Returns \c true if the class furnishes an infinite norm.
  bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a reference to the communicator of \c this object.
  const Epetra_Comm& Comm() const
  {
    return(A_->Comm());
  }

  //! Returns a reference to the operator domain map.
  const Epetra_Map& OperatorDomainMap() const
  {
    return(A_->OperatorDomainMap());
  }

  //! Returns a reference to the operator range map.
  const Epetra_Map& OperatorRangeMap() const
  {
    return(A_->OperatorRangeMap());
  }

private:
  //! Pointer to the matrix used to build the preconditioner.
  Epetra_RowMatrix* A_;
  //! Pointer to the ML preconditioner.
  ML_Epetra::MultiLevelPreconditioner* MLPrec_;
  //! Copy of the input parameter list.
  Teuchos::ParameterList MLList_;
}; // class Ifpack_ML

} // namespace ML_Epetra

#endif
#endif
