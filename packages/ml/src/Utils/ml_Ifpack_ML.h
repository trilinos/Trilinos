/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

/*!
 * \file ml_Ifpack_ML.h
 *
 * \class Ifpack_ML
 *
 * \brief Wrapper for Ifpack_Preconditioner
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last update do Doxygen: 08-Mar-05.
 *
 */

#ifndef ML_IFPACK_ML_H
#define ML_IFPACK_ML_H

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
  virtual int SetParameters(Teuchos::ParameterList& MLList)
  {
    MLList_ = MLList;
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
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
                         Epetra_RowMatrix* Matrix = 0)
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
    ML_RETURN(MLPrec_->ApplyInverse(X, Y));
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
  virtual ostream& Print(std::ostream& os) const
  {
    return(os);
  }

  //! Sets the use of transpose 9NOT SUPPORTED)
  int SetUseTranspose(bool)
  {
    ML_CHK_ERR(-1);
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

}; // namespace ML_Epetra

#endif
#endif
