/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef IFPACK_AMESOS_H
#define IFPACK_AMESOS_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Map;
class Epetra_Time;
class Epetra_Comm;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_RowMatrix;

//! Ifpack_Amesos: a class to use Amesos' factorizations as preconditioners.
/*!
Class Ifpack_Amesos enables the use of Amesos' factorizations as
Ifpack_Preconditioners.

Ifpack_Amesos is just a bare-bone wrap to Amesos. Currently, the
only parameter required recognized by SetParameters() is
\c "amesos: solver type" (defaulted to \c
"Amesos_Klu"), which defined the Amesos solver. The Teuchos list
in input to SetParameters() is copied, then the copied list is
used to set the parameters of the Amesos object.

This class works with matrices whose communicator contains only one
process, that is, either serial matrices, or Ifpack_LocalFilter'd matrices.

\warning The number of flops is NOT updated.

\author Marzio Sala, SNL 9214.

\date Last update Oct-04.

*/
class Ifpack_Amesos : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors.

  //! Constructor.
  Ifpack_Amesos(Epetra_RowMatrix* Matrix);

  //! Copy constructor.
  Ifpack_Amesos(const Ifpack_Amesos& rhs);

  //! Operator=.
  Ifpack_Amesos& operator=(const Ifpack_Amesos& rhs);

  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_Amesos() {};

  //@}

  //@{ \name Attribute set methods.

   //! If set true, transpose of this operator will be applied (not implemented).
    /*! This flag allows the transpose of the given operator to be used
     * implicitly.

    \param UseTranspose_in - (In) If true, multiply by the transpose of operator,
           otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */

  virtual int SetUseTranspose(bool UseTranspose_in);
  //@}

  //@{ \name Mathematical functions.

    //! Applies the matrix to an Epetra_MultiVector.
  /*!
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Applies the preconditioner to X, returns the result in Y.
  /*!
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method
    must support the case where X and Y are the same object.
    */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix (not implemented)
    virtual double NormInf() const;
  //@}

  //@{ \name Attribute access functions

    //! Returns a character string describing the operator
    virtual const char * Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm & Comm() const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map & OperatorDomainMap() const;

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map & OperatorRangeMap() const;

  //@}

  //@{ \name Construction and application methods.

  //! Returns \c true is the preconditioner has been successfully initialized.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Initializes the preconditioners.
  /*! \return
   * 0 if successful, 1 if problems occurred.
   */
  virtual int Initialize();

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  /*! \return
   * 0 if successful, 1 if problems occurred.
   */
  virtual int Compute();

  //! Sets all the parameters for the preconditioner.
  /*! Parameters currently supported:
   * - \c "amesos: solver type" : Specifies the solver type
   *   for Amesos. Default: \c Amesos_Klu.
   *
   * The input list will be copied, then passed to the Amesos
   * object through Amesos::SetParameters().
   */
  virtual int SetParameters(Teuchos::ParameterList& List);

  //@}

  //@{ \name Query methods.

  //! Returns a const reference to the internally stored matrix.
  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Returns the estimated condition number, computes it if necessary.
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
                         Epetra_RowMatrix* Matrix_in= 0);

  //! Returns the estimated condition number, never computes it.
  virtual double Condest() const
  {
    return(Condest_);
  }

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

  //! Returns the total time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the total time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the total time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  //! Returns the total number of flops to computate the preconditioner.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the total number of flops to apply the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  // Returns a constant reference to the internally stored
  virtual const Teuchos::ParameterList& List() const
  {
    return(List_);
  }

  //! Prints on ostream basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os) const;

  //@}

protected:

  //@{ \name Methods to get/set private data

  //! Sets the label.
  inline void SetLabel(const char* Label_in)
  {
    Label_ = Label_in;
  }

  //! Sets \c IsInitialized_.
  inline void SetIsInitialized(const bool IsInitialized_in)
  {
    IsInitialized_ = IsInitialized_in;
  }

  //! Sets \c IsComputed_.
  inline void SetIsComputed(const int IsComputed_in)
  {
    IsComputed_ = IsComputed_in;
  }

  //! Sets \c NumInitialize_.
  inline void SetNumInitialize(const int NumInitialize_in)
  {
    NumInitialize_ = NumInitialize_in;
  }

  //! Sets \c NumCompute_.
  inline void SetNumCompute(const int NumCompute_in)
  {
    NumCompute_ = NumCompute_in;
  }

  //! Sets \c NumApplyInverse_.
  inline void SetNumApplyInverse(const int NumApplyInverse_in)
  {
    NumApplyInverse_ = NumApplyInverse_in;
  }

  //! Sets \c InitializeTime_.
  inline void SetInitializeTime(const double InitializeTime_in)
  {
    InitializeTime_ = InitializeTime_in;
  }

  //! Sets \c ComputeTime_.
  inline void SetComputeTime(const double ComputeTime_in)
  {
    ComputeTime_ = ComputeTime_in;
  }

  //! Sets \c ApplyInverseTime_.
  inline void SetApplyInverseTime(const double ApplyInverseTime_in)
  {
    ApplyInverseTime_ = ApplyInverseTime_in;
  }

  //! Sets \c ComputeFlops_.
  inline void SetComputeFlops(const double ComputeFlops_in)
  {
    ComputeFlops_ = ComputeFlops_in;
  }

  //! Sets \c ComputeFlops_.
  inline void SetApplyInverseFlops(const double ApplyInverseFlops_in)
  {
    ApplyInverseFlops_ = ApplyInverseFlops_in;
  }

  //! Set \c List_.
  inline void SetList(const Teuchos::ParameterList& List_in)
  {
    List_ = List_in;
  }
  //@}

private:

  //! Pointers to the matrix to be preconditioned.
  Teuchos::RefCountPtr<const Epetra_RowMatrix> Matrix_;

  //! Linear problem required by Solver_.
  Teuchos::RefCountPtr<Epetra_LinearProblem> Problem_;
  //! Amesos solver, use to apply the inverse of the local matrix.
  Teuchos::RefCountPtr<Amesos_BaseSolver> Solver_;
  //! Contains a copy of the input parameter list.
  Teuchos::ParameterList List_;

  //! Contains the label of \c this object.
  std::string Label_;
  //! If true, the linear system on this processor is empty, thus the preconditioner is null operation.
  bool IsEmpty_;
  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! If true, the preconditioner solves for the transpose of the matrix.
  bool UseTranspose_;

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
  //! Time object.
  Teuchos::RefCountPtr<Epetra_Time> Time_;

  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  double ApplyInverseFlops_;

  //! Contains the estimated condition number.
  double Condest_;
};

#endif // IFPACK_AMESOS_H
