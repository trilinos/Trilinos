/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef IFPACK_ILUT_H
#define IFPACK_ILUT_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_CondestType.h"
#include "Ifpack_ScalingType.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
class Epetra_RowMatrix;
class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;

#ifdef HAVE_IFPACK_TEUCHOS
namespace Teuchos {
  class ParameterList;
}
#endif

//! Ifpack_ILUT: A class for constructing and using an incomplete Cholesky factorization of a given Epetra_RowMatrix.
// FIXME: role of Threshold ?? Athresh ? Rthresh??

class Ifpack_ILUT: public Ifpack_Preconditioner {
      
 public:
  //! Ifpack_ILUT constuctor with variable number of indices per row.
  Ifpack_ILUT(const Epetra_RowMatrix* A);
  
  Ifpack_ILUT(const Ifpack_ILUT& rhs);

  Ifpack_ILUT& operator=(const Ifpack_ILUT& rhs);

  //! Ifpack_ILUT Destructor
  virtual ~Ifpack_ILUT();

  //! Set absolute threshold value
  void SetLevelOfFill(double LevelOfFill) {
    LevelOfFill_ = LevelOfFill; return;
  }

  inline double LevelOfFill() const {
    return(LevelOfFill_);
  }

#ifdef FIXME
  //! Set relative threshold value
  int SetRelax(double Relax) {
    Relax_ = Relax; return;
    return(0);
  }

  //! Set relative threshold value
  double Relax() const {
    return(Relax_);
  }
#endif

#ifdef HAVE_IFPACK_TEUCHOS
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes five parameter names: level_fill, drop_tolerance,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive. For level_fill the ParameterEntry must have type int, the 
     threshold entries must have type double and overlap_mode must have type
     Epetra_CombineMode.
  */
  int SetParameters(Teuchos::ParameterList& parameterlis);
#endif

  const Epetra_RowMatrix& Matrix() const
  {
    return(A_);
  }

  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
    \param In 
           A - User matrix to be factored.
    \warning The graph of A must be identical to the graph passed in to Ifpack_IlukGraph constructor.
             
   */
  int Initialize();

  //! Compute IC factor U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the RILU(k) factors L and U using the current:
    <ol>
    <li> Ifpack_IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const {return(IsComputed_);};

  // Mathematical functions.
  
  //! Returns the result of a Ifpack_ILUT forward/back solve on a Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the maximum over all the condition number estimate for each local ILUT set of factors.
  /*! This functions computes a local condition number estimate on each processor and return the
      maximum over all processor of the estimate.
   \param In
    Trans -If true, solve transpose problem.
    \param Out
    ConditionNumberEstimate - The maximum across all processors of 
    the infinity-norm estimate of the condition number of the inverse of LDU.
  */
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap, 
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
		 Epetra_RowMatrix* Matrix = 0);

  double Condest() const
  {
    return(Condest_);
  }

  // Atribute access functions
  
  //! Get absolute threshold value
  double AbsoluteThreshold() const
  {
    return Athresh_;
  }

  //! Get relative threshold value
  double RelativeThreshold() const
  {
    return Rthresh_;
  }
    
  //! Set absolute threshold value
  int SetAbsoluteThreshold(const double Athresh)
  {
    Athresh_ = Athresh;
    return(0);
  }

  //! Set relative threshold value
  int SetRelativeThreshold(const double Rthresh)
  {
    Rthresh_ = Rthresh;
    return(0);
  }
    
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {
    return(L().NumGlobalNonzeros() + U().NumGlobalNonzeros());
  }
 
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {
    return(L().NumMyNonzeros() + U().NumMyNonzeros());
  }

  const Epetra_CrsMatrix & L() const {return(*L_);};
  const Epetra_CrsMatrix & U() const {return(*U_);};
    
  //@{ \name Additional methods required to support the Epetra_Operator interface.

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

    //! Returns 0.0 because this class cannot compute Inf-norm.
    double NormInf() const {return(0.0);};

    //! Returns false because this class cannot compute an Inf-norm.
    bool HasNormInf() const {return(false);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(UseTranspose_);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(A_.OperatorDomainMap());};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const{return(A_.OperatorRangeMap());};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    const Epetra_Comm & Comm() const{return(Comm_);};
  //@}

    char* Label() const
    {
      return((char*)Label_);
    }

    int SetLabel(const char* Label)
    {
      strcpy(Label_,Label);
      return(0);
    }
 
  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const;

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

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }


 private:
  
  const Epetra_RowMatrix& A_;
  const Epetra_Comm& Comm_;
  Epetra_CrsMatrix* L_;
  Epetra_CrsMatrix* U_;

  double Condest_;
  double Relax_;
  double Threshold_;
  double Athresh_;
  double Rthresh_;
  double LevelOfFill_;
  
  char Label_[160];

  bool IsInitialized_;
  bool IsComputed_;
  bool UseTranspose_;

  int NumMyRows_;

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
  double ApplyInverseFlops_;

mutable Epetra_Time Time_;
};

#endif /* IFPACK_ILUT_H */
