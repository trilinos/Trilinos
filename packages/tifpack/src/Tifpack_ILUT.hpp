/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef TIFPACK_ILUT_HPP
#define TIFPACK_ILUT_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_CondestType.hpp"
#include "Tifpack_ScalingType.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Tifpack {

//! Tifpack::ILUT: A class for constructing and using an ILUT factorization
// of a given Tpetra_RowMatrix.

/*! The Tifpack::ILUT class computes a "Relaxed" ILUT factorization with level k fill 
    of a given Tpetra::RowMatrix. 

    <P> Please refer to \ref ifp_ilu for a general description of the ILU algorithm.

    <P>The complete list of supported parameters is reported in page \ref ifp_params. 

    \author Michael Heroux, SNL 9214.

    \date Last modified on 22-Jan-05.
*/    
class ILUT: public Tifpack::Preconditioner {
      
public:
  // @{ Constructors and Destructors
  //! ILUT constuctor with variable number of indices per row.
  ILUT(const Tpetra_RowMatrix* A);
  
  //! ILUT Destructor
  virtual ~ILUT();

  // @}
  // @{ Construction methods
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes five parameter names: level_fill, drop_tolerance,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive. For level_fill the ParameterEntry must have type int, the 
     threshold entries must have type double and overlap_mode must have type
     Tpetra_CombineMode.
  */
  int SetParameters(Teuchos::ParameterList& parameterlis);

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
    \param In 
           A - User matrix to be factored.
    \warning The graph of A must be identical to the graph passed in to Tifpack_IlukGraph constructor.
             
   */
  int Initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Compute IC factor U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the RILU(k) factors L and U using the current:
    <ol>
    <li> Tifpack_IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const {return(IsComputed_);};

  // Mathematical functions.
  
  //! Returns the result of a Tifpack_ILUT forward/back solve on a Tpetra_MultiVector X in Y.
  /*! 
    \param 
    X - (In) A Tpetra_MultiVector of dimension NumVectors to solve for.
    \param 
    Y - (Out) A Tpetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

  int Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

  //! Computed the estimated condition number and returns the value.
  double Condest(const Tifpack_CondestType CT = Tifpack_Cheap, 
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
		 Tpetra_RowMatrix* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  double Condest() const
  {
    return(Condest_);
  }

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
      does not support transpose use, this method should return a value of -1.
      
     \param
     UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

     \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Tpetra_Map object associated with the domain of this operator.
  const Tpetra_Map & OperatorDomainMap() const {return(A_.OperatorDomainMap());};

  //! Returns the Tpetra_Map object associated with the range of this operator.
  const Tpetra_Map & OperatorRangeMap() const{return(A_.OperatorRangeMap());};

  //! Returns the Tpetra_BlockMap object associated with the range of this matrix operator.
  const Tpetra_Comm & Comm() const{return(Comm_);};

  //! Returns a reference to the matrix to be preconditioned.
  const Tpetra_RowMatrix& Matrix() const
  {
    return(A_);
  }

  //! Returns a reference to the L factor.
  const Tpetra_CrsMatrix & L() const {return(*L_);};
  
  //! Returns a reference to the U factor.
  const Tpetra_CrsMatrix & U() const {return(*U_);};
    
  //! Returns the label of \c this object.
  const char* Label() const
  {
    return(Label_.c_str());
  }

  //! Sets the label for \c this object
  int SetLabel(const char* Label_in)
  {
    Label_ = Label_in;
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

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  inline double LevelOfFill() const {
    return(LevelOfFill_);
  }

  //! Set relative threshold value
  inline double RelaxValue() const {
    return(Relax_);
  }

  //! Get absolute threshold value
  inline double AbsoluteThreshold() const
  {
    return(Athresh_);
  }

  //! Get relative threshold value
  inline double RelativeThreshold() const
  {
    return(Rthresh_);
  }
    
  //! Gets the dropping tolerance
  inline double DropTolerance() const
  {
    return(DropTolerance_);
  }
    
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {
    // FIXME: diagonal of L_ should not be stored
    return(L().NumGlobalNonzeros() + U().NumGlobalNonzeros() - L().NumGlobalRows());
  }
 
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {
    return(L().NumMyNonzeros() + U().NumMyNonzeros());
  }

private:
  
  // @}
  // @{ Internal methods

  //! Copy constructor (should never be used)
  Tifpack_ILUT(const Tifpack_ILUT& RHS) :
    A_(RHS.Matrix()),
    Comm_(RHS.Comm()),
    Time_(Comm())
  {};

  //! operator= (should never be used)
  Tifpack_ILUT& operator=(const Tifpack_ILUT& RHS)
  {
    return(*this);
  }

  //! Releases all allocated memory.
  void Destroy();

  // @}
  // @{ Internal data

  //! reference to the matrix to be preconditioned.
  const Tpetra_RowMatrix& A_;
  //! Reference to the communicator object.
  const Tpetra_Comm& Comm_;
  //! L factor
  Teuchos::RefCountPtr<Tpetra_CrsMatrix> L_;
  //! U factor
  Teuchos::RefCountPtr<Tpetra_CrsMatrix> U_;
  //! Condition number estimate.
  double Condest_;
  //! relaxation value
  double Relax_;
  //! Absolute threshold
  double Athresh_;
  //! Relative threshold
  double Rthresh_;
  //! Level-of-fill
  double LevelOfFill_;
  //! Discards all elements below this tolerance
  double DropTolerance_;
  //! Label for \c this object
  string Label_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;
  //! \c true if transpose has to be used.
  bool UseTranspose_;
  //! Number of local rows.
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
  mutable double ApplyInverseFlops_;
  //! Used for timing purposes
  mutable Teuchos::Time Time_;
  //! Global number of nonzeros in L and U factors
  int GlobalNonzeros_;
  Teuchos::RefCountPtr<Tpetra_SerialComm> SerialComm_;
  Teuchos::RefCountPtr<Tpetra_Map> SerialMap_;
}; // Tifpack_ILUT

}//namespace Tifpack

#endif /* TIFPACK_ILUT_HPP */
