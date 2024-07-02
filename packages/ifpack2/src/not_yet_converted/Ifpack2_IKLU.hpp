// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_IKLU_HPP
#define IFPACK2_IKLU_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_CondestType.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_IKLU_Utils.hpp"	
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Time.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_RowMatrix;
class Tpetra_SerialComm;
class Tpetra_Comm;
class Tpetra_Map;
class Tpetra_MultiVector;

namespace Teuchos {
  class ParameterList;
}

//! Ifpack2_IKLU: A class for constructing and using an incomplete Cholesky factorization of a given Tpetra_RowMatrix.

/*! The Ifpack2_IKLU class computes a "Relaxed" IKLU factorization with level k fill 
    of a given Tpetra_RowMatrix. 

    <P> Please refer to \ref ifp_ilu for a general description of the ILU algorithm.

    <P>The complete list of supported parameters is reported in page \ref ifp_params. 

    \author Heidi Thornquist, Org. 1437

    \date Last modified on 28-Nov-06.
*/    
class Ifpack2_IKLU: public Ifpack2_Preconditioner {
      
public:
  // @{ Constructors and Destructors
  //! Ifpack2_IKLU constuctor with variable number of indices per row.
  Ifpack2_IKLU(const Tpetra_RowMatrix* A);
  
  //! Ifpack2_IKLU Destructor 
  virtual ~Ifpack2_IKLU();

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
    \warning The graph of A must be identical to the graph passed in to Ifpack2_IlukGraph constructor.
             
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
    <li> Ifpack2_IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const {return(IsComputed_);};

  // Mathematical functions.
  
  //! Returns the result of a Ifpack2_IKLU forward/back solve on a Tpetra_MultiVector X in Y.
  /*! 
    \param 
    X - (In) A Tpetra_MultiVector of dimension NumVectors to solve for.
    \param 
    Y - (Out) A Tpetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Computed the estimated condition number and returns the value.
  double Condest(const Ifpack2_CondestType CT = Ifpack2_Cheap, 
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
  Ifpack2_IKLU(const Ifpack2_IKLU& RHS) :
    A_(RHS.Matrix()),
    Comm_(RHS.Comm()),
    Time_(Comm())
  {};

  //! operator= (should never be used)
  Ifpack2_IKLU& operator=(const Ifpack2_IKLU& RHS)
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
  Teuchos::RCP<Tpetra_CrsMatrix> L_;
  //! U factor
  Teuchos::RCP<Tpetra_CrsMatrix> U_;
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
  //! Number of local nonzeros
  int NumMyNonzeros_;
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
  //! Used for timing purposed
  mutable Tpetra_Time Time_;
  //! Global number of nonzeros in L and U factors
  int GlobalNonzeros_;
  Teuchos::RCP<Tpetra_SerialComm> SerialComm_;
  Teuchos::RCP<Tpetra_Map> SerialMap_;

  //! Containers for the matrix storage and permutation
  csr* csrA_;
  css* cssS_;
  //! Container for the L and U factor
  csrn* csrnN_;

}; // Ifpack2_IKLU

#endif /* IFPACK2_IKLU_HPP */
