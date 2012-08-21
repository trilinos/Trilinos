#ifndef IFPACK_SUPPORTGRAPH_H
#define IFPACK_SUPPORTGRAPH_H

#include "Ifpack_Preconditioner.h"
#include "Ifpack_CondestType.h"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ScalingType.h"

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"

#include "Teuchos_RefCountPtr.hpp"

#include "Amesos.h"



class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
namespace Teuchos{
  class ParameterList;
}


class Ifpack_SupportGraph: public Ifpack_Preconditioner {

 public:

  //@{ \name Constructors/Destructors.
 
  //! Constructor
  Ifpack_SupportGraph(Epetra_RowMatrix* A);
 
  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_SupportGraph();

  //@}

  //@{ \name Attribute set methods.
  //! If set true, transpose of this operator will be applied (not implemented).              
  /*! This flag allows the transpose of the given operator to be used                         
   * implicitly.                                                                                 

   \param    
   UseTranspose_in - (In) If true, multiply by the transpose of operator, 
   otherwise just use operator.   
   
   \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does
  */

  int SetUseTranspose(bool UseTranspose_in);

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
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

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
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)                        
  virtual double NormInf() const {return(0.0);}
  
  //@} 

  //@{ \name Attribute access functions.
  
  //! Returns a character string describing the operator.
  virtual const char* Label() const;

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns true if this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const {return(false);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(Comm_);};

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(A_->OperatorDomainMap());};

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(A_->OperatorRangeMap());};

  //@}



  //@{ \name Construction and application methods.

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Initialize the preconditioner
  /*! \return
   * 0 if successful, 1 if problems occured.
   */
  int Initialize();

  //! Returns \c true if the preconditioner has been successfully computed.
  bool IsComputed() const {return(IsComputed_);};

  //! Computes the preconditioners.
  /*! \return
   * 0 if successful, 1 if problems occurred.
   */
  int Compute();

  //! Sets all the parameters for the preconditioner.            
  /*! Parameters currently supported:                 
   * - \c "amesos: solver type" : Specifies the solver type    
   *   for Amesos. Default: \c Amesos_Klu.        
   *                                             
   * The input list will be copied, then passed to the Amesos    
   * object through Amesos::SetParameters().       
   *                           
   * - \c "MST: forest number" : Specified the number of      
   *   times Kruskal's algorithm adds another forest to   
   *   the preconditioner              
   */
  virtual int SetParameters(Teuchos::ParameterList& parameterlis);

  //@}  


  //@{ \name Query methods.

  
  //! Returns a const reference to the internally stored matrix.
  const Epetra_RowMatrix& Matrix() const
  {
    return(*A_);
  }

  //! Returns the estimated conditioner number, computes it if necessary.
  /*!
   * not implemented
   */
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap, const int MaxIters = 1550,
		 const double Tol = 1e-9, Epetra_RowMatrix* Matrix_in = 0);

  //! Returns the computed condition number.
  double Condest() const
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

  //! Returns the total number of flops to compute the preconditioner.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the total number of flops to apply the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  //! Prints on ostream basic information about \c this object.
  virtual ostream& Print(std::ostream& os) const;


  //@}

 private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RefCountPtr<Epetra_RowMatrix> A_;
  
  //! Pointer to the matrix of the support graph.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> B_;

  //! Linear Problem required by Solver_.
  Teuchos::RefCountPtr<Epetra_LinearProblem> Problem_;
  
  //! Amesos solver, use to apply the inverse of the local matrix.
  Teuchos::RefCountPtr<Amesos_BaseSolver> Solver_;

  //! Contains a copy of the input parameter list.
  Teuchos::ParameterList List_;

  //! Epetra_Comm associated with matrix to be preconditioned.
  const Epetra_Comm & Comm_;

  //! If true, the preconditioner solves for the transpose of the matrix.
  bool UseTranspose_;

  //! Contains the estimated conditioner number.
  double Condest_;

  //! Contains the label of \c this object.
  string Label_;

  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;

  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;

  //! Contains the number of successful calls to Compute().
  int NumCompute_;

  //! Contains the number of successful calls to ApplyInverse().
  mutable int NumApplyInverse_;

  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;

  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;

  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;

  //! Contains the number of flops for Compute().
  double ComputeFlops_;

  //! Contain the number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;

  //! Time object.    
  Teuchos::RefCountPtr<Epetra_Time> Time_;

};

#endif /* IFPACK_SUPPORTGRAPH_H */
