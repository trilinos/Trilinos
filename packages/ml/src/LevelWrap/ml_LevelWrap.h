/*!
 * \file ml_LevelWrap.h
 *
 * \class LevelWrap
 *
 * \brief Class for wrapping a single level of a multilevel method.  Usually this is used for manual grid transfers on the finest level
 *
 * \date Last update to Doxygen: 31-Aug-11
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#ifndef ML_LEVELWRAP_H
#define ML_LEVELWRAP_H
#include "ml_common.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "ml_Preconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_RCP.hpp"
#include "Ifpack_Preconditioner.h"
#include "ml_MultiLevelPreconditioner.h"

namespace ML_Epetra
{

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsLevelWrap(Teuchos::ParameterList & inList,bool OverWrite=true);

  
  /*! A general preconditioner wrapper that allows a user to specify prolongators and
    potentially a coarse grid.  This can be used for adding a manual level on top of ML
    or (recursively) for a geometric multigrid algorithm, if you wanted.

  */
  class LevelWrap: public virtual ML_Preconditioner
  {
  public:
    //@{ \name Constructors.
    
    //! Constructs a Level Wrap (using R=P^T)
    LevelWrap(Teuchos::RCP<Epetra_CrsMatrix> A0,
	      Teuchos::RCP<Epetra_CrsMatrix> P0,
	      const Teuchos::ParameterList& List,
	      const bool ComputePrec = true);

    //! Constructs a Level Wrap (using different R)
    LevelWrap(Teuchos::RCP<Epetra_CrsMatrix> A0,
	      Teuchos::RCP<Epetra_CrsMatrix> P0,
	      Teuchos::RCP<Epetra_CrsMatrix> R0,
	      const Teuchos::ParameterList& List,
	      const bool ComputePrec = true);
    //@}
    

    //! @name Destructor
    //@{ 
    //! Destructor
    ~LevelWrap();
    //@}

    
    //@{ \name Mathematical functions.
    
    //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
      return(-1);}
    
    //! Apply the preconditioner w/ RHS B and get result X
    int ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;
    
    //@}

    
    //@{ \name Atribute access functions

    // Manually set the A1 operator
    //    void SetA1(Teuchos::RCP<Epetra_Operator> A1){A1_=A1;}

    //! Computes the preconditioner
    int ComputePreconditioner(const bool CheckFiltering = false);

    //! Recomputes the preconditioner
    int ReComputePreconditioner(){return(-1);}

    //! Print the individual operators in the multigrid hierarchy.
    void Print(int level);

    //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
    int DestroyPreconditioner();

    //! Sets use transpose (not implemented).
    int SetUseTranspose(bool UseTranspose){return(-1);}

    //! Returns the infinity norm (not implemented).
    double NormInf() const {return(0.0);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(false);};
  
    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const{return(false);};

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const{return(A0_->Comm());};
  
    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const {return(A0_->OperatorDomainMap());};
  
    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {return(A0_->OperatorRangeMap());};
    //@}


  private:
    
    //@{ \name Internal data    
    //!  A0 matrix
    Teuchos::RCP<const Epetra_CrsMatrix> A0_;

    //!  A1 matrix
    Teuchos::RCP<Epetra_CrsMatrix> A1_;

    //!  P0 matrix
    Teuchos::RCP<const Epetra_CrsMatrix> P0_;

    //!  R0 matrix
    Teuchos::RCP<const Epetra_CrsMatrix> R0_;

    //! Use P^T instead of R0.
    bool use_pt_;

    //! User provided A1
    bool user_A1_;

    //! Smoother
    Teuchos::RCP<Ifpack_Preconditioner> Smoother_;

    //! Smoother pre or post
    int pre_or_post;

    //! A1 preconditioner
    Teuchos::RCP<MultiLevelPreconditioner> A1prec_;

    //! Verbosity flag
    bool verbose_;

    //! Teuchos list
    Teuchos::ParameterList List_;
    //@}


    //@{ \name Variables for Timing
    //! Number of applications
    int NumApplications_;
    //! CPU time for all applications of the preconditioner
    mutable double ApplicationTime_;
    bool FirstApplication_;
    //@ CPU time for first application
    double FirstApplicationTime_;
    //! Number of construction phases
    int NumConstructions_;
    //! CPU time for construction of the preconditioner.
    double ConstructionTime_;
    //@}        
  };// end LevelWrap
}//end namespace ML_Epetra

#endif
#endif
