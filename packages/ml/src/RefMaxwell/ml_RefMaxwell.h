/*!
 * \file ml_RefMaxwell.h
 *
 * \class RefMaxwellPreconditioner
 *
 * \brief Class for Reformulated Maxwell's Equations solvers.  Inherited from ML_Epetra_Operator
 *
 * \date Last update to Doxygen: 25-Jan-07
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#ifndef ML_REFMAXWELL_H
#define ML_REFMAXWELL_H

// Some compilers think this name is too long...
#define RefMaxwellPreconditioner ML_RMP

#define ENABLE_MS_MATRIX

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)  
#include "ml_common.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_Preconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_ML_EPETRAEXT
#include "EpetraExt_SolverMap_CrsMatrix.h"
#endif

namespace ML_Epetra
{

  int UpdateList(Teuchos::ParameterList &source, Teuchos::ParameterList &dest, bool OverWrite=true);
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsRefMaxwell(Teuchos::ParameterList & inList,bool OverWrite=true);


  
  /*! The preconditioner(s) for the Eddy Current Maxwell's equations using
    compatible discretizations (edge elements).  Since the preconditioners
    involve both an edge and nodal component, different combinations of additive
    and/or multiplicative preconditioners can be used.  These options are
    controlled with the "refmaxwell: mode" option in the  Techos::ParameterList.
    The sublists "refmaxwell: 11list" and "refmaxwell: 22list" will be passed on
    to the edge and nodal sub-solvers appropriately.

    
    Detail on this preconditioner can be found in Bochev, Hu, Siefert and
    Tuminaro, 2007.
  */
  class RefMaxwellPreconditioner: public virtual ML_Preconditioner
  {
  public:

    //@{ \name Constructors.
    
    //! Constructs a RefMaxwellPreconditioner.
    // WARNING: All of these matrices will be shallow pointed to.  Please be
    //sure the matrices last until after RefMaxwellPreconditioner does.
    RefMaxwellPreconditioner(const Epetra_CrsMatrix& SM_Matrix,      //S+M
                             const Epetra_CrsMatrix& D0_Clean_Matrix,//T or D0 w/ nothing zero'd
#ifdef ENABLE_MS_MATRIX
                             const Epetra_CrsMatrix& Ms_Matrix,      //M1(sigma)
#endif
                             const Epetra_CrsMatrix& M0inv_Matrix,   //M0^{-1}
                             const Epetra_CrsMatrix& M1_Matrix,      //M1(1)
                             //                             const Epetra_CrsMatrix& TMT_Matrix,     //T' M1(sigma) T
                             const Teuchos::ParameterList& List,
                             const bool ComputePrec = true);
    //@}
    

    //! @name Destructor
    //@{ 
    //! Destructor
    ~RefMaxwellPreconditioner();
    //@}

    
    //@{ \name Mathematical functions.
    
    //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
      return(-1);}
    
    //! Apply the preconditioner w/ RHS B and get result X
    int ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;
    
    //@}

    
    //@{ \name Attribute access functions

    //! Computes the multilevel hierarchy.
    /*! Computes the multilevel hierarchy. This function retrives the user's defines parameters (as
      specified in the input ParameterList), or takes default values otherwise, and creates the ML
      objects for aggregation and hierarchy. Allocated data can be freed used DestroyPreconditioner(),
      or by the destructor,
      
      In a Newton-type procedure, several linear systems have to be solved, Often, these systems
      are not too different. In this case, it might be convenient to keep the already 
      computed preconditioner (with hierarchy, coarse solver, smoothers), and use it to
      precondition the next linear system. ML offers a way to determine whether the 
      already available preconditioner is "good enough" for the next linear system. 
      The user should proceed as follows:
      - define \c "reuse: enable" == \c true
      - solve the first linear system. ML tries to estimate the rate of convergence, and record it;
      - change the values of the linear system matrix (but NOT its structure)
      - compute the new preconditioner as \c ComputePreconditioner(true)
      It is supposed that the pointer to the Epetra_RowMatrix remains constant. Currently,
      it is not possible to modify this pointer (other than creating a new preconditioner)
    */

    //! Computes the preconditioner
    int ComputePreconditioner(const bool CheckFiltering = false);

    //! Recomputes the preconditioner
    int ReComputePreconditioner(){return(-1);}

    //! Print the individual operators in the multigrid hierarchy.
    void Print(int whichHierarchy = 11);

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
    const Epetra_Comm& Comm() const{return(*Comm_);};
  
    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const {return(*DomainMap_);};
  
    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {return(*RangeMap_);};
    //@}


  private:

    //@{ \name Internal functions   

    //! Sets the Edge Smoother (if needed)
    int SetEdgeSmoother(Teuchos::ParameterList &List_);
    
    //! Implicitly applies in the inverse in a 2-1-2 format
    int ApplyInverse_Implicit_212(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;

    //! Implicitly applies in the inverse in an additive format
    int ApplyInverse_Implicit_Additive(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;

    //! Implicitly applies in the inverse in an 1-2-1 format
    int ApplyInverse_Implicit_121(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;


    //@}
    
    //@{ \name Internal data    
    
    //! Matrix: S+M1(sigma)
    const Epetra_CrsMatrix * SM_Matrix_;
    //! Matrix: T or D0 w/ dirichlet nodes and edges zero'd
    Epetra_CrsMatrix * D0_Matrix_;
    //! Matrix: D0 w/ nothing zero'd
    const Epetra_CrsMatrix * D0_Clean_Matrix_;      
#ifdef ENABLE_MS_MATRIX
    //! Matrix: M1(sigma)
    const Epetra_CrsMatrix * Ms_Matrix_;
#endif
    //! Matrix: M0^{-1}
    const Epetra_CrsMatrix * M0inv_Matrix_;
    //! Matrix: M1(1)
    Epetra_CrsMatrix * M1_Matrix_;
    //! Matrix: D0' M1(sigma) D0
    Epetra_CrsMatrix * TMT_Matrix_;
    //! Matrix: D0' M1(1) D0
    Epetra_CrsMatrix * TMT_Agg_Matrix_;


#ifdef HAVE_ML_EPETRAEXT
    //! Structure for compatibility between Epetra and ML column maps.
    EpetraExt::CrsMatrix_SolverMap SM_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.    
    EpetraExt::CrsMatrix_SolverMap D0_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.    
    EpetraExt::CrsMatrix_SolverMap D0_Clean_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.    
#ifdef ENABLE_MS_MATRIX
    EpetraExt::CrsMatrix_SolverMap Ms_Matrix_Trans_;
#endif
    //! Structure for compatibility between Epetra and ML column maps.
    EpetraExt::CrsMatrix_SolverMap M0inv_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.    
    EpetraExt::CrsMatrix_SolverMap M1_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.   
    EpetraExt::CrsMatrix_SolverMap TMT_Matrix_Trans_;
    //! Structure for compatibility between Epetra and ML column maps.
    EpetraExt::CrsMatrix_SolverMap TMT_Agg_Matrix_Trans_;    
#endif


    
    //! Dirichelt Edges
    int* BCrows; 
    int numBCrows;
    bool HasOnlyDirichletNodes;
    
    //! Vector: Diagonal of reformulated operator
    Epetra_Vector* Diagonal_;
    //! (1,1) Block Operator
    Teuchos::RCP<Epetra_Operator> Operator11_;
    
    //! (1,1) Block Preconditioner
    ML_Preconditioner * EdgePC;
    //! (2,2) Block Preconditioner
    //    ML_Preconditioner * NodePC;
    MultiLevelPreconditioner * NodePC; // This is a HAQ

    //! Outer Edge Smoother(s)
    MultiLevelPreconditioner *PreEdgeSmoother;
    MultiLevelPreconditioner *PostEdgeSmoother;
#ifdef HAVE_ML_IFPACK
    Epetra_Operator *IfSmoother;
#endif    

    //! Solver mode
    string mode;

    //! Aggregation info
    bool aggregate_with_sigma;
   
    //! Mass lumping
    bool lump_m1;   

    // EXPERIMENTAL: Local nodal subsolver
    bool use_local_nodal_solver;
    MultiLevelPreconditioner *LocalNodalSolver;
    const Epetra_CrsMatrix* NodesToLocalNodes;
    Epetra_CrsMatrix* LocalNodalMatrix;
    
    //! Domain Map
    const Epetra_Map* DomainMap_;
    //! Range Map
    const Epetra_Map* RangeMap_;
    //! Node Map
    const Epetra_Map* NodeMap_;

    //! Epetra communicator object
    const Epetra_Comm* Comm_;

    //! Verbosity flag
    bool verbose_;

    //! Extreme Verbosity flag
    mutable bool very_verbose_;    //HAQ

    //! Print hierarchy flag
    int print_hierarchy;
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
  };// end RefMaxwellPreconditioner
}//end namespace ML_Epetra

#endif
#endif
