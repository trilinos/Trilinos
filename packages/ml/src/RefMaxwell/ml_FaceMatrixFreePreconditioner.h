/*!
 * \file ml_FaceMatrixFreePreconditioner.h
 *
 * \class FaceMatrixFreePreconditioner
 *
 * \brief Matrix-Free preconditioning class for face problems.
 *
 * \date Last update to Doxygen: 2-Feb-10
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#ifndef ML_FACE_MATRIX_FREE_PRECONDITIONER_H
#define ML_FACE_MATRIX_FREE_PRECONDITIONER_H
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "ml_Preconditioner.h"

#ifdef HAVE_ML_IFPACK
#include "Ifpack_Chebyshev.h"
#endif

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#ifdef HAVE_ML_EPETRAEXT
#include "EpetraExt_SolverMap_CrsMatrix.h"
#endif

#ifdef HAVE_ML_IFPACK
class Ifpack_Chebyshev;
#endif


#include "ml_epetra.h"
#include "ml_mat_formats.h"

namespace ML_Epetra
{

  /*! Matrix-Free preconditioning class for edge Maxwell Problems. 
  */
  class FaceMatrixFreePreconditioner: public virtual ML_Preconditioner
  {
  public:
    //@{ \name Constructor
    //! Constructs an FaceMatrixFreePreconditioner.
    FaceMatrixFreePreconditioner(Teuchos::RCP<const Epetra_Operator> Operator, 
				 Teuchos::RCP<const Epetra_Vector> Diagonal,
				 Teuchos::RCP<const Epetra_CrsMatrix> FaceNode_Matrix,
                                 Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix,
                                 Teuchos::ArrayRCP<int> BCfaces,
                                 const Teuchos::ParameterList &List,
				 const bool ComputePrec = true);
    //@}
    
    
    //@{ 
    //! Destructor
    ~FaceMatrixFreePreconditioner();
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
      it is not possible to modify this pointer (other than creating a new preconditioner)  */


    //! Computes the preconditioner
    int ComputePreconditioner(const bool CheckFiltering = false);
    
    //! Recomputes the preconditioner
    int ReComputePreconditioner(){return(-1);}

    //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
      return(-1);}
    
    //! Apply the preconditioner to RHS B to get result X (X is also initial guess)
    int ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;
    
    //! Print the individual operators in the multigrid hierarchy.
    void Print(int whichHierarchy = -2);
    
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
    const Epetra_Map& OperatorDomainMap() const {return(*FaceDomainMap_);};
  
    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {return(*FaceRangeMap_);};
    
    //@}
    
  private:
    //@{ \name Internal functions

    //! Sets up the Chebyshev smoother 
    int SetupSmoother();
    
    //! Build pi nullspace described by Bochev, Siefert, Tuminaro, Xu and Zhu (2007).
    int BuildNullspace(Epetra_MultiVector *& nullspace);

    //! Build the face-to-node prolongator described by Bochev, Siefert, Tuminaro, Xu and Zhu (2007).
    int BuildProlongator();

    //! Forms the coarse matrix, given the build prolongator.
    int FormCoarseMatrix();

    //! Build the 1-unknown sparsity pattern
    int PBuildSparsity(ML_Operator *P, Epetra_CrsMatrix *&Psparse);

    //@}

    
    //@{ \name Internal data
    //! Dimension of space
    int dim;
    
    //! ML Communicator
    ML_Comm* ml_comm_;

    //! Fine-level operator
    Teuchos::RCP<const Epetra_Operator> Operator_;

    //! Matrix diagonal
    Teuchos::RCP<const Epetra_Vector> Diagonal_;

    //! Face-node indicence matrix.  Needed for aggregation
    Teuchos::RCP<const Epetra_CrsMatrix> FaceNode_Matrix_;

    //! TMT_Matrix.  Needed for nodal maps
    Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix_;    

    //! Dirichlet edges
    Teuchos::ArrayRCP<int> BCfaces_;
    
    //! Prolongator
    Epetra_CrsMatrix * Prolongator_;
#ifdef HAVE_ML_EPETRAEXT
    EpetraExt::CrsMatrix_SolverMap ProlongatorColMapTrans_;
#endif
    
    //! Inverse Diagonal
    Epetra_Vector * InvDiagonal_;

    //! Coarse Matrix
    Epetra_CrsMatrix * CoarseMatrix;
    ML_Operator * CoarseMat_ML;
    
    //! Level 2+ Preconditioner
    MultiLevelPreconditioner * CoarsePC; 

#ifdef HAVE_ML_IFPACK    
    //! Ifpack Chebyshev Smoother
    Epetra_Operator* Smoother_;
#endif
    
    //! Face Domain Map
    const Epetra_Map* FaceDomainMap_;
    //! Face Range Map
    const Epetra_Map* FaceRangeMap_;
    //! Epetra communicator object
    const Epetra_Comm* Comm_;


    //! Nodal Domain Map
    const Epetra_Map* NodeDomainMap_;
    //! Nodal Range Map
    const Epetra_Map* NodeRangeMap_;    
    //! Coarse Domain/Range Map
    Epetra_Map* CoarseMap_;
    
    //! Number of V-cycles to run
    int num_cycles;
    int MaxLevels;
    bool verbose_;
    bool very_verbose_;
    bool print_hierarchy;
    
    
    //@}  
  };//ML_FaceMatrixFreePreconditioner

}//end namespace ML_Epetra  


#endif

#endif
