/*!
 * \file ml_EdgeMatrixFreePreconditioner.h
 *
 * \class EdgeMatrixFreePreconditioner
 *
 * \brief Matrix-Free preconditioning class for edge Maxwell Problems.
 *
 * \date Last update to Doxygen: 8-Feb-07
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_EDGE_MATRIX_FREE_PRECONDITIONER_H
#define ML_EDGE_MATRIX_FREE_PRECONDITIONER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "ml_Preconditioner.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "ml_RefMaxwell_Utils.h"

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

namespace ML_Epetra
{

  /*! Matrix-Free preconditioning class for edge Maxwell Problems.
  */
  class EdgeMatrixFreePreconditioner: public virtual ML_Preconditioner
  {
  public:
    //@{ \name Constructor
    //! Constructs an EdgeMatrixFreePreconditioner.
    EdgeMatrixFreePreconditioner(Teuchos::RCP<const Epetra_Operator> Operator,
				 Teuchos::RCP<const Epetra_Vector> Diagonal,
                                 Teuchos::RCP<const Epetra_CrsMatrix> D0_Matrix,
				 Teuchos::RCP<const Epetra_CrsMatrix> D0_Clean_Matrix,
                                 Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix,
                                 Teuchos::ArrayRCP<int> BCedges,
                                 const Teuchos::ParameterList &List,const bool ComputePrec = true);
    //@}


    //@{
    //! Destructor
    ~EdgeMatrixFreePreconditioner();
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
    int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {
      return(-1);}

    //! Apply the preconditioner to RHS B to get result X (X is also initial guess)
    int ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const;

    //! Print the individual operators in the multigrid hierarchy.
    void Print(int whichHierarchy = -2);

    //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
    int DestroyPreconditioner();

    //! Sets use transpose (not implemented).
    int SetUseTranspose(bool /* UseTranspose */){return(-1);}

    //! Returns the infinity norm (not implemented).
    double NormInf() const {return(0.0);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(false);};

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const{return(false);};

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const{return(*Comm_);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const {return(*EdgeDomainMap_);};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {return(*EdgeRangeMap_);};

    //! Return operator complexity and #nonzeros in fine grid matrix.
    void Complexities(double &complexity, double &fineNnz);

    //@}

  private:
    //@{ \name Internal functions

    //! Sets up the Chebyshev smoother
    int SetupSmoother();

    //! Build the edge-to-vector-node prolongator described in Bochev, Hu, Siefert and Tuminaro (2006).
    int BuildProlongator(const Epetra_MultiVector & nullspace);

    //! Forms the coarse matrix, given the prolongator
    int FormCoarseMatrix();


    //@}


    //@{ \name Internal data
    //! Dimension of space
    int dim;

    //! ML Communicator
    ML_Comm* ml_comm_;

    //! Fine-level operator
    Teuchos::RCP<const Epetra_Operator> Operator_;

    //! D0 or T matrix w/ dirichlet nodes and edges zero'd.  Used to generate prolongator.
    Teuchos::RCP<const Epetra_CrsMatrix> D0_Matrix_;

    //! D0 or T matrix w/ nothing zero'd.  Needed to generate the nullspace
    Teuchos::RCP<const Epetra_CrsMatrix> D0_Clean_Matrix_;

    //! TMT_Matrix.  Needed for nodal maps
    Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix_;

    //! Dirichlet edges
    Teuchos::ArrayRCP<int> BCedges_;

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

    //! Coarse data objects
    ML_Epetra::CoordPack CoarseCoord_;

    //! Level 2+ Preconditioner
    MultiLevelPreconditioner * CoarsePC;

    //! Coarse nullspace (slung in consecutive ML ordering)
    double *CoarseNullspace_;

#ifdef HAVE_ML_IFPACK
    //! Ifpack Chebyshev Smoother
    Epetra_Operator* Smoother_;
#endif

    //! Edge Domain Map
    const Epetra_Map* EdgeDomainMap_;
    //! Edge Range Map
    const Epetra_Map* EdgeRangeMap_;
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
  };//ML_EdgeMatrixFreePreconditioner

}//end namespace ML_Epetra


#endif

#endif
