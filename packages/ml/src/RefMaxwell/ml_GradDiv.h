/*!
 * \file ml_GradDiv.h
 *
 * \class GradDivPreconditioner
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

#ifndef ML_GRADDIV_H
#define ML_GRADDIV_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

// Some compilers think this name is too long...
#define GradDivPreconditioner GDP

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_common.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "ml_Preconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_ML_IFPACK
#include "Ifpack_Preconditioner.h"
#endif

namespace ML_Epetra
{


  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsGradDiv(Teuchos::ParameterList & inList,bool OverWrite=true);

  /*! The preconditioner(s) for the equations:
    - grad \lambda div B + \mu^{-1} B = f
    discretized using compatible discretizations (face elements).
    Since the preconditioners involve both a face and nodal component, different
    combinations of additive and/or multiplicative preconditioners can be used.  These options are
    controlled with the "refmaxwell: mode" option in the  Techos::ParameterList.
    The sublists "graddiv: 11list" and "graddiv: 22list" will be passed on
    to the edge and nodal sub-solvers appropriately.


    Details on this preconditioner can be found in Bochev, Siefert, Tuminaro, Xu and Zhu, 2007.
  */
  class GradDivPreconditioner: public virtual ML_Preconditioner
  {
  public:

    //@{ \name Constructors.

    //! Constructs a GradDivPreconditioner.
    // WARNING: All of these matrices will be shallow pointed to.  Please be
    //sure the matrices last until after RefMaxwellPreconditioner does.
    GradDivPreconditioner(const Epetra_CrsMatrix & K2_Matrix,                  // Face Grad-div + Mass
                          const Epetra_CrsMatrix & FaceNode_Matrix,            // Face-to-node interpolation matrix
                          const Epetra_CrsMatrix & D1_Clean_Matrix,            // Discrete curl w/o BCs
                          const Epetra_CrsMatrix & D0_Clean_Matrix,            // Discrete gradient w/o BCs
                          const Epetra_CrsMatrix & K0_Matrix,                  // Nodal Laplacian (for aggregation)
                          const Teuchos::ParameterList& List,
                          const bool ComputePrec = true);

    // Parameter-list driven version of main constructor
    // Paramters of relevance (below) should be RCP's to matrices:
    // FaceNode Face-to-node interpolation matrix
    // D0    - Discrete gradient w/o BCs
    // D1    - Discrete curl w/o BCs
    // K0    - Nodal Laplacian (for aggregation)
    GradDivPreconditioner(const Epetra_CrsMatrix& K2_Matrix,      // Face Grad-div + Mass
                          const Teuchos::ParameterList& List,
                          const bool ComputePrec = true);
    //@}


    //! @name Destructor
    //@{
    //! Destructor
    ~GradDivPreconditioner();
    //@}


    //@{ \name Mathematical functions.

    //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
    int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {
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
    const Epetra_Map& OperatorDomainMap() const {return(*DomainMap_);};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {return(*RangeMap_);};

    //! Return operator complexity and #nonzeros in fine grid matrix.
    void Complexities(double &complexity, double &fineNnz);
    //@}


  private:


    //@{ \name Internal data

    //! Matrix: Face Grad-div + Mass
    const Epetra_CrsMatrix * K2_Matrix_;
    //! Matrix: Edge Curl-curl
    Teuchos::RCP<Epetra_CrsMatrix> K1_Matrix_;
    //! Matrix: Face-to-node interpolation
    const Epetra_CrsMatrix * FaceNode_Matrix_;
    //! Matrix: Discrete curl w/ BCs
    Teuchos::RCP<Epetra_CrsMatrix> D1_Matrix_;
    //! Matrix: Discrete curl w/o BCs
    const Epetra_CrsMatrix * D1_Clean_Matrix_;
    //! Matrix: Discrete gradient w/ BCs
    Teuchos::RCP<Epetra_CrsMatrix> D0_Matrix_;
    //! Matrix: Discrete gradient w/o BCs
    const Epetra_CrsMatrix * D0_Clean_Matrix_;
    //! Nodal Laplacian (for aggregation)
    const Epetra_CrsMatrix * TMT_Matrix_;

#ifdef HAVE_ML_EPETRAEXT
    //! Structure for compatibility between Epetra and ML column maps.
    EpetraExt::CrsMatrix_SolverMap K2_Matrix_Trans_;
    EpetraExt::CrsMatrix_SolverMap K1_Matrix_Trans_;
    EpetraExt::CrsMatrix_SolverMap FaceNode_Matrix_Trans_;
    EpetraExt::CrsMatrix_SolverMap D0_Clean_Matrix_Trans_;
    EpetraExt::CrsMatrix_SolverMap D0_Matrix_Trans_;
    EpetraExt::CrsMatrix_SolverMap TMT_Matrix_Trans_;
#endif

    //! (1,1) Block Preconditioner
    ML_Preconditioner * FacePC;

    //! (2,2) Block Preconditioner
    ML_Preconditioner * EdgePC;

#ifdef HAVE_ML_IFPACK
    // Outer Smoother
    Epetra_Operator *IfSmoother;
#endif

    //! Domain Map
    const Epetra_Map* DomainMap_;
    //! Range Map
    const Epetra_Map* RangeMap_;
    //! Edge Map
    const Epetra_Map* EdgeMap_;

    //! Epetra communicator object
    const Epetra_Comm* Comm_;

    //! Verbosity flag
    bool verbose_;

    //! Extreme Verbosity flag
    mutable bool very_verbose_;

    //! Print hierarchy flag
    int print_hierarchy;
    //@}

#ifdef ML_TIMING
    //@{ \name Variables for Timing
    //! Number of applications
    mutable int NumApplications_;

    //! CPU time for all applications of the preconditioner
    mutable double ApplicationTime_;
    mutable bool FirstApplication_;
    //@ CPU time for first application
    mutable double FirstApplicationTime_;
    //! Number of construction phases
    int NumConstructions_;
    //! CPU time for construction of the preconditioner.
    double ConstructionTime_;
    //@}
#endif // ML_TIMING
  };// end GradDivPreconditioner
}//end namespace ML_Epetra

#endif
#endif
