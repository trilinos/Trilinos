/*!
 * \file ml_MultiLevelPreconditioner.h
 *
 * \class MultiLevelPreconditioner
 *
 * \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#ifndef ML_MULTILEVELPRECONDITIONER_H
#define ML_MULTILEVELPRECONDITIONER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
// define the following to allow compilation without AztecOO
#ifndef HAVE_ML_AZTECOO
#ifndef AZ_PROC_SIZE
#define AZ_PROC_SIZE 1
#endif
#ifndef AZ_OPTIONS_SIZE
#define AZ_OPTIONS_SIZE 1
#endif
#ifndef AZ_PARAMS_SIZE
#define AZ_PARAMS_SIZE 1
#endif
#ifndef AZ_STATUS_SIZE
#define AZ_STATUS_SIZE 1
#endif
#endif

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_FECrsMatrix;
class Epetra_VbrMatrix;

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

#define ML_MEM_SIZE      20
#define ML_MEM_INITIAL    0
#define ML_MEM_FINAL      1
#define ML_MEM_SMOOTHER   2
#define ML_MEM_COARSE     3
#define ML_MEM_HIERARCHY  4
#define ML_MEM_PREC_FIRST 5
#define ML_MEM_PREC_OTHER 6
#define ML_MEM_TOT1       7
#define ML_MEM_TOT2       8
#define ML_MEM_INITIAL_MALLOC    10
#define ML_MEM_FINAL_MALLOC      11
#define ML_MEM_SMOOTHER_MALLOC   12
#define ML_MEM_COARSE_MALLOC     13
#define ML_MEM_HIERARCHY_MALLOC  14
#define ML_MEM_PREC_FIRST_MALLOC 15
#define ML_MEM_PREC_OTHER_MALLOC 16
#define ML_MEM_TOT1_MALLOC       17
#define ML_MEM_TOT2_MALLOC       18

#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef HAVE_ML_AZTECOO
#include "Epetra_MultiVector.h"
#include "Epetra_MsrMatrix.h"
#endif
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_ML_EPETRAEXT
#include "EpetraExt_SolverMap_CrsMatrix.h"
#endif
#include "ml_epetra_utils.h"

namespace ML_Epetra
{

  //! Sets default parameters for aggregation-based preconditioners.
  /*! This function, defined in the namespace ML_Epetra, can be used to set
    default values in a user's defined Teuchos::ParameterList.
    \param ProblemType (In) : a std::string, whose possible values are:
       - "SA" : classical smoothed aggregation preconditioners;
       - "NSSA" : default values for Petrov-Galerkin preconditioner for nonsymmetric systems
       - "maxwell" : default values for aggregation preconditioner for eddy current systems
       - "DD" : defaults for 2-level domain decomposition preconditioners based
       on aggregation;
       - "DD-LU" : Like "DD", but use exact LU decompositions on each subdomain;
       - "DD-ML" : 3-level domain decomposition preconditioners, with coarser
       spaces defined by aggregation;
      - "DD-ML-LU" : Like "DD-ML", but with LU decompositions on each subdomain.
    \param List (Out) : list which will populated by the default parameters
    \param options (In/Out) : integer array, of size \c AZ_OPTIONS_SIZE, that will be
    populated with suitable values. A pointer to \c options will be stick into
    the parameters list. Note that this array is still required to apply the
    preconditioner! Do not delete options, nor let it go out of scope. The default value is
    0, meaning that \c SetDefaults() will allocate the array.
    \param params (In/Out) : double array, of size \c AZ_PARAMS_SIZE. See comments
    for \c options.
    \param OverWrite (In) : boolean.  If false, any pre-existing values in the
    parameter list will be preserved.  Default value is true, i.e., any
    pre-existing values may be overwritten.
   */
  int SetDefaults(std::string ProblemType, Teuchos::ParameterList & List,
		  int * options = 0, double * params = 0, const bool OverWrite=true);

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsDD(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners, using LU on each subdomain
  int SetDefaultsDD_LU(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.
  int SetDefaultsDD_3Levels(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners with LU.
  int SetDefaultsDD_3Levels_LU(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets default parameters for the eddy current equations equations.
  int SetDefaultsMaxwell(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets default parameters for classical smoothed aggregation.
  int SetDefaultsSA(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets defaults for energy minimization preconditioning for nonsymmetric problems.
  int SetDefaultsNSSA(Teuchos::ParameterList & List,
		               Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Sets defaults for classical amg
  int SetDefaultsClassicalAMG(Teuchos::ParameterList & List,
                      Teuchos::RCP<std::vector<int> > &options,
                       Teuchos::RCP<std::vector<double> > &params,
                       bool Overwrite=true);

  //! Reads in parameter list options from file.
  int ReadXML(const std::string &FileName, Teuchos::ParameterList &List,
                   const Epetra_Comm &Comm);

  //! Enumerated type indicating the type of AMG solver to be used.

  enum AMGType {

  ML_SA_FAMILY,         /*< Smoothed aggregation solver including EMIN,NSA */
  ML_MAXWELL,           /*< Old Maxwell solver */
  ML_COMPOSITE,         /*< Composite AMG: block diagonal prolongator */
  ML_CLASSICAL_FAMILY   /*< Classical AMG */
  };

/*!

   \brief MultiLevelPreconditioner: a class to define black-box multilevel preconditioners using aggregation methods.

   Class ML_Epetra::MultiLevelPreconditioner defined black-box algebraic
   multilevel preconditioners of matrices defined as Epetra_RowMatrix derived
   objects. The resulting preconditioner can be used in AztecOO, and in any
   other solver that accepts Epetra_Operator derived objects, and apply the
   action of the given Epetra_Operator using ApplyInverse().

   Please refer to the user's guide for a detailed introduction to
   this class, examples, and description of input parameters.

    This file requires ML to be configured with the following options:
    - \c --enable-epetra
    - \c --enable-teuchos

    The following option is suggested:
    - \c --enable-amesos
    - \c --enable-ifpack

    Some part of this class needs the following options:
    - \c --enable-aztecoo
    - \c --enable-anasazi

    It is important to note that ML is more restrictive than Epetra for
    the definition of maps. It is required that RowMatrixRowMap() is equal
    to OperatorRangeMap(). This is because ML needs to perform matrix-std::vector
    product, as well as getrow() functions, on the same data distribution.

    Also, for square matrices, OperatorDomainMap() must be as
    OperatorRangeMap().

    Several examples are provided in the \c examples subdirectories:
    - \ref ml_preconditioner_cpp is an introductory
      example;
    - \ref ml_2level_DD_cpp shows how to
      define a 2-level domain decomposition preconditioner using
      this class;
    - \ref ml_viz_cpp details how to visualize the aggregates;
    - \ref ml_maxwell_cpp reports how to
      use this class for Maxwell problems.

   \note
   Namespace ML_Epetra contains another Epetra_Operator derived class,
   ML_Epetra::MultiLevelOperator.
   - you should use MultiLevelOperator
     when your code already defines the required ML objects, with the optimal
     choice of parameters, and you just want to wrap the already defined ML
     preconditioners for AztecOO problems;
   - you should use MultiLevelPreconditioner
     when you have an Epetra_RowMatrix, and you don't want to code the
     conversion to ML_Operator, the creation of the hierarchy and the
     aggregates, and/or you want to experiment various combinations of the
     parameters, simply changing some parameters in a Teuchos::ParameterList.

   Defaults parameters can be specified using function SetDefaults().

    \author Marzio Sala, SNL 9214
*/
class MultiLevelPreconditioner : public virtual Epetra_Operator {

public:

  //@{ \name Constructors.

  //! Constructs a MultiLevelPreconditioner with default values.

  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
                           const bool ComputePrec = true);

  //! Constructs a MultiLevelPreconditioner. Retrieves parameters from \c List.

  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
			   const Teuchos::ParameterList & List,
			   const bool ComputePrec = true);

  //! Constructs a MultiLevelPreconditioner from an ML_Operator. Retrieves parameters from \c List.

  MultiLevelPreconditioner(ML_Operator* Operator,
			   const Teuchos::ParameterList& List,
			   const bool ComputePrec = true);

  //! Constructs a MultiLevelPreconditioner which is actually a composite AMG hierarchy using an array of ML_Operator's and an array of parameter lists.

  MultiLevelPreconditioner(ML_Operator *Operator,
                           const Teuchos::ParameterList& List,
                           Epetra_RowMatrix **DiagOperators,
			   Teuchos::ParameterList *DiagLists,
                           int NBlocks = 1,
			   const bool ComputePrec = true);

  //! \brief MultiLevelPreconditioner constructor for Maxwell's equations.
  /*! Takes the stiffness and mass terms of the matrix combined.

      \param EdgeMatrix - (In) Linear matrix to be solved.
      \param GradMatrix - (In) Node-to-edge connectivity matrix, a.k.a,
                               topological gradient
      \param NodeMatrix - (In) Auxiliary nodal finite element matrix
      \param List - (In) Teuchos parameter list containing solver options.
      \param ComputePrec - (In) Optional argument that specifies whether to
                                create preconditioner immediately.
                                Default is true.
      \param UseNodeMatrixForSmoother - (In) Use the nodal matrix for the nodal
                            portion of the Hipmair smoother (if used).
  */

  MultiLevelPreconditioner(const Epetra_RowMatrix& EdgeMatrix,
			   const Epetra_RowMatrix& GradMatrix,
			   const Epetra_RowMatrix& NodeMatrix,
			   const Teuchos::ParameterList& List,
			   const bool ComputePrec = true,
                           const bool UseNodeMatrixForSmoother = false);

  //! \brief MultiLevelPreconditioner constructor for Maxwell's equations.
  /*! Takes the stiffness and mass terms of the matrix separately.

      \param CurlCurlMatrix - (In) The curl-curl (stiffness) term of the
                                   matrix to be solved.
      \param MassMatrix - (In) The mass term of the matrix to be solved.
      \param GradMatrix - (In) Node-to-edge connectivity matrix, a.k.a,
                               topological gradient
      \param NodeMatrix - (In) Auxiliary nodal finite element matrix
      \param List - (In) Teuchos parameter list containing solver options.
      \param ComputePrec - (In) Optional argument that specifies whether to
                                create preconditioner immediately.
                                Default is true.
  */

  MultiLevelPreconditioner(const Epetra_RowMatrix & CurlCurlMatrix,
             const Epetra_RowMatrix & MassMatrix,
             const Epetra_RowMatrix & TMatrix,
             const Epetra_RowMatrix & NodeMatrix,
             const Teuchos::ParameterList & List,
             const bool ComputePrec = true);

#define NewStuff
#ifdef NewStuff
  //! Constructs a MultiLevelPreconditioner for multiphysics with variable dofs per node

  MultiLevelPreconditioner(Epetra_RowMatrix & RowMatrix,
			   const Teuchos::ParameterList & List,
                           const int & nNodes,
                           const int & maxDofPerNode,
                           bool * dofPresent,
                           Epetra_Vector & Lhs, 
                           Epetra_Vector & Rhs, 
                           const bool  rhsAndsolProvided,
			   const bool ComputePrec = true);
// ================================================ ====== ==== ==== == =
/*! Constructor for scalar PDE problems based on applying AMG to the distance
 *  Laplacian operator when constructing grid transfers. The main unique
 *  feature is that there may be some dofs that correspond to the same node
 *  location. These shared dofs fall into two categories. If these dofs are
 *  strongly connected to each other (as determined by tol), they are
 *  explicitly elminated from the Laplacian (merged into a supernode). Once
 *  a P is obtained, this P is then expanded to account for shared nodes
 *  by simply duplicating the supernodes row of P for each of the individual
 *  vertices that contribute to the supernode. If share dofs are weakly
 *  connected (or not connected at all), nothing special is done (other than
 *  the ususal ignoring of weak connections). One last trick is employed, 
 *  connections between supernodes and non-supernodes (i.e., regular nodes)
 *  are always assumed to be weak. Shared nodes are often used to capture
 *  interfaces or other features. By breaking these connections, the AMG
 *  can better maintain these features throughout the hierarchy. Finally, the
 *  constructor also allows for *  the removal of column nonzeros associated
 *  with Dirichlet points. To use this option the rhs and initial guess must be
 *  provided. Modification of the matrix, rhs, and initial guess must be 
 *  allowable to use this option.
 */

  MultiLevelPreconditioner(Epetra_RowMatrix & RowMatrix,
    const Teuchos::ParameterList & List,
    const double distTol, // two points are at the same location when
                          // || (x_1,y_1,z_1) -  (x_2,y_2,z_2)||_2 < distTol
    const double tol,     // ignore values when
                          //       A(i,j)^2 < A(i,i)*A(j,j)*tol^2
    Epetra_Vector & Lhs,
    Epetra_Vector & Rhs,
    const bool  rhsAndsolProvided,
    const bool ComputePrec = true);
#endif
#ifdef HAVE_ML_AZTECOO
  //! MultiLevelPreconditioner constructor for Maxwell's equations.
  /*! Takes the stiffness and mass terms of the matrix combined.  The edge
      matrix is of type Epetra_Msr, a light-weight wrapper for old-style Aztec
      MSR matrices.  This is intended as transition code for Aztec users.

      \param EdgeMatrix - (In) Linear matrix to be solved.
      \param GradMatrix - (In) Node-to-edge connectivity matrix, a.k.a,
                               topological gradient
      \param NodeMatrix - (In) Auxiliary nodal finite element matrix
      \param proc_config - (In) Aztec array specifying processor layout.
      \param List - (In) Teuchos parameter list containing solver options.
      \param ComputePrec - (In) Optional argument that specifies whether to
                                create preconditioner immediately.
                                Default is true.
  */

  MultiLevelPreconditioner(const Epetra_MsrMatrix & EdgeMatrix,
                         ML_Operator * GradMatrix,
                         AZ_MATRIX * NodeMatrix,
                         int       * proc_config,
                         const Teuchos::ParameterList & List,
                         const bool ComputePrec = true);
#endif

  //@}

  //@{ \name Destructor.

  //! Destroys the preconditioner.
  virtual ~MultiLevelPreconditioner() {
    if (IsComputePreconditionerOK_)
      DestroyPreconditioner();
  }

  //@}

  //@{ \name Query functions

  //! Prints label associated to this object.
  const char* Label() const{return(Label_);};

  //! Prints unused parameters in the input ParameterList on standard output.
  void PrintUnused() const
  {
    List_.unused(std::cout);
  }

  //! Prints unused parameters in the input ParameterList on the specified stream.
  void PrintUnused(std::ostream & os) const
  {
    List_.unused(os);
  }

  //! Prints unused parameters in the input ParameterList to std::cout on proc \c MyPID.
  /*! Mispelled parameters are simply ignored. Therefore, it is often the best
   * choice to print out the parameters that have not been used in the
   * construction phase.
   * - \param MyPID (In) : ID of process that should print the unused parameters.
   */
  void PrintUnused(const int MyPID) const;

  //! Gets a reference to the internally stored parameters' list.
  Teuchos::ParameterList& GetList()
  {
    return List_;
  }

  // Get a copy of the internally stored output list.
  Teuchos::ParameterList GetOutputList()
  {
    return OutputList_;
  }

  //! Prints on \c std::cout the values of the internally stored parameter list
  void PrintList();

  //! Copies \c List into the internally stored parameter list object.
  int SetParameterList(const Teuchos::ParameterList& List);

  //@}

  //@{ \name Mathematical functions.

  //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
  int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {
    return(-1);}

  //! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

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

  int ComputePreconditioner(const bool CheckFiltering = false);

  /*! @brief Recompute the preconditioner (not implemented for Maxwell).

    @param keepFineLevelSmoother (In) : If true, the fine level
    smoother is not recomputed.  This is useful if the smoother is
    expensive to create, e.g., an incomplete factorization, and the
    fine level matrix has not changed.
  */
  int ReComputePreconditioner(bool keepFineLevelSmoother=false);

  //! Print the individual operators in the multigrid hierarchy.
  void Print(int level = -2);

  int ComputeAdaptivePreconditioner(int TentativeNullSpaceSize,
                                    double* TentativeNullSpace);

  //! Queries whether multilevel hierarchy has been computed or not.
  int IsPreconditionerComputed()  const
  {
    return(IsComputePreconditionerOK_);
  }

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool /* useTranspose */){return(-1);}

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

  //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
  int DestroyPreconditioner();

  //! Returns a reference to the internally stored RowMatrix.
  const Epetra_RowMatrix& RowMatrix() const
  {
    return(*RowMatrix_);
  }

  //! Returns a reference to RowMatrix->Map().
  const Epetra_BlockMap& Map() const
  {
    return(RowMatrix_->Map());
  }

  //! Returns the global number of rows in the matrix.
  int NumGlobalRows() const
  {
    return(RowMatrix_->NumGlobalRows());
  }

  //! Returns the global number of columns in the matrix.
  int NumGlobalCols() const
  {
    return(RowMatrix_->NumGlobalCols());
  }

  //! Returns the local number of rows in the matrix.
  int NumMyRows() const
  {
    return(RowMatrix_->NumMyRows());
  }

  //! Returns the local number of columns in the matrix.
  int NumMyCols() const
  {
    return(RowMatrix_->NumMyCols());
  }

  //! Prints the computational stencil for the specified row and equation (for 2D Cartesian grids only)
  /*! For problems defined on 2D Cartesian grids (with node numbering increasing
   * along the x-axis), this function prints out the stencil in an intelligible
   * form.
   * \param nx (In) : number of nodes along the X-axis
   * \param ny (In) : number of nodes along the Y-axis
   * \param NodeID (In) : (local) ID of node that will be used to print the
   *   stencil. If set to -1, the code will automatically chose an internal node.
   *   Default: -1.
   * \param EquationID (In) : ID of the equation that will be used to print the
   *   stencil (default = 0)
   */
  int PrintStencil2D(const int nx, const int ny,
		     int NodeID = -1,
		     const int EquationID = 0);

  //! Cheap analysis of each level matrix.
  int AnalyzeHierarchy(const bool AnalyzeMatrices,
                       const int PreCycles, const int PostCycles,
                       const int MLCycles);

  //! Analyze the effect of each level's smoother on a random std::vector.
  int AnalyzeSmoothers(const int NumPreCycles = 1,
                       const int NumPostCycles = 1);

  //! Analyze the effect of the coarse solver on a random std::vector.
  int AnalyzeCoarse();

  //! Analyze the effect of the ML cycle on a random std::vector.
  int AnalyzeCycle(const int NumCycles = 1);

  //! Test several smoothers on fine-level matrix.
  int TestSmoothers(Teuchos::ParameterList& InputList,
		    const bool IsSymmetric = false);

  //! Test several smoothers on fine-level matrix using the current parameters.
  int TestSmoothers(const bool IsSymmetric = false) {
    return(TestSmoothers(List_,IsSymmetric));
  }

  //! Returns a pointer to the internally stored ml pointer
  const ML* GetML(const int WhichML= -1) const
  {
    if (WhichML < 0)
      return ml_;
    else if (WhichML == 0)
      return ml_nodes_;
    else
      return(0);
  }

  //! Returns a pointer to the internally stored agg pointer
  const ML_Aggregate* GetML_Aggregate() const
  {
    return agg_;
  }

  //! Generic interface to visualization methods.
  int Visualize(bool VizAggre, bool VizPreSmoother,
		bool VizPostSmoother, bool VizCycle,
		int NumApplPreSmoother, int NumApplPostSmoother,
		int NumCycleSmoother);

  //! Visualizes the shape of the aggregates.
  int VisualizeAggregates();

  //! Visualizes the effect of smoothers on a random std::vector.
  int VisualizeSmoothers(int NumPrecCycles = 1,
			 int NumPostCycles = 1);

  //! Visualizes the effect of the ML cycle on a random std::vector.
  int VisualizeCycle(int NumCycles = 1);

  /*! Creates label for this object (printed out by AztecOO).  This does not
      allocate/reallocate any memory.
  */
  int CreateLabel();

  void ReportTime();

  //! Return operator complexity and #nonzeros in fine grid matrix.
  void Complexities(double &complexity, double &fineNnz);

//@}

private:

  //! Copy constructor (NOT DEFINED)
  MultiLevelPreconditioner(const MultiLevelPreconditioner & /* rhs */)
  {};

  //! operator = (NOT DEFINED)
  MultiLevelPreconditioner & operator = (const MultiLevelPreconditioner & /* rhs */)
  {
    return *this;
  };

  //@{ \name Internal setting functions
  //! Initializes object with defauls values.
  int Initialize();

  /*! Sets smoothers.
    @param[in] skipFineLevelSmoother : If true, the fine level smoother is not set.  This is intended to be used in
    combination with ReComputePreconditioner.
  */
  int SetSmoothers(bool skipFineLevelSmoother=false);

  //! Sets coarse level solvers.
  int SetCoarse();

  //! Sets aggregation schemes.
  int SetAggregation();

  //! Sets preconditioner type (usually, V-cycle).
  int SetPreconditioner();

  //! Sets the null space for non-Maxwell problems.
  int SetNullSpace();

  //! Checks correctness of null space (discrete gradient) for Maxwell problems.
  //! The curl-curl and mass matrices must be supplied separately.
  void CheckNullSpace();

  //! Applies boundary conditions to gradient matrix.  (Maxwell's equations)
  void Apply_BCsToGradient( const Epetra_RowMatrix & EdgeMatrix,
                            const Epetra_RowMatrix & T);

  //! Transforms Epetra matrix column map (if necessary) to be compatible with
  /*! how ML handles column indices.  Any matrix that cannot be dynamically
      cast to an Epetra_CrsMatrix will not be changed.

      \param A - (In) Matrix that is to be transformed.
      \param transform - (In) EpetraExt widget that does the transformation.
      \param matrixName - (In) Optional label for the incoming matrix.
   */

#ifdef HAVE_ML_EPETRAEXT
  Epetra_RowMatrix* ModifyEpetraMatrixColMap( const Epetra_RowMatrix &A,
                                   EpetraExt::CrsMatrix_SolverMap &transform,
                                   const char* matrixName );
#endif
  //! Destroys Preconditioner if it not needed anymore. This includes some 'filtering' checks.

  int ConditionallyDestroyPreconditioner(const bool CheckPreconditioner);

  //! Set the finest level matrix in the MG hierarchy

  int SetFinestLevelMatrix();

  //! Set pointers indicating correspondence between array entries and MG levels

  int SetLevelIds(int Direction);

  //! Set eigenvalue scheme to be used by ML for spectral radius

  int SetEigenScheme();

  //! Dump various output matrices for debugging

  int MatrixDumper();

  //! Recompute complexities and print them.

  int ComputeAndPrintComplexities();

  //! Sets prolongator smoother parameters.
  int SetSmoothingDamping();

  //! Sets damping parameter for classical smoothed aggregation.
  int SetSmoothingDampingClassic();

#define OLD_AUX
#ifdef OLD_AUX
  int CreateAuxiliaryMatrixCrs(Epetra_FECrsMatrix * & FakeMatrix);

  int CreateAuxiliaryMatrixVbr(Epetra_VbrMatrix * & FakeMatrix);
#endif

  int SetupCoordinates();

  void PrintMem(char *fmt, int size, int, int);

  void PrintMemoryUsage();

  int SetFiltering();

  void RandomAndZero(double *, double *, int);

  //! Checks whether the previously computed preconditioner is still valuable for the newly available linear system.
  /*! Used only when \c "reuse: enable" is \c true, and
   * ComputePreconditioner(true) is called. */
  bool CheckPreconditionerKrylov();

  void VectorNorms(double*, int, double*,double*);

  //@}

  //@{ \name Internal data

  //! Pointer to ML_Struct
  ML* ml_;
  //! ML communicator, convenient to have separately from ml_,
  //  ml_nodes_, one or all of which may be null.
  ML_Comm* ml_comm_;

  //! indicates the type of AMG solver to be used: ML_SA_FAMILY, ML_MAXWELL, ML_COMPOSITE
  AMGType AMGSolver_;

  //! ML_Aggregate, contains aggregate information
  ML_Aggregate* agg_;

  //! ML_AMG, contains amg family informatin
  ML_AMG* amg_;

  //! Label for this object
  char* Label_;
  //! User-provided label for identifying preconditioner ctor/dtor, in the case
  //  of multiple instances of ML_Epetra::MultiLevelPreconditioner.
  std::string mlpLabel_;

  //! pointer to linear system matrix
  const Epetra_RowMatrix* RowMatrix_;

  //! AfineML_ points to the original ML operator passed in to the block
  //  matrix/composite version of the constructor.
  ML_Operator *AfineML_;

  //! Multigrid hierarchies applied to submatrices and used in a composite
  //  form to define the overall AMG hierarchy

  ML_Epetra::MultiLevelPreconditioner **SubMatMLPrec_;

  //! specifies whether a hierarchy already exists or not.
  bool IsComputePreconditionerOK_;

  //! Number of levels
  int NumLevels_;
  //! Domain Map
  const Epetra_Map* DomainMap_;
  //! Range Map
  const Epetra_Map* RangeMap_;
  //! Epetra communicator object
  const Epetra_Comm* Comm_;
  bool  ownership_;
  //! proc_config for Aztec smoothers
  int   ProcConfig_[AZ_PROC_SIZE];
  //! options for Aztec smoothers
  Teuchos::RCP<std::vector<int> >    SmootherOptions_;
  //! params for Aztec smoothers
  Teuchos::RCP<std::vector<double> > SmootherParams_;
  //! status for Aztec smoothers
  double SmootherStatus_[AZ_STATUS_SIZE];

  //! List containing all input parameters.
  Teuchos::ParameterList List_;
  //! List containing all output parameters
  Teuchos::ParameterList OutputList_;

  //! Maximum number of levels
  int MaxLevels_;

  //! Number of applications of the ML cycle
  int CycleApplications_;

  //! If \c true, zero starting solution is used in the application of the cycle.
  bool ZeroStartingSolution_;

  //! Integer array used to easily handle ML_INCREASING and ML_DECREASING
  /*! Integer array, of size MaxLevels_, that contain the ML level ID
    for the first logical level, and so on for all levels. The ML level ID
    of logical level L is LevelID_[L].
    In this interface, all levels move from 0 to MaxLevels-1.
    ML's level for interface's level i is LevelID_[i]
  */
  std::vector<int> LevelID_;

  //! If not NULL, contains the allocated null space std::vector
  double* NullSpaceToFree_;

  //! all std::cout's have this prefix (default'd in Initialize() )
  std::string PrintMsg_;
  //! all std::cerr's have this prefix (default'd in Initialize() )
  char ErrorMsg_[80];
  //! true if information has to be printed on this process
  bool verbose_;
  //! Number of PDE equations.
  int NumPDEEqns_;
  //! Number of iterations to use in profiling
  int profileIterations_;

  //@}

  //@{ \name Composite AMG variables

  //! Number of blocks making up composite operator

  int NBlocks_;

  //! Array of Diagonal Operators

  Epetra_RowMatrix **DiagOperators_;

  //! Array of Parameter lists for diagonal operators

  Teuchos::ParameterList *DiagLists_;

  //! special flag (currently for variable dof per node multiphysics AMG) skips smoother when building Laplacian hierarchy. 
  bool DontSetSmoothers_;
  //@}

  //@{ \name Maxwell variables

  //! Main matrix for Maxwell
  const Epetra_RowMatrix* EdgeMatrix_;
  //! stiffness and mass matrices
  const Epetra_RowMatrix* CurlCurlMatrix_;
  //! true if we summed curl-curl and mass
  bool CreatedEdgeMatrix_;
  const Epetra_RowMatrix* MassMatrix_;
  //! aux matrix for Maxwell
  const Epetra_RowMatrix* NodeMatrix_;
  //! T^T A T Matrix for use with Maxwell
  ML_Operator* TtATMatrixML_;
  bool UseNodeMatrixForSmoother_;
  bool CreatedNodeMatrix_;
  //! Auxiliary matrix used in intermediate step
  ML_Operator* ML_Kn_;
  bool CreatedML_Kn_;
  //! T matrix for Maxwell
  const Epetra_RowMatrix* TMatrix_;
#ifdef HAVE_ML_EPETRAEXT
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap RowMatrixColMapTrans_;
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap NodeMatrixColMapTrans_;
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap TMatrixColMapTrans_;
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap CurlCurlMatrixColMapTrans_;
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap MassMatrixColMapTrans_;
  //! Structure for compatibility between Epetra and ML column maps.
  EpetraExt::CrsMatrix_SolverMap TtATMatrixColMapTrans_;
#endif
  bool CreatedTMatrix_;
  ML_Operator* TMatrixML_;
  ML_Operator* TMatrixTransposeML_;
  ML_Operator** Tmat_array, ** Tmat_trans_array;
  ML_Operator** MassMatrix_array; // If curlcurl & mass are separate
  ML_Operator** CurlCurlMatrix_array;  // If curlcurl & mass are separate
  //! Auxiliary ML structure for Maxwell's equations.
  ML* ml_nodes_;

  void** nodal_args_,** edge_args_;

  //@}

  //@{ \name Variables for Timing
  //! Number of applications
  int NumApplications_;
  //! CPU time for all applications of the preconditioner
  double ApplicationTime_;
  bool FirstApplication_;
  //@ CPU time for first application
  double FirstApplicationTime_;
  //! Number of construction phases
  int NumConstructions_;
  //! CPU time for construction of the preconditioner.
  double ConstructionTime_;

  //@}

  // other stuff for old ML's compatibility
  Epetra_CrsMatrix* RowMatrixAllocated_;
  bool  AllocatedRowMatrix_; // used for composite constructor only

  bool AnalyzeMemory_;

  int memory_[ML_MEM_SIZE];

  // filtering stuff
  std::vector<double> flt_NullSpace_;
  ML* flt_ml_;
  ML_Aggregate* flt_agg_;

  // for reuse of preconditioning
  double RateOfConvergence_;

}; // class MultiLevelPreconditioner

} // namespace ML_Epetra

#ifdef NewStuff
/******************************************************************************
  Vector class that is used for for multiphysics multigrid on matrices with
  varaible dofs per node. Didn't want anything specific to epetra, ML, and 
  MueLu; didn't like new/malloc features of std::vector, and I have a hard
  time figuring out teuchos stuff ... so I implemented my own. The syntax
  is pretty similar to std::vector, but it is more efficient from a memory
  allocation perspective (though less flexible). 
******************************************************************************/

enum free_s {useFree, useDelete, neither};
template <class T> class MLVec {
  public:

    // different constructors & destructors

    MLVec();              
    MLVec(int mysize);
    MLVec(T *startptr, T *endptr, bool okToChangePtr = false, free_s freeType = neither);
    ~MLVec();

    int     size() const;        // return vector size
    T       *getptr();           // return pointer to start of raw vector data
    void    wrap(T *startptr, T *endptr,// same as above construct but allows
            bool okToChangePtr = false, // one to first create an empty MLVec
            free_s freeType = neither); // and then wrap it with data.
                                 // is not invoked when MLVec goes away
    void    relinquishData();    // disassociates data from MLVec, so free 
                                 // is not invoked when MLVec goes away
    void    resize(int newSize); // change vector size. If newSize is smaller
                                 // than size(), keep first newSize values
    bool    empty() const;       // indicates whether vector is empty

    inline T& operator [] (int index) { return data_[index]; }
    inline const T& operator [] (int index) const { return data_[index]; }

  private:
    int    size_;
    T      *data_;
    bool   okToChangePtr_;   // These last few indicate the status of vector.
    free_s freeType_;        // This is used to see if memory can be freed, 
    bool   getPtrInvoked_;   // reallocated, or pointers can be moved. 
};

// Create an empty vector 

template <class T> MLVec<T>::MLVec() {
  data_          = NULL;      size_          = 0;
  okToChangePtr_ = true;      getPtrInvoked_ = false;
  freeType_      = neither;
}

// Create an vector without any data but allocate space
template <class T> MLVec<T>::MLVec(int mysize) {

  TEUCHOS_TEST_FOR_EXCEPTION(mysize < 0,std::logic_error,
  "MLVec error, cannot create a negative size (= " << mysize << ") vector\n");

  data_          = NULL;        size_          = mysize;
  okToChangePtr_ = true;        getPtrInvoked_ = false;
  freeType_      = neither;

  if (size_ > 0) {
     data_ = (T *) ML_allocate(sizeof(T)*size_); 

     TEUCHOS_TEST_FOR_EXCEPTION(data_ == NULL,std::logic_error,
     "MLVec error, not enough space for vector of size " << size_ << "\n");

     freeType_ = useFree;
  }
}

// Create a vector that wraps existing data that starts and ends at 'start'
// and 'end' respectively. In most cases, this means that the pointers
// associated with the raw data cannot be changed, and that the raw data
// should remain when this vector goes away 
template <class T> MLVec<T>::MLVec(T *start, T *end, bool okToChangePtr,
                                   free_s freeType) { 

  TEUCHOS_TEST_FOR_EXCEPTION(end < start,std::logic_error,
  "MLVec error, cannot create a negative size (= " << end-start<< ") vector\n");

  freeType_      = freeType;   okToChangePtr_ = okToChangePtr;
  data_          = start;      size_          = end - start;
  if (size_ == 0) data_ = NULL;
  getPtrInvoked_ = false;
}
template <class T> void MLVec<T>::wrap(T *start, T *end, bool okToChangePtr,
                                   free_s freeType) { 

  TEUCHOS_TEST_FOR_EXCEPTION((data_ != NULL)||(size_  != 0)||(okToChangePtr_ != true)||
                             (getPtrInvoked_ != false)||(freeType_!= neither),
                             std::logic_error,
  "MLVec wrap error, wrap() must be invoked on vectors created with default constructor\n");

  TEUCHOS_TEST_FOR_EXCEPTION(end < start,std::logic_error,
  "MLVec error, cannot create a negative size (= " << end-start<< ") vector\n");

  freeType_      = freeType;   okToChangePtr_ = okToChangePtr;
  data_          = start;      size_          = end - start;
  if (size_ == 0) data_ = NULL;
  getPtrInvoked_ = false;
}

template <class T> MLVec<T>::~MLVec()              { 
   if (data_ != NULL) {
      TEUCHOS_TEST_FOR_TERMINATION(size_ == 0,
      "MLVec error, data pointer should be null for 0 length vectors\n");
      if (okToChangePtr_ && (freeType_ == useFree))   ML_free(data_); 
      if (okToChangePtr_ && (freeType_ == useDelete)) delete data_;
   }
}

// data_ is no longer associated with this vector. The raw data continues to
// exist.  Can be used in conjunction with getptr() to first get the raw
// data out of MLVec and then effectively empty MLVec.
template <class T> void MLVec<T>::relinquishData() { 

  TEUCHOS_TEST_FOR_EXCEPTION((data_ != NULL) && (getPtrInvoked_ == false),
  std::logic_error,"relinquishData without invoking getPtr() should lead"
  " to memory leak\n");

  data_          = NULL;     size_          = 0; 
  okToChangePtr_ = true;     getPtrInvoked_ = false;
  freeType_      = neither;
}
template <class T> int MLVec<T>::size() const  {   return size_; }
template <class T> T* MLVec<T>::getptr() {getPtrInvoked_ =true; return data_;}
template <class T> bool MLVec<T>::empty() const {if (size_ > 0) return false; 
                                                 else return true;}

// change the size of the vector. There are lots of restrictions as to when 
// this is allowed. There are also lots of cases depending on whether the 
// original vector was zero length (i.e. just allocate memory), the new vector
// is zero length (i.e. just free memory), the new vector is longer or shorter
// than the original vector.
template <class T> void MLVec<T>::resize(int newsize) {

  if (newsize == size_) return;

  TEUCHOS_TEST_FOR_EXCEPTION( newsize < 0, 
  std::logic_error,"MLVec error, cannot resize to a negative length (= " << 
  newsize << ") vector \n");

  TEUCHOS_TEST_FOR_EXCEPTION((newsize > size_) && (size_ > 0) && 
                             (freeType_ == useDelete), std::logic_error,
  "MLVec error, cannot resize vector allocated with new\n");

  if ( (newsize <= size_) && (size_ > 0) && (freeType_ == useDelete) &&
       (okToChangePtr_) ) { size_ = newsize; return; }

  if (okToChangePtr_) {
     if (size_ == 0) {
       if ( newsize > 0) {
         data_ = (T *) ML_allocate(sizeof(T)*newsize);

         TEUCHOS_TEST_FOR_EXCEPTION( data_ == NULL, std::logic_error,
         "MLVec error, not enough resize space for " << newsize <<
         " length vector\n");

         freeType_ = useFree;
       }
       else data_ = NULL;
     }
     else {
       if ( newsize > 0) {
          T* newStuff = (T *) ML_allocate(sizeof(T)*newsize);
          if (size_ < newsize) 
             for (int i = 0; i < size_; i++) newStuff[i] = data_[i];
          else 
             for (int i = 0; i < newsize; i++) newStuff[i] = data_[i];
          ML_free(data_);
          data_ = newStuff;
         TEUCHOS_TEST_FOR_EXCEPTION( data_ == NULL, std::logic_error,
         "MLVec error, not enough resize realloc space for " << newsize <<
         " length vector\n");

       }
       else {
          if (freeType_ == useFree) { ML_free(data_); data_ = NULL; }
          else if (freeType_ == useDelete) { delete data_; data_ = NULL; }

         TEUCHOS_TEST_FOR_EXCEPTION( freeType_ == neither, std::logic_error,
         "MLVec error, do not know how to free this vector\n"); 
        
       }
     }
     size_ = newsize;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(okToChangePtr_ == false, std::logic_error,
  "MLVec error, not allowed to resize this type of wrapped vector\n"); 
}


/******************************************************************************
  Data structures & functions to hide communication details associated with 
  epetra matrices, ML matrices, or MueLu matrices (i.e.  xpetra operators). 
******************************************************************************/

enum epetraOrMLOrMueLu {epetraType , mlType, mueluType}; 

struct wrappedCommStruct {
    epetraOrMLOrMueLu whichOne;
    void                 *data;  // matrix on which communication is based
    int                 nProcs;  // total number of processors participating
    int                  myPid;  
    int          maxDofPerNode;
    int                vecSize;  
};

template <class T> int nodalComm(MLVec<T>& vector, MLVec<int>& myLocalNodeIds,
                                 struct wrappedCommStruct& framework)
{
   /***************************************************************************
    Performs communication to update an amalgamated/nodal vector using either
    an ML style or an epetra style unamalgamated/dof matrix to define the
    communication/import pattern. To do this, it utilzes myLocalNodeIds[i]
    which indicates that the ith dof lies within the myLocalNodeIds[i]th
    node. The ghost portion of myLocalNodeIds[] was computed by 
    assignGhostLocalNodeIds().

    Note: As several dofs lie within each node, several of myLocalNodeIds[]'s
    entries will be duplicates. 
   ***************************************************************************/

    MLVec<T> temp(myLocalNodeIds.size());

    // copy from nodal vector to dof vector
    for (int i = 0; i < myLocalNodeIds.size(); i++)
       temp[i] = vector[ myLocalNodeIds[i]];

    dofComm(temp, framework);

    // copy from dof vector to nodal vector
    for (int i = 0; i < myLocalNodeIds.size(); i++)
       vector[ myLocalNodeIds[i]] = temp[i];

    return 0;
}

template <class T> int dofComm(MLVec<T>& vector,
                               struct wrappedCommStruct& framework)
{
   /***************************************************************************
    Performs communication to update vector using either an ML style matrix
    or an epetra style matrix to define the communication/import pattern.

    Note: Directly utilizes the communication of the underlying matrix
    as opposed to nodalComm() and dofCommUsingMlNodalMatrix(). nodalComm() uses
    an unamalgamated matrix to define the communication associated with an
    amalgamated matrix. dofCommUsingMlNodalMatrix() uses an amalgamated/nodal
    matrix maxDofPerNode times to interpolate each dof of an unamalgamated
    matrix. dofCommUsingMlNodalMatrix() assumes all nodes have maxDofPerNode
    degrees-of-freedom (so it is only appropriate on coarser levels).
   ***************************************************************************/

   MLVec<double> aCopy(vector.size());

   // need to make a copy of the data as we will only use framework comm()
   // functions that work with doubles

   for (int i = 0; i < vector.size(); i++) aCopy[i] = (double) vector[i];

   double *vectorData = aCopy.getptr();

   if (framework.whichOne == mlType) {

      ML_Operator *Amat = (ML_Operator *) framework.data;

      // grab the data pointer

      if (Amat == NULL) return 0;
      if (Amat->getrow == NULL) return 0;
      if (Amat->getrow->pre_comm == NULL) return 0;

      ML_exchange_bdry(vectorData,Amat->getrow->pre_comm,
                 Amat->invec_leng,Amat->comm,ML_OVERWRITE,NULL);
   }
   else {
      Epetra_CrsMatrix *Amat = (Epetra_CrsMatrix *) framework.data;
      Epetra_RowMatrix *ArMat = (Epetra_RowMatrix *) Amat;

      ML_Epetra_comm_wrapper(vectorData, (void *) ArMat);
   }

   for (int i = 0; i < vector.size(); i++) vector[i] = (T) aCopy[i];

   return(0);
}

extern "C" int dofCommUsingMlNodalMatrix(double *data, void *widget);

extern int MLsortCols(const MLVec<int>& ARowPtr, MLVec<int>& ACols, 
                    MLVec<double>& AVals);

extern int MLnMyGhost(struct wrappedCommStruct& framework);

extern int MLfillNodalMaps(MLVec<int> &amalgRowMap, MLVec<int> &amalgColMap,
        MLVec<int> &myLocalNodeIds, int nLocalDofs, 
        struct wrappedCommStruct &framework,
        int nLocalNodes, int nLocalPlusGhostNodes);

extern int MLcolGlobalIds(struct wrappedCommStruct& framework, MLVec<int>& myGids);

extern int MLassignGhostLocalNodeIds(MLVec<int> &myLocalNodeIds, int nLocalDofs,
                   int nLocalPlusGhostDofs, struct wrappedCommStruct &framework,
                   int &nLocalNodes, int &nLocalPlusGhostNodes);

extern int MLextractDiag(const MLVec<int>& rowPtr, const MLVec<int>& cols,
                       const MLVec<double>& , MLVec<double>& diagonal,
                       struct wrappedCommStruct &framework);

extern int MLfindDirichlets(const MLVec<int>& rowPtr, const MLVec<int>& cols,
                          const MLVec<double>& vals,
                          const MLVec<double>& diagonal,
                          double tol, MLVec<bool>& dirOrNot,
                          struct wrappedCommStruct &framework);

extern int MLrmDirichletCols(MLVec<int>& rowPtr, MLVec<int>& cols,
                           MLVec<double>& vals, const MLVec<double>& diagonal,
                           bool squeeze, MLVec<double>& solution,
                           MLVec<double>& rhs, const MLVec<bool>& dirOrNot, 
                           struct wrappedCommStruct &framework);

extern int MLsqueezeOutNnzs(MLVec<int>& rowPtr, MLVec<int>& cols,
                          MLVec<double>& vals, const MLVec<bool>& keep);

extern int MLbuildMap(const MLVec<bool>& dofPresent, MLVec<int>& map, int NDof);

extern int MLvariableDofAmalg(int nCols, const MLVec<int>& rowPtr,
                          const MLVec<int>& cols, const MLVec<double>& vals,
                          int nNodes, int maxDofPerNode, const MLVec<int>& map,
                          const MLVec<double>& diag, double tol,
                          MLVec<int>& amalgRowPtr, MLVec<int>& amalgCols,
                          struct wrappedCommStruct &framework,
                          MLVec<int>& myLocalNodeIds);


extern int MLrmDifferentDofsCrossings(const MLVec<bool>& dofPresent,
                                    int maxDofPerNode, MLVec<int>& rowPtr,
                                    MLVec<int>& cols, int nCols,
                                    struct wrappedCommStruct& framework, 
                                    MLVec<int> &myLocalNodeIds);

extern int MLbuildLaplacian(const MLVec<int>& rowPtr, const MLVec<int>& cols,
                          MLVec<double>& vals, const MLVec<double>& x,
                          const MLVec<double>& y, const MLVec<double>& z);

extern int MLunamalgP(const MLVec<int>& amalgRowPtr,
                    const MLVec<int>& amalgCols,
                    const MLVec<double>& amalgVals, int maxDofPerNode,
                    const MLVec<char>& status, bool fineIsPadded,
                  MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals);

extern int MLfindEmptyRows(const MLVec<int>& rowPtr, const MLVec<int>& cols,
                         const MLVec<double>& vals, MLVec<bool>& rowEmpty);

extern int MLreplaceEmptyByDirichlet(MLVec<int>& rowPtr, MLVec<int>& cols,
                        MLVec<double>& vals, const MLVec<bool>& colEmpty);

extern int MLfineStatus(const MLVec<bool>& dofPresent,
                           const MLVec<int>& map, const MLVec<bool>& dirOrNot,
                           MLVec<char>& status);

extern int MLcoarseStatus(const MLVec<bool>& rowEmpty,
                        const MLVec<bool>& dirOrNot, MLVec<char>& status);



extern int MLShove(ML_Operator *Mat, MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals, int invec_leng, int (*commfunc  )(double *vec, void *data), struct wrappedCommStruct& framework, int nGhost);


extern int ZeroDist(MLVec<double>& xxx, MLVec<double>& yyy, MLVec<double>& zzz,
             MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals,
             MLVec<double>& diagonal, double tol, MLVec<int>& rowZeros, MLVec<int>& colZeros,
             double disttol);

extern int MergeShared(MLVec<int>& cols, MLVec<int>& rowZeros, MLVec<int>& colZeros, 
                       MLVec<int>& groupHead, MLVec<int>& groupNext);

extern int BuildNonSharedMap(MLVec<int>& newMap, MLVec<int>& groupHead, MLVec<int>& groupNext);

extern int buildCompressedA(MLVec<int>& inputRowPtr, MLVec<int>& inputCols,
                     MLVec<double>& inputVals, MLVec<double>& diagonal,
                     MLVec<int>& groupHead, MLVec<int>& groupNext, 
                     MLVec<int>& outputRowptr, MLVec<int>& outputCols,
                     MLVec<int>& map, int newN, double tol);


#endif /* NewStuff */
#endif /* defined HAVE_ML_EPETRA and HAVE_ML_TEUCHOS */

#endif /* define ML_MULTILEVELPRECONDITIONER_H */
