// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file ml_nox_preconditioner.H
 *
 * \class ML_Nox_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_PRECONDITIONER_H
#define ML_NOX_PRECONDITIONER_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

// NOX-headers
#include "NOX_Epetra_Interface.H"

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"
// a matrixfree level
#include "ml_nox_matrixfreelevel.H"
// a nonlinear level
#include "ml_nox_nonlinearlevel.H"
// utils
#include "ml_nox_preconditioner_utils.H"
// this is the wrapper for the Epetra_CrsGraph used by ml in matrix-free method
// it is derived from Epetra_RowMatrix so the ml - Epetra_RowMatrix - wrapper will take it
#include "ml_Epetra_wrap_CrsGraph_as_RowMatrix.H"

// interface to user code on fine level
#include "ml_nox_fineinterface.H"



// ---------- Forward Declarations ----------
class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;
class Nox_Problem_Interface;
//! ML_NOX: default namespace for all ML/NOX nonlinear methods
namespace ML_NOX
{
/*!
 
   \brief ML_Nox_Preconditioner: a class to define a multilevel nonlinear preconditioner and solver

   Class ML_NOX::ML_Nox_Preconditioner defines an algebraic nonlinear preconditioner
   or solver that is based on a full approximation scheme (FAS) approach and
   is implemented making use of NOX' nonlinear solver capabilities.
   It operates on Epetra objects and implements an Epetra_Operator to be easily
   used as a preconditioner to some outer nonlinear iteration.
   It also implements the NOX::EpetraNew::Interface::Preconditioner from NOX
   so it can serve as it's own interface (see NOX documentation for details).
   In this class several methods and options can be chosen, which are described in the
   parameter set functions
   
  
    To make use of the full functionality of this class it 
    requires ML to be configured with the following options:
    - \c --with-ml_nox
    - \c --enable-epetra
    - \c --enable-epetraext
    - \c --enable-nox
    - \c --enable-nox-epetra
    - \c --enable-prerelease
    - \c --enable-aztecoo
    - \c --enable-amesos
    
   \warning
    Currently, using this class makes it necessary to configure and make install
    of Trilinos WITHOUT the --with-ml_nox, and then repeating the configure
    and make install process WITH the --with-ml_nox option to be able to resolve the
    dependencies of this class on NOX header files
    
    \note
    an example is given in the directory ml/examples/NonlinML

    \author Michael Gee, SNL 9214
*/  
class  ML_Nox_Preconditioner : public virtual Epetra_Operator, 
                               public NOX::EpetraNew::Interface::Preconditioner
{
public:

  

  //! Constructs an ML_Nox_Preconditioner class with default values.
  /*! 
   Constructs a ML_Nox_Preconditioner class and sets a set of default values.
   This constructor does NOT actually create the Multilevel hierarchy and the
   methods, this is done by a call to computePreconditioner.
   \param interface (In) : ref to an ML_NOX::Ml_Nox_Fineinterface derived class
   \param comm      (In) : the Epetra_Comm to be used
   */
  explicit ML_Nox_Preconditioner(ML_NOX::Ml_Nox_Fineinterface& interface,
                                 const Epetra_Comm&            comm);

  //! Constructs an ML_Nox_Preconditioner class with parameter list.
  /*! 
   Constructs a ML_Nox_Preconditioner class and sets a set of default values.
   This constructor does NOT actually create the Multilevel hierarchy and the
   methods, this is done by a call to computePreconditioner.<br>
   
   <b>Teuchos parameters recognized and there default values:</b>
   \code
   Teuchos::ParameterList mlparams;
   mlparams.set("nlnML output",6);                                        // Level of output to be generated (0-10), 6 recommended
   mlparams.set("nlnML max levels",10);                                   // Max # of levels to produce
   mlparams.set("nlnML coarse: max size",128);                            // Size of coarse grid where to stop coarsening
   mlparams.set("nlnML is linear preconditioner",false);                  // Act as linear preconditioner, else run FAS-cycle as nonlinear preconditioner
   mlparams.set("nlnML is matrixfree",true);                              // Finite difference for Jacobian
   mlparams.set("nlnML finite difference fine level",true);               // Finite difference for Jacobian on fine level only
   
   mlparams.set("nlnML finite difference alpha",1.0e-07);                 // Finite differencing perturbation parameter
   mlparams.set("nlnML finite difference beta",1.0e-06);                  // ---"---
   mlparams.set("nlnML finite difference centered",false);                // Finite differnencing scheme, true=centered, false=forward

   mlparams.set("nlnML absolute residual tolerance",1.0e-05);             // Internal absolute residual tolerance for FAS-cycle
   mlparams.set("nlnML max cycles",250);                                  // Max. # FAS-cycles before giving up
   mlparams.set("nlnML adaptive recompute",0.0);                          // Recompute preconditioner when outer residual is above this tolerance
   mlparams.set("nlnML offset recompute",0);                              // Recompute preconditioner every offset outer iterations
   mlparams.set("nlnML additional adaptive nullspace",0);                 // # additional adaptive nullspace candidates to compute
   
   mlparams.set("nlnML PDE equations",3);                                 // # PDE equations (degrees of freedom per node)
   mlparams.set("nlnML null space: dimension",3);                         // nullspace dimension
   mlparams.set("nlnML spatial dimension",3);                             // spatial dimension of problem
   
   mlparams.set("nlnML coarse: type","Uncoupled");                        // aggregation scheme: Uncoupled METIS VBMETIS
   mlparams.set("nlnML nodes per aggregate",27);                          // # nodes per aggregate for aggregation schemes "METIS" / "VBMETIS"

   mlparams.set("nlnML use nlncg on fine level",true);                    // use nonlinear CG as smoother on fine level, else use Newton-Krylov
   mlparams.set("nlnML use nlncg on medium level",true);                  // use nonlinear CG as smoother, else use Newton-Krylov
   mlparams.set("nlnML use nlncg on coarsest level",true);                // use nonlinear CG as smoother on coarsest level, else use Newton-Krylov
   
   mlparams.set("nlnML max iterations newton-krylov fine level",5);       // # iterations LINEAR CG in Newton-Krylov smoother fine level
   mlparams.set("nlnML max iterations newton-krylov medium level" ,5);    // # iterations LINEAR CG in Newton-Krylov smoother
   mlparams.set("nlnML max iterations newton-krylov coarsest level",5);   // # iterations LINEAR CG in Newton-Krylov smoother coarsest level 

   mlparams.set("nlnML linear smoother type fine level","SGS");           // type of linear fine smoother  : SGS BSGS Jacobi MLS Bcheby AmesosKLU
   mlparams.set("nlnML linear smoother type medium level","SGS");         // type of linear smoother       : SGS BSGS Jacobi MLS Bcheby AmesosKLU
   mlparams.set("nlnML linear smoother type coarsest level","AmesosKLU"); // type of linear coarse smoother: SGS BSGS Jacobi MLS Bcheby AmesosKLU 
   mlparams.set("nlnML linear smoother sweeps fine level",2);             // # of linear smoother sweeps/polynomial order for linear fine smoother
   mlparams.set("nlnML linear smoother sweeps medium level",2);           // # of linear smoother sweeps/polynomial order for linear smoother
   mlparams.set("nlnML linear smoother sweeps coarsest level",1);         // # of linear smoother sweeps/polynomial order for linear coarse smoother

   mlparams.set("nlnML nonlinear presmoothing sweeps fine level",0);      // # nonlinear presmoothing sweeps in FAS-cycle fine level
   mlparams.set("nlnML nonlinear presmoothing sweeps medium level",0);    // # nonlinear presmoothing sweeps in FAS-cycle
   mlparams.set("nlnML nonlinear smoothing sweeps coarse level",3);       // # nonlinear sweeps in FAS-cycle coarsest level
   mlparams.set("nlnML nonlinear postsmoothing sweeps medium level",3);   // # nonlinear postsmoothing sweeps in FAS-cycle
   mlparams.set("nlnML nonlinear postsmoothing sweeps fine level",3);     // # nonlinear postsmoothing sweeps in FAS-cycle fine level
   \endcode
   
   \param interface (In) : ref to an ML_NOX::Ml_Nox_Fineinterface derived class
   \param mlparams  (In) : Teuchos parameter list 
   \param comm      (In) : the Epetra_Comm to be used
   */
  explicit ML_Nox_Preconditioner(ML_NOX::Ml_Nox_Fineinterface& interface,
                                 Teuchos::ParameterList& mlparams,
                                 const Epetra_Comm& comm);
  

  //! Destroys the preconditioner.
  /*! 
   Destroys the ML_Nox_Preconditioner.
   */
  ~ML_Nox_Preconditioner();

  

  //! Set the type of nonlinear method to be used.
  /*! 
   Set the type of nonlinear method to be used on each level of the hierarchy
   \param islinPrec (In) :     if true, compute a linear multigrid operator,
                               if false, create a FAS multilevel hierarchy
   \param maxlevel   (In) :    maximum number of levels.
   \param ismatrixfree (In) :  if true, use Finite Differencing to construct the Jacobians
                               on all levels. In this case, the preconditioner will
                               ask for a Epetra_CrsGraph through the fine level interface
                               and will build a plain aggregation hierarchy based on it.
                               If Jacobi is chosen as a smoother on any level, on that level
                               only the main diagonal of the Jacobian will be computed.
                               If false, the preconditioner will ask for an Epetra_CrsMatrix
                               and an Epetra_CrsGraph through the fine level interface and will
                               construct a smoothed aggregation MG hierarchy based on it.
                               Coarse level operators will then be Galerkin.
   \param ismatfreelev0 (In) : if true AND ismatrixfree==true, the preconditioner
                               will construct a finite difference Jacobian on the fine level and
                               then will use it to construct a smoothed aggregation hierarchy. 
                               Coarse level operators will be Galerkin. Option has no effect when
                               ismatrixfree==false.
   */
  bool SetNonlinearMethod(bool islinPrec, int maxlevel, bool ismatrixfree, 
                          bool ismatfreelev0, bool fixdiagonal);

  //! Set type of coarsening algorithm to be used.
  /*! 
   Set type of coarsening algorithm to be used. Currently, this class recognizes
   the ML-coarsening schemes Uncoupled, METIS and VBMETIS.
   \param coarsentype   (In) : string giving the coarsening scheme, 
                               currently "Uncoupled", "METIS" and "VBMETIS" are
                               recognized
   \param maxlevel      (In) : maximum number of levels to generate
   \param maxcoarsesize (In) : coarse grid size to stop coarsening when reached
   \param nnodeperagg   (In) : objective how many nodes to put in one aggregate for
                               coarsening schemes "METIS" and "VBMETIS".
   */
  bool SetCoarsenType(string coarsentype,int maxlevel,int maxcoarsesize,
                      int nnodeperagg);

  //! Set type of smoother.
  /*! 
   Set type of smoother. Currently, this class recognizes
   the ML-smoothers Jacobi, SGS and AmesosKLU.
   \param finesmoothertype (In) : string giving the fine level smoother, 
                               currently "Jacobi", "SGS" and "AmesosKLU" are
                               recognized.
   \param smoothertype     (In) : string giving the smoother on all
                                  intermediate levels,  currently "Jacobi", "SGS" and "AmesosKLU" are
                               recognized.
   \param coarsesolve      (In) : string giving the smoother on the coarsest level,  
                                  currently "Jacobi", "SGS" and "AmesosKLU" are
                                 recognized.
   */
  bool SetSmoothers(string finesmoothertype,string smoothertype,
                    string coarsesolve);

  //! Set type of nonlinear solver for each level.
  /*! 
   Set type of nonlinear solver for each level and number of LINEAR
   CG-iterations allowed inside Newton's method. 
   \param usenlnCG_fine   (In) : If true, use nlnCG on finest level, if
                                 false, use Newton's method.
   \param usenlnCG        (In) : If true, use nlnCG on all intermediate levels, 
                                 if false use Newton's method.
   \param usenlnCG_coarse (In) : If true, use nlnCG on coarsest level, if
                                 false, use Newton's method.
   \param useBroyden      (In) : If true, use Broyden update for Jacobian, 
                                 will only affect levels on which Newton's
                                 method is used.
   \param nitersCG_fine   (In) : If Newton's method is used on finest level, 
                                 this is the max. number of LINEAR CG iterations
                                 allowed for the linear solve inside.
   \param nitersCG        (In) : If Newton's method is used on intermediate levels, 
                                 this is the max. number of LINEAR CG iterations
                                 allowed for the linear solve inside.
   \param nitersCG_fine   (In) : If Newton's method is used on coarsest level, 
                                 this is the max. number of LINEAR CG iterations
                                 allowed for the linear solve inside.
    \note
    When using Newton's method with or without Broyden update, Aztec's PCG will
    be used as linear solver inside. The preconditioner to this LINEAR CG will then
    be the smoother chosen via the method SetSmoothers.
   */
  bool SetNonlinearSolvers(bool usenlnCG_fine, bool usenlnCG, bool usenlnCG_coarse,
                           bool useBroyden, int nitersCG_fine, int nitersCG, int nitersCG_coarse);


  //! Set number of sweeps of smoother.
  /*! 
   Set number of sweeps of smoother.
   \param nsmooth_fine   (In) : Number of smoothing sweeps on fine level
   \param nsmooth        (In) : Number of smoothing sweeps on intermediate levels
   \param nsmooth_coarse (In) : Number of smoothing sweeps on coarsest level
   \note
   When using nlnCG on a level, this is the number of sweeps of the
   preconditioner to nlnCG. When using Newton's method on a level, this is
   the number of sweeps of the preconditioner to the LINEAR CG inside 
   the Newton iteration on that level.
   */
  bool SetSmootherSweeps(int nsmooth_fine, int nsmooth, int nsmooth_coarse);


  //! Set the ML print level.
  /*! 
   Set the ML print level between 0 and 10.
   \param printlevel   (In) : (0-10), amount of output generated by NOX/ML
   */
  bool SetPrintLevel(int printlevel); 


  //! Set the offset for recomputing this preconditioner.
  /*! 
   Set the offset for recomputing this preconditioner. Every offset times, the
   preconditioner is applied it will recreate itself and will use the
   fine level interface to obtain either a new fine level Jacobian or
   to compute a Finite Differencing - Jacobian.
   \param offset (In) : How often to apply this preconditioner before recomputing it.
   */
  bool SetRecomputeOffset(int offset);

  //! Set the offset for recomputing this preconditioner.
  /*! 
   Set the offset for recomputing this preconditioner. Every offset times, the
   preconditioner is applied it will recreate itself and will use the
   fine level interface to obtain either a new fine level Jacobian or
   to compute a Finite Differencing - Jacobian.
   Optionally, the MG hierarchy and Jacobian can be recomputed after the
   recomputestep iteration
   \param offset (In) : How often to apply this preconditioner before recomputing it.
   */
  bool SetRecomputeOffset(int offset, int recomputestep, double adaptrecompute, int adaptns);

  //! Set dimensions for the PDE and the Nullspace.
  /*! 
   Set dimensions for the PDE and the Nullspace. 
   \param spatialDimension (In) : spatial dimension of the problem.
   \param numPDE           (In) : Number of PDE equations, normally number of dofs per node.
   \param dimNS           (In)  : Dimension of the nullspace. The preconditioner will ask for 
                                  the Nullspace through the fine level interface. If know
                                  Nullspace is supplied there, the default Nullspace is used 
                                  and a warning is printed.
   */
  bool SetDimensions(int spatialDimension, int numPDE, int dimNS);

  //! Set convergence criteria for nonlinear preconditioner.
  /*! 
   Set convergence criteria for nonlinear preconditioner. Takes effect for
   islinPrec==false only.
   \param FAS_normF   (In) : Residual norm a  nonlinear iteration on some level is
                             considered converged.
   \param FAS_nupdate (In)  : Norm of stepsize a nonlinear iteration on some level is
                             considered converged.
   */
  bool SetConvergenceCriteria(double FAS_normF, double FAS_nupdate);

  //! Set cycle for nonlinear preconditioner. 
  /*! 
   Set cycle for nonlinear preconditioner. Takes effect for
   islinPrec==false only.
   \param prefsmooth   (In) : maximum number of nonlinear iterations to take 
                             as a presmoothing step on finest level
   \param presmooth (In)  : maximum number of nonlinear iterations to take 
                             as a presmoothing step on coarse levels
   \param coarsesmooth (In)  : maximum number of nonlinear iterations to take 
                             on coarsest level
   \param postsmooth (In)  : maximum number of nonlinear iterations to take 
                             as a postsmoothing step on coarse levels
   \param postfsmooth (In)  : maximum number of nonlinear iterations to take 
                             as a postsmoothing step on finest level
   \param maxcycle (In)  : maximum number of nonlinear MG cycles to take
                           before we give up
   */

  bool SetFAScycle(int prefsmooth,int presmooth,int coarsesmooth, int postsmooth,
                   int postfsmooth,int maxcycle);

  //! Set finite differencing method.
  /*! 
   Set finite differencing method. Takes effect for
   ismatrixfree==true only.
   \param centered   (In) : if true, use centered finite differencing, else use
                            forward finite differencing (true is twice the cost)
   \param alpha (In)  : perturbation parameter, see EpetraExt manual
   \param beta (In)  : perturbation parameter, see EpetraExt manual
   */
  bool SetFiniteDifferencing(bool centered,double alpha,double beta);

  //! Register the outer Nox iterative method with the nonlinear preconditioner.
  /*! 
   Register the outer Nox iterative method with the nonlinear preconditioner.
   This is necessary, because the outer NOX iterative method will provide the
   nonlinear preconditioner with current solution and residual. Necessary only
   in the islinearPrec==false case and when used as a preconditioner. If used as
   stand-alone solver, this is not necessary.
   \param noxsolver   (In) : ptr to the NOX::Solver::Manager of the outer
                             Nox iterative method.
   */
  bool set_nox_solver(NOX::Solver::Manager* noxsolver) {noxsolver_ = noxsolver; return true;} 
  
  

  //! Returns the label of this class.
  /*! 
   Returns the label of this class.
   */
  const char* Label() const {return(&(label_[0]));}

  //! get Comm of this class
  /*! 
   Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class 
   */
  const Epetra_Comm & Comm() const {return(comm_);}

  //! Get fine level OperatorDomainMap.
  /*! 
   Derived from Epetra_Operator, get fine level OperatorDomainMap.
   */
  const Epetra_Map & OperatorDomainMap() const {return interface_.getMap();}

  //! Get fine level OperatorRangeMap.
  /*! 
   Derived from Epetra_Operator, get fine level OperatorRangeMap.
   */
  const Epetra_Map & OperatorRangeMap() const {return interface_.getMap();}
  
  //! Check initialization flag of the instant of this class
  /*! 
   check the init-flag of this instant, if false, ApplyInverse will throw an error.
   */
  bool isinit() const {return(isinit_);}
  
  //! Set the initialization flag of the instant of this class
  /*! 
   Set the init-flag of this instant.
   */
  bool setinit(bool flag) 
  {isinit_=flag; ncalls_NewPrec_ = 0; return(true);}

  //! Get the offset to recompute this preconditioner.
  /*! 
   Get the offset to recompute this preconditioner.
   */
  int getoffset() {return(offset_newPrec_);} 

  

  // derived from NOX::EpetraNew::Interface::Preconditioner
  //! Construct the preconditioner.
  /*! 
   Construct the preconditioner. This is derived from the 
   NOX::EpetraNew::Interface::Preconditioner class. Nox will use this to
   (re)compute the preconditioner.
   \param x   (In) : Vector of the current solution
   \param precParams (In)  : ignored
   */
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     NOX::Parameter::List* precParams = 0);
  
  //! ApplyInverse the preconditioner
  /*! 
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator. 
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
   */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  // Apply the nonlinear multigrid FAS-solver
  //! Apply the nonlinear multigrid FAS-solver. 
  /*! 
   Apply the nonlinear multigrid FAS-solver. 
   This preconditioner can act as a stand-alone FAS-Solver. Will throw an error,
   if islinearPrec==true; 
   */
  int solve();

  //! Apply the preconditioner, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { cout << "ML_Nox_Preconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";return(-1);}
    
  //! SetUseTranspose, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  int SetUseTranspose(bool UseTranspose);

  //! NormInf, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  double NormInf() const {return(-1.0);}


  //! UseTranspose, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  bool UseTranspose() const {return(usetranspose_);}

  //! HasNormInf, not implemented.
  /*! 
   Derived from Epetra_Operator, always returns false. 
   */
  bool HasNormInf() const {return(false);}

private:
  // don't want copy-ctor and = operator
  ML_Nox_Preconditioner(ML_Nox_Preconditioner& old);
  ML_Nox_Preconditioner operator = (const ML_Nox_Preconditioner& old);
  
  // Init the preconditioner to default values
  bool Initialize();
  
  // Compute Preconditioner functions, 
  // depending on type of preconditioner and whether the Jacobian was supplied or not
  bool compPrec(const Epetra_Vector& x);
  bool ML_Nox_compute_Jacobian_Linearpreconditioner(const Epetra_Vector& xfine);      
  bool ML_Nox_compute_Matrixfree_Linearpreconditioner(const Epetra_Vector& xfine);    
  bool ML_Nox_compute_Jacobian_Nonlinearpreconditioner(const Epetra_Vector& xfine);   
  bool ML_Nox_compute_Matrixfree_Nonlinearpreconditioner(const Epetra_Vector& xfine); 

  // utils
  bool Set_Smoothers();
  bool fix_MainDiagonal(Epetra_CrsMatrix** A, int level) const;
  bool test_nullspace(int dimns, double* nullsp, Epetra_CrsMatrix* Jac) const;
  Epetra_CrsMatrix* ML_Nox_computeFineLevelJacobian(const Epetra_Vector& x);
 
  // ApplyInverse functions, 
  // depending on type of preconditioner
  int ML_Nox_ApplyInverse_Linear(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ML_Nox_ApplyInverse_NonLinear(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  // nonlinear multigrid functions
  //the FAS-preconditioner
  bool ML_Nox_FAS_cycle(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged, double* finenorm) const;
  // the FAS-solver
  bool ML_Nox_FAS_cycle1(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged);

  // destructor of vector of matfreelevels
  bool destroy_matfreelevels(ML_NOX::ML_Nox_MatrixfreeLevel*** m, int nlevel);
  // destructor of vector of nonlinearlevels
  bool destroy_nonlinearlevels(ML_NOX::ML_Nox_NonlinearLevel*** nln, int nlevel);

  // run adatpve setup
  bool Ml_Nox_adaptivesetup(double** oldns, Epetra_CrsMatrix* Jac,int oldnumpde,int olddimns);

private:
  // ---general data---
  char*                                 label_;             // the label of this class
  Teuchos::ParameterList*               params_;           // holding the parameters list
  bool                                  usetranspose_;     // use transpose of operator (not supported)
  bool                                  isinit_;           // flag whether preconditioner is valid
  bool                                  ismatrixfree_;     // use Graph or Jacobian in building the preconditioner
  bool                                  matfreelev0_;
  bool                                  fixdiagonal_;
  bool                                  islinearPrec_;     // do standard linear ML-SA-Preconditioner (not matrixfree!)) 
  ML_NOX::Ml_Nox_Fineinterface&         interface_;        // the user interface
  const Epetra_Comm&                    comm_;             // the communicator;
  bool                                  destroyfineJac_;   // flag to indicate this class destroys the fineJac_
  Epetra_CrsMatrix*                     fineJac_;          // fine grid matrix (only used when ismatrixfree_=false)
  int                                   offset_newPrec_;   // how often to recompute the precond from scratch
  int                                   recompute_newPrec_;// initial recomputing step
  double                                adaptive_NewPrec_; // recompute preconditioner if residual is larger than this
  int                                   adaptns_;          // run adaptive setup to create additional adaptns_ kernel vectors
  int                                   ncalls_NewPrec_;   // how often computePreconditioner is called

  // ---data related to the ML-hierarchy---
  ML*                                   ml_;               // the ml handle
  ML_Aggregate*                         ag_;               // the aggregate handle
  int                                   ml_N_levels_;      // max. number of levels
  int                                   ml_printlevel_;    // ML's printlevel (0 - 10)
  int                                   ml_spatialDimension_; // spatial dimension of the problem
  int                                   ml_numPDE_;        // number of pdes per node
  int                                   ml_dim_nullsp_;    // dimension of nullspace
  int                                   ml_dim_nullsp2_;   // dimension of nullspace
  string                                ml_coarsentype_;   // type of coarsening algorithm: Uncoupled MIS METIS VBMETIS
  int                                   ml_nnodeperagg_;   // objective number of nodes per aggregate with coarsentypes METIS & VBMETIS
  string                                ml_smoothertype_;  // type of smoothing algorithm: SGS BSGS Jacobi
  string                                ml_fsmoothertype_; // type of smoothing algorithm: SGS BSGS Jacobi on level 0
  string                                ml_coarsesolve_;   // type of coarse solver: SGS BSGS AmesosKLU Jacobi
  int                                   ml_maxcoarsesize_; // max. coarse grid size
  int                                   ml_nlevel_;        // number of level generated by ML_Gen_MGHierarchy_UsingAggregation
  int                                   ml_coarsestlev_;   // index of the coarsest level;
  int                                   ml_nblocks_;
  int*                                  ml_blocks_;        // block information for Block-smoothers and aggregation scheme VBMETIS
  int*                                  ml_block_pde_;     // information which pde each row belongs to
  int                                   nsmooth_fine_;     // number of smoothing steps on finest level
  int                                   nsmooth_;          // number of smoothing steps on intermediate levels
  int                                   nsmooth_coarse_;   // number of smoothing steps on coarsest level

  // ---the linear preconditioner---
  ML_Epetra::MultiLevelOperator*        ml_linPrec_;       // the linear preconditioner class
  Epetra_Map*                           ml_map_;           // the linear operators range and domain map map

  // ---matrix-free method--- 
  const Epetra_CrsGraph*                fineGraph_;        // the fine grid graph (only used when ismatrixfree_=true)
  ML_Epetra::CrsGraphWrapper*           ml_graphwrap_;     // in the matrix-free case this wraps the input fine level graph as an Epetra_RowMatrix
  int                                   nmatfreelevel_;    // number of matfreelevels, should be equal to ml_nlevel_
  ML_NOX::ML_Nox_MatrixfreeLevel**      ml_matfreelevel_;  // vector of matrix-free levels (including fine level)
  double                                fd_alpha_;         // FiniteDifferencing parameters (for matrixfree only)
  double                                fd_beta_;
  bool                                  fd_centered_;      // use centered finite differencing

  // --- data special to the nonlinear preconditioner ---
  int                                   n_nlnlevel_;
  ML_NOX::ML_Nox_NonlinearLevel**       nlnLevel_;         // ptr to vector of coarse interfaces
  double                                FAS_normF_;       // coarse NOX's convergence criteria
  double                                FAS_nupdate_;
  int                                   FAS_maxiter_;
  int                                   FAS_prefinesmooth_;
  int                                   FAS_presmooth_;
  int                                   FAS_coarsesmooth_;
  int                                   FAS_postsmooth_;
  int                                   FAS_postfinesmooth_;
  int                                   FAS_maxcycle_; 
  NOX::Solver::Manager*                 noxsolver_;

  bool                                  useBroyden_;       // use Broyden update for Newton's method
  bool                                  usenlnCG_fine_;
  bool                                  usenlnCG_; 
  bool                                  usenlnCG_coarse_;
  int                                   nitersCG_fine_; 
  int                                   nitersCG_; 
  int                                   nitersCG_coarse_;

}; // class  ML_Nox_Preconditioner : public virtual Epetra_Operator
}  // namespace ML_NOX


#endif // ML_NOX_PRECONDITIONER_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
