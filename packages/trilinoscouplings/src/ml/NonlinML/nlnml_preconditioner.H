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
 * \file nlnml_preconditioner.H
 *
 * \class NLNML_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 6-April-06
 *
 */
#include "ml_common.h"
#include "ml_include.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_PRECONDITIONER_H
#define NLNML_PRECONDITIONER_H

// ML-headers
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"

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
#include "Epetra_IntVector.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

// Teuchos
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"


// NOX-headers
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Group.H"
//JJH
//#include "NOX_Solver_Manager.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_Vector.H"

// own headers
#include "nlnml_preconditioner_utils.H"
#include "nlnml_finelevelnoxinterface.H"
#include "nlnml_nonlinearlevel.H"

using namespace Teuchos;
//! NLNML: default namespace for all ML/NOX nonlinear methods
namespace NLNML
{

// ---------- Forward Declarations ----------
class NLNML_FineLevelNoxInterface;

/*!

   \brief NLNML_Preconditioner: a class to define a multilevel nonlinear preconditioner and solver

   Class NLNML::NLNML_Preconditioner defines an algebraic nonlinear preconditioner
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

    \note
    an example is given in the directory ml/examples/NonlinML

    \author Michael Gee, SNL 9214
*/
class NLNML_Preconditioner : public virtual Epetra_Operator,
                             public virtual NOX::Epetra::Interface::Preconditioner
{
public:

  //! Constructs an NLNML_Preconditioner class with given parameters
  /*!
   Constructs a NLNML_Preconditioner class and sets parameters from a parameter list
   This constructor does NOT actually create the Multilevel hierarchy and the
   methods, this is done by a call to computePreconditioner.

   The complete list of parameters is:
   \code
   Teuchos::ParameterList mlparams;
   mlparams.set("nlnML output",                                      6          ); // ML-output-level (0-10)
   mlparams.set("nlnML max levels",                                  10         ); // max. # levels (minimum = 2 !)
   mlparams.set("nlnML coarse: max size",                            5000       ); // the size ML stops generating coarser levels
   mlparams.set("nlnML is linear preconditioner",                    false      ); // linear or nonlinear preconditioner
   mlparams.set("nlnML apply constraints",                           true       ); // make application apply constraints
   mlparams.set("nlnML is matrixfree",                               true       ); // create Jacobian yourself or ask application for it
   mlparams.set("nlnML finite difference fine level",                true       ); // create jacobian on fine grid only
   mlparams.set("nlnML finite difference alpha",                     1.0e-08    ); // finite differencing perturbation value
   mlparams.set("nlnML finite difference beta",                      1.0e-07    ); // finite differencing perturbation value
   mlparams.set("nlnML finite difference centered",                  false      ); // use centered or forward finite differencing
   mlparams.set("nlnML Jacobian fix diagonal",                       true       ); // check Jacobian for zero diagonal entries and fix them

   mlparams.set("nlnML absolute residual tolerance",                 1.0e-06    ); // internal tolerance for nonlinear preconditioner
   mlparams.set("nlnML max cycles",                                  500        ); // max # cycles before we give up
   mlparams.set("nlnML adaptive recompute",                          0.0        ); // recompute if residual is larger then this value
   mlparams.set("nlnML offset recompute",                            0          ); // every offset this preconditioner is recomputed
   mlparams.set("nlnML additional adaptive nullspace",               0          ); // compute adaptive nullspace (additional kernel vectors)
   mlparams.set("nlnML PDE equations",                               1          ); // dof per node
   mlparams.set("nlnML null space: dimension",                       1          ); // dimension of nullspace
   mlparams.set("nlnML spatial dimension",                           1          ); // spatial dimension of problem: 1, 2, 3
   mlparams.set("nlnML coarse: type",                                "Uncoupled"); // Uncoupled METIS VBMETIS
   mlparams.set("nlnML nodes per aggregate",                         3          ); // # nodes per agg for coarsening METIS and VBMETIS

   mlparams.set("nlnML use nlncg on fine level",                     true       ); // use nlnCG or mod. Newton's method
   mlparams.set("nlnML use nlncg on medium level",                   true       );
   mlparams.set("nlnML use nlncg on coarsest level",                 false      );

   mlparams.set("nlnML max iterations newton-krylov fine level",     0          ); // # iterations of lin. CG in mod. Newton's method
   mlparams.set("nlnML max iterations newton-krylov medium level" ,  0          );
   mlparams.set("nlnML max iterations newton-krylov coarsest level", 5          );

   mlparams.set("nlnML linear smoother type fine level",             "MLS"      ); // SGS BSGS Jacobi MLS Bcheby AmesosKLU
   mlparams.set("nlnML linear smoother type medium level",           "SGS"      );
   mlparams.set("nlnML linear smoother type coarsest level",         "AmesosKLU");
   mlparams.set("nlnML linear smoother sweeps fine level",           -3         ); // number of sweeps or polynomial order
   mlparams.set("nlnML linear smoother sweeps medium level",         3          );
   mlparams.set("nlnML linear smoother sweeps coarsest level",       1          );

   mlparams.set("nlnML nonlinear presmoothing sweeps fine level",    1          ); // # nonlinear presmoothing sweeps in FAS-cycle fine level
   mlparams.set("nlnML nonlinear presmoothing sweeps medium level",  0          );
   mlparams.set("nlnML nonlinear smoothing sweeps coarse level",     3          );
   mlparams.set("nlnML nonlinear postsmoothing sweeps medium level", 2          );
   mlparams.set("nlnML nonlinear postsmoothing sweeps fine level",   3          );
   \endcode

   \param interface (In) : rcptr to an NLNML::NLNML_FineLevelNoxInterface derived class
   \param mlparams  (In) : paramter list
   \param comm      (In) : the Epetra_Comm to be used
   */
  explicit NLNML_Preconditioner(
                     RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> interface,
                     ParameterList& mlparams,
                     const Epetra_Comm& comm);

  //! Destroys the preconditioner.
  /*!
   Destroys the NLNML_Preconditioner.
   */
  virtual ~NLNML_Preconditioner();

  //! Register the outer Nox iterative method with the nonlinear preconditioner.
  /*!
   Register the outer Nox iterative method with the nonlinear preconditioner.
   This is necessary, because the outer NOX iterative method will provide the
   nonlinear preconditioner with current solution and residual. Necessary only
   in the nonlinear preconditioner case and when used as a preconditioner. If used as
   stand-alone solver, this is not necessary.
   \param s   (In) : rcptr to the NOX::Solver::Generic of the outer
                     Nox iterative method.
  */
  //JJH
  //void SetNoxSolver(RefCountPtr<NOX::Solver::Manager> s) { noxsolver_ = s; }
  void SetNoxSolver(RefCountPtr<NOX::Solver::Generic> s) { noxsolver_ = s; }

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
  const Epetra_Map & OperatorDomainMap() const {return interface_->getMap();}

  //! Get fine level OperatorRangeMap.
  /*!
   Derived from Epetra_Operator, get fine level OperatorRangeMap.
   */
  const Epetra_Map & OperatorRangeMap() const {return interface_->getMap();}

  //! Returns the level of output to be generated to stdout ( 0 - 10 )
  inline int OutLevel() const { return params_->get("nlnML output",0); }

  //! Check initialization flag
  /*!
   check the init-flag, if false, ApplyInverse will call computePreconditioner
   */
  inline bool isinit() const {return(isinit_);}

  //! Set the initialization flag
  bool setinit(bool flag)
  {isinit_=flag; ncalls_NewPrec_ = 0; return(true);}

  //! Get the offset to recompute this preconditioner.
  int getoffset() const {return params_->get("nlnML offset recompute",0);}

  // derived from NOX::Epetra::Interface::Preconditioner
  //! Construct the preconditioner.
  /*!
   Construct the preconditioner. This is derived from the
   NOX::Epetra::Interface::Preconditioner class. Nox will use this to
   (re)compute the preconditioner.
   \param x          (In) : Vector of the current solution
   \param M          (In) : this
   \param precParams (In) : ignored
   */
  virtual bool computePreconditioner(const Epetra_Vector& x,
             Epetra_Operator& M,
             Teuchos::ParameterList* precParams = 0);

  //! ApplyInverse the preconditioner
  /*!
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator.
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
   */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  //! Apply the nonlinear multigrid FAS-solver.
  /*!
   Apply the nonlinear multigrid FAS-solver.
   This preconditioner can act as a stand-alone FAS-Solver. Will throw an error,
   if the preconditioner is set up as a linear preconditioner
   */
  int solve() const ;

  //! Apply the preconditioner, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  { std::cout << "ML_Nox_Preconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";return(-1);}

  //! SetUseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  int SetUseTranspose(bool UseTranspose) { return -1; }

  //! NormInf, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  double NormInf() const {return(-1.0);}

  //! UseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  bool UseTranspose() const {return false;}

  //! HasNormInf, not implemented.
  /*!
   Derived from Epetra_Operator, always returns false.
   */
  bool HasNormInf() const {return(false);}

private:

  explicit NLNML_Preconditioner(NLNML::NLNML_Preconditioner& old);
  NLNML_Preconditioner operator = (const NLNML_Preconditioner& old);

  void Initialize();

  bool compPrec(const Epetra_Vector& x);
  bool Compute_Linearpreconditioner(const Epetra_Vector& x);
  bool Compute_Nonlinearpreconditioner(const Epetra_Vector& x);

  void Set_Smoothers();
  void fix_MainDiagonal(RefCountPtr<Epetra_CrsMatrix> A, int level);
  RefCountPtr<Epetra_CrsMatrix> ComputeFineLevelJacobian(const Epetra_Vector& x);

  int ApplyInverse_Linear(const Epetra_MultiVector& X,
                          Epetra_MultiVector& Y) const;
  int ApplyInverse_NonLinear(const Epetra_MultiVector& X,
                             Epetra_MultiVector& Y) const;

  //the FAS-preconditioner
  bool FAS_cycle(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged, double* finenorm) const;
  // the FAS-solver
  bool FAS_cycle_solver() const;

  bool Adaptivesetup(double** oldns,Epetra_CrsMatrix* Jac,int oldnumpde,int olddimns);

  template<typename kind>
  kind getParameter(string name,kind dvalue) const
  { return params_->get(name,dvalue); }

  template<typename kind>
  void setParameter(string name,kind dvalue) const
  { params_->set(name,dvalue); return; }

  bool CheckInputParameters(ParameterList& params);

private:

  const char*                                                     label_;
  bool                                                            isinit_;
  int                                                             ncalls_NewPrec_;

  const Epetra_Comm&                                              comm_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface>                 interface_;
  RefCountPtr<Teuchos::ParameterList>                             params_;
  //JJH
  //RefCountPtr<NOX::Solver::Manager>                               noxsolver_;
  RCP<NOX::Solver::Generic>                                       noxsolver_;
  RefCountPtr<Epetra_CrsMatrix>                                   fineJac_;
  RefCountPtr<Epetra_Map>                                         rowmap_;
  RefCountPtr<ML_Epetra::MultiLevelOperator>                      linPrec_;
  RefCountPtr<std::vector<RefCountPtr<NLNML::NLNML_NonlinearLevel> > > nlnlevel_;
  ML*                                                             ml_;
  ML_Aggregate*                                                   ag_;


}; // class NLNML_Preconditioner
} // namespace NLNML

#endif
#endif
