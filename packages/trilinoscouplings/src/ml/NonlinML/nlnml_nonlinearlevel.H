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
 * \file nlnml_nonlinearlevel.H
 *
 * \class NLNML_NonlinearLevel
 *
 * \brief a nonlinear coarse grid class
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "ml_common.h"
#include "ml_include.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_NONLINEARLEVEL_H
#define NLNML_NONLINEARLEVEL_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_MapColoring.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX-headers
#include <NOX.H>
#include <NOX_Epetra.H>
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"

// own headers
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"
#include "nlnml_prepostoperator.H"
#include "nlnml_ConstrainedMultiLevelOperator.H"

// ---------- Forward Declarations ----------
namespace NLNML {
class NLNML_FineLevelNoxInterface;
}

using namespace Teuchos;
namespace NLNML
{

class NLNML_NonlinearLevel : public virtual NOX::Epetra::Interface::Jacobian,
                             public virtual NOX::Epetra::Interface::Preconditioner
{
public:

  explicit NLNML_NonlinearLevel(
                    int level, RefCountPtr<Teuchos::ParameterList> params,
                    ML* ml, ML_Aggregate* ag,
                    RefCountPtr< std::vector< RefCountPtr<Epetra_CrsMatrix> > > P,
                    RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface,
                    const Epetra_Comm& comm,
                    RefCountPtr<Epetra_CrsMatrix> fineJac,
                    const Epetra_Vector& xfine,
                    bool isnlncg, int niterscg,
                    int numpde, int dimns);

  virtual ~NLNML_NonlinearLevel();

  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  bool computePreconditioner(const Epetra_Vector& x,
           Epetra_Operator& M,
           Teuchos::ParameterList* precParams = 0);

  const Epetra_Comm & Comm() const {return(comm_);}

  inline int OutLevel() const { return params_->get("nlnML output",0); }

  bool setModifiedSystem(bool flag, Epetra_Vector* fbar, Epetra_Vector* fxbar)
  { return coarseinterface_->setModifiedSystem(flag,fbar,fxbar);}

  // evaluate the nonlinear function on this level
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
                NOX::Epetra::Interface::Required::FillType fillFlag)
  { return coarseinterface_->computeF(x,F,fillFlag);}

  bool Iterate(Epetra_Vector* f, Epetra_Vector* x, int numiter, double* norm);

  // restrict vector to next coarser level
  Epetra_Vector* restrict_to_next_coarser_level(Epetra_Vector* finevec, int current, int next)
  { return coarseinterface_->restrict_to_next_coarser_level(finevec,current,next);}

  // prolong vector from next coarser level to current level
  Epetra_Vector* prolong_to_this_level(Epetra_Vector* coarsevec, int current, int next)
  { return coarseinterface_->prolong_to_this_level(coarsevec,current,next);}

private:

  // don't want copy-ctor and = operator
  NLNML_NonlinearLevel(NLNML::NLNML_NonlinearLevel& old);
  NLNML_NonlinearLevel operator = (const NLNML::NLNML_NonlinearLevel& old);

  inline int Level() const { return level_; }

  bool Set_Smoother(ML* ml, ML_Aggregate* ag, int level, int nlevel,
                    ML* thislevel_ml, ML_Aggregate* thislevel_ag,
                    string smoothertype, int nsmooth);

  void create_Nox_Convergencetest(double normf, double norm_update, int maxiter);

  template<typename kind>
  kind getParameter(string name,kind dvalue) const
  { return (*params_).get(name,dvalue); }

  template<typename kind>
  void setParameter(string name,kind dvalue) const
  { (*params_).set(name,dvalue); return; }

  int                                                     level_;
  bool                                                    isnlncg_;
  ML*                                                     ml_;
  ML_Aggregate*                                           ag_;
  ML*                                                     thislevel_ml_;     // this level's 1-Level ML as a smoother to this level's nlnCG
  ML_Aggregate*                                           thislevel_ag_;     // this level's 1-Level ML as a smoother to this level's nlnCG
  RefCountPtr<Teuchos::ParameterList>                     params_;
  RefCountPtr< std::vector< RefCountPtr<Epetra_CrsMatrix> > >  P_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface>         fineinterface_;
  const Epetra_Comm&                                      comm_;
  RefCountPtr<Epetra_CrsMatrix>                           fineJac_;
  RefCountPtr<Epetra_CrsMatrix>                           SmootherA_;
  RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface>       coarseinterface_;
  RefCountPtr<NLNML::NLNML_PrePostOperator>               prepost_;
  RefCountPtr<Epetra_Vector>                              xthis_;
  RefCountPtr<NLNML::NLNML_ConstrainedMultiLevelOperator> thislevel_prec_;
  RefCountPtr<Teuchos::ParameterList>                     nlparams_;
  RefCountPtr<NOX::Epetra::MatrixFree>                    thislevel_A_;
  RefCountPtr<NLNML::NLNML_LinearSystem>                  thislevel_linSys_;
  RefCountPtr<NOX::Epetra::LinearSystemAztecOO>           azlinSys_;
  RefCountPtr<NOX::Epetra::Group>                         group_;
  //JJH
  //RefCountPtr<NOX::Solver::Manager>                       solver_;
  RCP<NOX::Solver::Generic>                               solver_;

  RefCountPtr<NOX::StatusTest::NormF>                     absresid_;
  RefCountPtr<NOX::StatusTest::NormUpdate>                nupdate_;
  RefCountPtr<NOX::StatusTest::Combo>                     combo1_;
  RefCountPtr<NOX::StatusTest::FiniteValue>               fv_;
  RefCountPtr<NOX::StatusTest::MaxIters>                  maxiters_;
  RefCountPtr<NOX::StatusTest::Combo>                     combo2_;

}; // class NLNML_NonlinearLevel
} // namespace NLNML

#endif
#endif
