// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_NOXSOLVER_H
#define PIRO_EPETRA_NOXSOLVER_H

#include <iostream>

#include "NOX.H"
#include "NOX_Epetra.H"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"
#include <NOX_Epetra_MultiVector.H>
#include <NOX_Epetra_Observer.H>

#include "LOCA_GlobalData.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"

namespace Piro {
namespace Epetra {

/** \brief Epetra-based Model Evaluator for NOX solves
 *  \ingroup Piro_Epetra_solver_grp
 * */
class NOXSolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  NOXSolver(Teuchos::RCP<Teuchos::ParameterList> piroParams,
            Teuchos::RCP<EpetraExt::ModelEvaluator> model,
            Teuchos::RCP<NOX::Epetra::Observer> observer = Teuchos::null,
	    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface = 
	    Teuchos::null,
	    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys = Teuchos::null
            );


  //@}

  ~NOXSolver();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_lower_bounds(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_upper_bounds(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_DgDp_op( int j, int l ) const;
  /** \brief . */
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  void setProblemParamDefaults(Teuchos::ParameterList* piroParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* piroParams_);

  //@}

  void resetCounters();

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> piroParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<NOX::Epetra::Observer> observer;
   NOX::Utils utils;

   Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
   Teuchos::RCP<NOX::Solver::Generic> solver;
   Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
   Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface;
   int num_p;
   int num_g;

   Teuchos::RCP<NOX::Epetra::Group> grp;
   Teuchos::RCP<LOCA::GlobalData> globalData;
   Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy;
   bool writeOnlyConvergedSol;
   bool exitUponFailedNOXSolve; 

   //Store current iteration of Analysis solver 
   mutable int current_iteration;

   enum DerivativeLayout { OP, COL, ROW };

  mutable int totalNewtonIters;
  mutable int totalKrylovIters;
  mutable int stepNum;
};

}
}
#endif
