// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_LOCASOLVER_H
#define PIRO_EPETRA_LOCASOLVER_H

#include <iostream>

#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"
#include "LOCA_SaveEigenData_AbstractStrategy.H"

#include <NOX_Epetra_MultiVector.H>
#include <NOX_Epetra_Observer.H>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"

namespace Piro {
namespace Epetra {

/** \brief Epetra-based Model Evaluator for LOCA solves
 *  \ingroup Piro_Epetra_solver_grp
 * */
class LOCASolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  LOCASolver(Teuchos::RCP<Teuchos::ParameterList> piroParams,
            Teuchos::RCP<EpetraExt::ModelEvaluator> model,
            Teuchos::RCP<NOX::Epetra::Observer> observer = Teuchos::null,
            Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData = Teuchos::null,
            Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTest_ = Teuchos::null
            );


  //@}

  ~LOCASolver();


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
//  Teuchos::RCP<Epetra_Operator> create_W() const;
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
  
  public:
  
  /** \brief \c const getter for the LOCA stepper. Can be used to retrieve statistics. */
  Teuchos::RCP<const LOCA::Stepper> getLOCAStepper() const;
  /** \brief Nonconst getter for the LOCA stepper. Can be used to adapt the stepper configuration while running. */
  Teuchos::RCP<LOCA::Stepper> getLOCAStepperNonConst();
  
  /** \brief \c const getter for the \c LOCA::GlobalData object.*/
  Teuchos::RCP<const LOCA::GlobalData> getGlobalData() const;
  /** \brief Nonconst getter for the \c LOCA::GlobalData object.*/
  Teuchos::RCP<LOCA::GlobalData> getGlobalDataNonConst();

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> piroParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<NOX::Epetra::Observer> observer;
   Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData;
   Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTest;
   NOX::Utils utils;

   Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
   Teuchos::RCP<LOCA::Stepper> stepper;
   Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
   Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;
   Teuchos::RCP<LOCA::ParameterVector> pVector;
   int num_p;
   int num_g;

   Teuchos::RCP<LOCA::Epetra::Group> grp;
   Teuchos::RCP<LOCA::GlobalData> globalData;

  mutable int totalNewtonIters;
  mutable int totalKrylovIters;
  mutable int totalLinSolves;
  mutable int stepNum;
};

}
}
#endif
