// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_LOCAADAPTIVESOLVER_H
#define PIRO_EPETRA_LOCAADAPTIVESOLVER_H

#include <iostream>

#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"
#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "LOCA_Epetra_AdaptiveStepper.H"


#include <NOX_Epetra_MultiVector.H>
#include <NOX_Epetra_Observer.H>
#include <Piro_Epetra_AdaptiveSolutionManager.hpp>


#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"

namespace Piro {
namespace Epetra {

/** \ingroup Piro_Epetra_solver_grp */
class LOCAAdaptiveSolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  LOCAAdaptiveSolver(const Teuchos::RCP<Teuchos::ParameterList>& piroParams,
            const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
            const Teuchos::RCP<Piro::Epetra::AdaptiveSolutionManager>& solnManager,
            Teuchos::RCP<NOX::Epetra::Observer> observer = Teuchos::null,
            Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData = Teuchos::null
            );


  //@}

  ~LOCAAdaptiveSolver();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
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
  
  
  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Teuchos::ParameterList> piroParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<NOX::Epetra::Observer> observer;
   Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData;
   Teuchos::RCP<Piro::Epetra::AdaptiveSolutionManager> solnManager;
   NOX::Utils utils;

   Teuchos::RCP<LOCA::Epetra::AdaptiveStepper> stepper;
   Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;
   Teuchos::RCP<LOCA::ParameterVector> pVector;
   Teuchos::RCP<LOCA::GlobalData> globalData;

   int num_p;
   int num_g;

};

}
}
#endif
