// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_ADAPTIVESOLUTIONMANAGER_H
#define PIRO_EPETRA_ADAPTIVESOLUTIONMANAGER_H

#include "LOCA_Epetra_AdaptiveSolutionManager.H"


namespace Piro {
namespace Epetra {

  class AdaptiveSolutionManager : public LOCA::Epetra::AdaptiveSolutionManager {
  
  public:
     AdaptiveSolutionManager (const Teuchos::RCP<Teuchos::ParameterList>& appParams,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_);

    virtual ~AdaptiveSolutionManager();

    virtual void initialize(
       const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
       const Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface>& interface,
       const Teuchos::RCP<LOCA::ParameterVector>& pVector,
       const Teuchos::RCP<LOCA::GlobalData>& globalData,
       bool createPrec);

    //! Returns true if the user has configured a mesh adaptive problem
    virtual bool hasAdaptation(){ return adaptiveMesh; }

    //! Build the LOCA solution group
    virtual Teuchos::RCP<LOCA::Epetra::Group> 
       buildSolutionGroup();

    //! Cleanup and destroy solution group and support infrastructure
    virtual void 
       destroySolutionGroup();

    virtual void applyJacobianInverseMultiVector(
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result)
                        { grp->computeJacobian();
                          grp->applyJacobianInverseMultiVector(*piroParams, input, result); }

  protected:

      //! The interface to the model evaluator
      Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;

      //! The LOCA parameter set
      Teuchos::RCP<LOCA::ParameterVector> pVector;

      //! The LOCA global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! The model we are solving
      Teuchos::RCP<EpetraExt::ModelEvaluator> model;

      //! Preconditioning flag
      bool createPrec;

      //! Mesh adaptivity flag
      bool adaptiveMesh;

  };

}
}

#endif //PIRO_EPETRA_ADAPTIVESOLUTIONMANAGER


