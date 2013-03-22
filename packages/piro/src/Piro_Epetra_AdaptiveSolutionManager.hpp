// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Glen Hansen (gahanse@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
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


