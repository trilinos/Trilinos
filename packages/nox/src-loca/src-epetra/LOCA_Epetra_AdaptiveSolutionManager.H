//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2012) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact 
// Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_ADAPTIVESOLUTIONMANAGER
#define LOCA_EPETRA_ADAPTIVESOLUTIONMANAGER

#include "LOCA.H"
#include "LOCA_Epetra.H"

#include "NOX_Epetra_AdaptiveSolutionManager.H"


namespace LOCA {
namespace Epetra {

class AdaptiveSolutionManager : public NOX::Epetra::AdaptiveSolutionManager
{
public:

  AdaptiveSolutionManager (
           const Teuchos::RCP<Teuchos::ParameterList>& params,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_);

  virtual ~AdaptiveSolutionManager () {}

  //! Get initial solution

  virtual Teuchos::RCP<LOCA::Epetra::Group>& getSolutionGroup(){ return grp;}

  //! Build the LOCA solution group
  virtual Teuchos::RCP<LOCA::Epetra::Group> 
       buildSolutionGroup();

  virtual void projectCurrentSolution();

  virtual void initialize(
       Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
        Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface,
       Teuchos::RCP<Teuchos::ParameterList>& piroParams,
       Teuchos::RCP<LOCA::ParameterVector>& pVector,
       Teuchos::RCP<LOCA::GlobalData>& globalData,
       bool createPrec);

  //! Accessor functions

  virtual Teuchos::RCP<NOX::Epetra::LinearSystem>& getLinearSystem(){ return linsys;}

  virtual Teuchos::RCP<LOCA::GlobalData> getGlobalDataNonConst() { return globalData; }

  virtual Teuchos::RCP<const LOCA::GlobalData> getGlobalData() const { return globalData; }

protected:
 
      //! The solution group
      Teuchos::RCP<LOCA::Epetra::Group> grp;

      //! The interface to the model evaluator
      Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;

      //! The LOCA parameter set
      Teuchos::RCP<LOCA::ParameterVector> pVector;

      //! The LOCA global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! The linear system used for the solution in the group
      Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;

      //! The model we are solving
      Teuchos::RCP<EpetraExt::ModelEvaluator> model;

};

}
}

#endif //LOCA_EPETRA_ADAPTIVESOLUTIONMANAGER
