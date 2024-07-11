// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
           const int number_of_time_derivatives,
           const Teuchos::RCP<const Epetra_Map>& map_,
           const Teuchos::RCP<const Epetra_Map>& overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph>& overlapJacGraph_);

  virtual ~AdaptiveSolutionManager () {}

  virtual Teuchos::RCP<const Epetra_Vector> updateSolution();

  //! Build the LOCA solution group
  virtual Teuchos::RCP<LOCA::Epetra::Group>
       buildSolutionGroup() = 0;

  //! Remap "old" solution into new data structures
  virtual void projectCurrentSolution();

  //! Accessor functions

  virtual void getConvergenceData(int& KrylovIters, int& lastSolveKrylocIters, int& linSolves, double& tolAchieved) const;

  virtual Teuchos::RCP<Teuchos::ParameterList> getAdaptParamsNonConst() { return adaptParams; }

  virtual Teuchos::RCP<const Teuchos::ParameterList> getAdaptParams() const { return adaptParams; }

protected:

      //! The solution group
      Teuchos::RCP<LOCA::Epetra::Group> grp;

      //! The linear system used for the solution in the group
      Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;

      //! The adaptation parameter list
      Teuchos::RCP<Teuchos::ParameterList> adaptParams;

      //! The parent of the solution parameter list
      Teuchos::RCP<Teuchos::ParameterList> piroParams;

};

}
}

#endif //LOCA_EPETRA_ADAPTIVESOLUTIONMANAGER
