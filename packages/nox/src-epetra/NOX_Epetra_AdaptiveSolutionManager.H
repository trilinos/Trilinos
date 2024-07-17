// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_ADAPTIVESOLUTIONMANAGER
#define NOX_EPETRA_ADAPTIVESOLUTIONMANAGER

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Epetra_Vector.H"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "NOX_Epetra_AdaptManager.H"

namespace NOX {
namespace Epetra {

class AdaptiveSolutionManager
{
public:

  AdaptiveSolutionManager (
           const int number_of_time_derivatives,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_);

  virtual ~AdaptiveSolutionManager () {}

  //! Apply adaptation method to mesh and problem. Returns true if adaptation is performed successfully.
  virtual bool adaptProblem() = 0;

  //! Accessor function to Epetra_Import the solution from other PEs for output
  virtual Epetra_Vector* getOverlapSolution(const Epetra_Vector& solution);

  //! Function to get the active AdaptManager
  virtual Teuchos::RCP<NOX::Epetra::AdaptManager> getAdaptManager()
        { return adaptManager; }

  //! Function to resize the data structures when adapting
  virtual void resizeMeshDataArrays(
        const Teuchos::RCP<const Epetra_Map> &map,
        const Teuchos::RCP<const Epetra_Map> &overlapMap,
        const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph);


  //! Accessor functions

  virtual Teuchos::RCP<const Epetra_Vector> getInitialSolution() const { return initial_x; }

  virtual Teuchos::RCP<Epetra_Vector>& get_initial_x(){return initial_x;}

  virtual void setInitialSolution(const Teuchos::RCP<Epetra_Vector>& soln_);

  virtual Teuchos::RCP<const Epetra_Vector> getInitialSolutionDot() const { return initial_xdot; }

  virtual Teuchos::RCP<Epetra_Vector>& get_initial_xdot(){return initial_xdot;}

  virtual Teuchos::RCP<const Epetra_Vector> getInitialSolutionDotDot() const { return initial_xdotdot; }

  virtual Teuchos::RCP<Epetra_Vector>& get_initial_xdotdot(){return initial_xdotdot;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_x(){return overlapped_x;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_xdot(){return overlapped_xdot;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_xdotdot(){return overlapped_xdotdot;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_f(){return overlapped_f;}

  virtual Teuchos::RCP<Epetra_CrsMatrix>& get_overlapped_jac(){return overlapped_jac;}

  virtual Teuchos::RCP<Epetra_Import>& get_importer(){return importer;}

  virtual Teuchos::RCP<Epetra_Export>& get_exporter(){return exporter;}


protected:

     //! Initial solution vector
      Teuchos::RCP<Epetra_Vector> initial_x;

      //! Initial solution vector
      Teuchos::RCP<Epetra_Vector> initial_xdot;

      //! Initial solution vector
      Teuchos::RCP<Epetra_Vector> initial_xdotdot;

      //! Importer for overlapped data
      Teuchos::RCP<Epetra_Import> importer;

      //! Exporter for overlapped data
      Teuchos::RCP<Epetra_Export> exporter;

      //! Overlapped solution vector
      Teuchos::RCP<Epetra_Vector> overlapped_x;

      //! Overlapped time derivative vector
      Teuchos::RCP<Epetra_Vector> overlapped_xdot;

      //! Overlapped time acceleration vector
      Teuchos::RCP<Epetra_Vector> overlapped_xdotdot;

      //! Overlapped residual vector
      Teuchos::RCP<Epetra_Vector> overlapped_f;

      //! Overlapped Jacobian matrix
      Teuchos::RCP<Epetra_CrsMatrix> overlapped_jac;

      //! Temporary overlapped solution vector (for output)
      Teuchos::RCP<Epetra_Vector> tmp_ovlp_sol;

      //! The adaptation manager - makes the call to adapt the mesh
      Teuchos::RCP<NOX::Epetra::AdaptManager> adaptManager;

      //! Current solution vector
      Teuchos::RCP<NOX::Epetra::Vector> currentSolution;

      //! Number of time derivatives (x - 0, x_dot - 1, x_dotdot - 2)
      const int num_time_deriv;


};

}
}

#endif //NOX_EPETRA_ADAPTIVESOLUTIONMANAGER
