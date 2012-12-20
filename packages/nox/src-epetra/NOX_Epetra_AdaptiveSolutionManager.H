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

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_x(){return overlapped_x;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_xdot(){return overlapped_xdot;}

  virtual Teuchos::RCP<Epetra_Vector>& get_overlapped_f(){return overlapped_f;}

  virtual Teuchos::RCP<Epetra_CrsMatrix>& get_overlapped_jac(){return overlapped_jac;}

  virtual Teuchos::RCP<Epetra_Import>& get_importer(){return importer;}

  virtual Teuchos::RCP<Epetra_Export>& get_exporter(){return exporter;}


protected:

     //! Initial solution vector
      Teuchos::RCP<Epetra_Vector> initial_x;
  
      //! Initial solution vector
      Teuchos::RCP<Epetra_Vector> initial_xdot;
  
      //! Importer for overlapped data
      Teuchos::RCP<Epetra_Import> importer;
  
      //! Exporter for overlapped data
      Teuchos::RCP<Epetra_Export> exporter;
  
      //! Overlapped solution vector
      Teuchos::RCP<Epetra_Vector> overlapped_x;
  
      //! Overlapped time derivative vector
      Teuchos::RCP<Epetra_Vector> overlapped_xdot;
  
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


};

}
}

#endif //NOX_EPETRA_ADAPTIVESOLUTIONMANAGER
