// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_SOLVERFACTORY_H
#define PIRO_EPETRA_SOLVERFACTORY_H

#include "EpetraExt_ModelEvaluator.h"

#include "Piro_Provider.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "NOX_Epetra_Observer.H"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_LinearSystem.H"

#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "LOCA_StatusTest_Abstract.H"

#include "Piro_Epetra_AdaptiveSolutionManager.hpp"
#endif /* HAVE_PIRO_NOX */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <string>
#include <map>

namespace Piro {

namespace Epetra {

/*! \brief Factory for creating Epetra-based %Piro solvers
 *
 *  Piro::SolverFactory is the counterpart for Thyra-based models.
 */
class SolverFactory {
public:
  /*! \brief Create a solved model
   *
   *  The type of %Piro solver to instantiate is determined by the value of the string entry labeled <tt>"Solver Type"</tt>
   *  and located at the top level of parameter list \c piroParams.
   *
   *  Currently, the following solver types are available (each accompanied by the corresponding token value):
   *  - Piro::Epetra::NOXSolver (<tt>"NOX"</tt>)
   *  - Piro::Epetra::LOCASolver (<tt>"LOCA"</tt>)
   *  - Piro::Epetra::LOCAAdaptiveSolver (<tt>"LOCA Adaptive"</tt>)
   *  - Piro::Epetra::TrapezoidRuleSolver (<tt>"Trapezoid Rule"</tt>)
   *  - Piro::Epetra::VelocityVerletSolver (<tt>"Velocity Verlet"</tt>)
   */
  Teuchos::RCP<EpetraExt::ModelEvaluator>
  createSolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<EpetraExt::ModelEvaluator> &model);

  /*! \brief Set the source for auxiliary objects
   *
   *  While all %Piro solvers only require a model evaluator to be properly initialized,
   *  some of them also optionally accept so-called `auxiliary' objects which extend their functionality.
   *  For instance, observers allow to monitor the solution process and the intermediate steps taken by the solver.
   *
   *  Currently, the following types are acceptable values for the generic type \c T:
   *  - NOX::Epetra::Observer
   *  - NOX::Epetra::ModelEvaluatorInterface
   *  - NOX::Epetra::LinearSystem
   *  - LOCA::SaveEigenData::AbstractStrategy
   *  - LOCA::StatusTest::Abstract
   *  - Piro::Epetra::AdaptiveSolutionManager
   */
  template <typename T>
  void setSource(const Provider<T> &p);

private:
  template <typename T>
  static std::string getLabel();

  template <typename T>
  Teuchos::RCP<T> create(const Teuchos::RCP<Teuchos::ParameterList> &params);

  template <typename T>
  Provider<T> &getSource();

#ifdef HAVE_PIRO_NOX
  Provider<NOX::Epetra::Observer> noxObserverSource_;
  Provider<NOX::Epetra::ModelEvaluatorInterface> noxInterfaceSource_;
  Provider<NOX::Epetra::LinearSystem> noxLinearSystemSource_;
  Provider<LOCA::SaveEigenData::AbstractStrategy> saveEigSource_;
  Provider<LOCA::StatusTest::Abstract> statusTestSource_;
  Provider<Piro::Epetra::AdaptiveSolutionManager> adaptSolMgrSource_;
#endif /* HAVE_PIRO_NOX */
};

template <typename T>
inline
void
SolverFactory::setSource(const Piro::Provider<T> &p)
{
  this->getSource<T>() = p;
}

} // namespace Epetra

} // namespace Piro

#endif /* PIRO_EPETRA_SOLVERFACTORY_H */
