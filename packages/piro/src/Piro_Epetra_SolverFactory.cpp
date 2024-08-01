// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_SolverFactory.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#include "Piro_Epetra_LOCAAdaptiveSolver.hpp"
#include "Piro_Epetra_VelocityVerletSolver.hpp"
#include "Piro_Epetra_NewmarkSolver.hpp"
#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#endif

#include "Teuchos_TestForException.hpp"

#include <string>

namespace Piro {

namespace Epetra {

#ifdef HAVE_PIRO_NOX
template <>
std::string
SolverFactory::getLabel<NOX::Epetra::Observer>()
{
  return "NOX Observer";
}

template <>
Provider<NOX::Epetra::Observer> &
SolverFactory::getSource<NOX::Epetra::Observer>()
{
  return noxObserverSource_;
}


template <>
std::string
SolverFactory::getLabel<NOX::Epetra::ModelEvaluatorInterface>()
{
  return "Interface";
}

template <>
Provider<NOX::Epetra::ModelEvaluatorInterface> &
SolverFactory::getSource<NOX::Epetra::ModelEvaluatorInterface>()
{
  return noxInterfaceSource_;
}


template <>
std::string
SolverFactory::getLabel<NOX::Epetra::LinearSystem>()
{
  return "Linear System";
}

template <>
Provider<NOX::Epetra::LinearSystem> &
SolverFactory::getSource<NOX::Epetra::LinearSystem>()
{
  return noxLinearSystemSource_;
}


template <>
std::string
SolverFactory::getLabel<LOCA::SaveEigenData::AbstractStrategy>()
{
  return "Save Eigen Data Strategy";
}

template <>
Provider<LOCA::SaveEigenData::AbstractStrategy> &
SolverFactory::getSource<LOCA::SaveEigenData::AbstractStrategy>()
{
  return saveEigSource_;
}


template <>
std::string
SolverFactory::getLabel<LOCA::StatusTest::Abstract>()
{
  return "Status Test";
}

template <>
Provider<LOCA::StatusTest::Abstract> &
SolverFactory::getSource<LOCA::StatusTest::Abstract>()
{
  return statusTestSource_;
}


template <>
std::string
SolverFactory::getLabel<Piro::Epetra::AdaptiveSolutionManager>()
{
  return "Adaptive Solution Manager";
}

template <>
Provider<Piro::Epetra::AdaptiveSolutionManager> &
SolverFactory::getSource<Piro::Epetra::AdaptiveSolutionManager>()
{
  return adaptSolMgrSource_;
}
#endif /* HAVE_PIRO_NOX */

} // namespace Epetra

} // namespace Piro

template <typename T>
Teuchos::RCP<T>
Piro::Epetra::SolverFactory::create(const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  const std::string &label = this->getLabel<T>();
  return this->getSource<T>()(Teuchos::sublist(params, label));
}

Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::SolverFactory::createSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<EpetraExt::ModelEvaluator> &model)
{
  Teuchos::RCP<EpetraExt::ModelEvaluator> result;

  const std::string defaultType = "NOX";
  const std::string type = piroParams->get("Solver Type", defaultType);

#ifdef HAVE_PIRO_NOX
  if (type == "NOX") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    const Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface =
      this->create<NOX::Epetra::ModelEvaluatorInterface>(piroParams);
    const Teuchos::RCP<NOX::Epetra::LinearSystem> linsys =
      this->create<NOX::Epetra::LinearSystem>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::NOXSolver(piroParams, model, observer, interface, linsys));
  } else if (type == "LOCA") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    const Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData =
      this->create<LOCA::SaveEigenData::AbstractStrategy>(piroParams);
    const Teuchos::RCP<LOCA::StatusTest::Abstract> statusTest =
      this->create<LOCA::StatusTest::Abstract>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::LOCASolver(piroParams, model, observer, saveEigData, statusTest));
  } else if (type == "LOCA Adaptive") {
    const Teuchos::RCP<Piro::Epetra::AdaptiveSolutionManager> adaptMgr =
      this->create<Piro::Epetra::AdaptiveSolutionManager>(piroParams);
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    const Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData =
      this->create<LOCA::SaveEigenData::AbstractStrategy>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::LOCAAdaptiveSolver(piroParams, model, adaptMgr, observer, saveEigData));
  } else if (type == "Trapezoid Rule") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::TrapezoidRuleSolver(piroParams, model, observer));
  } else if (type == "Velocity Verlet") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::VelocityVerletSolver(piroParams, model, observer));
  } else if (type == "Newmark") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      this->create<NOX::Epetra::Observer>(piroParams);
    result = Teuchos::rcp(new Piro::Epetra::NewmarkSolver(piroParams, model, observer));
  } else
#endif

  {
    const bool solverTypeFound = false;
    TEUCHOS_TEST_FOR_EXCEPTION(
        !solverTypeFound, Teuchos::Exceptions::InvalidParameter,
        "Error in Piro::Epetra::SolverFactory(): Invalid Solver Type " << type);
  }

  return result;
}
