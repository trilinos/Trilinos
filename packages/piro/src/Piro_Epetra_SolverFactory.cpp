// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_Epetra_SolverFactory.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#include "Piro_Epetra_LOCAAdaptiveSolver.hpp"
#include "Piro_Epetra_VelocityVerletSolver.hpp"
#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#endif

#ifdef Piro_ENABLE_Rythmos
#include "Piro_Epetra_RythmosSolver.hpp"
#endif

#include "Teuchos_TestForException.hpp"

#include <string>

Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::SolverFactory::createSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<EpetraExt::ModelEvaluator> &model)
{
  Teuchos::RCP<EpetraExt::ModelEvaluator> result;

  const std::string defaultType = "NOX";
  const std::string type = piroParams->get("Solver Type", defaultType);

#ifdef Piro_ENABLE_NOX
  if (type == "NOX") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      noxObserverFactory_.create(Teuchos::sublist(piroParams, "NOX Observer"));
    const Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface =
      noxInterfaceFactory_.create(Teuchos::sublist(piroParams, "Interface"));
    const Teuchos::RCP<NOX::Epetra::LinearSystem> linsys =
      noxLinearSystemFactory_.create(Teuchos::sublist(piroParams, "Linear System"));
    result = Teuchos::rcp(new Piro::Epetra::NOXSolver(piroParams, model, observer, interface, linsys));
  } else if (type == "LOCA") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      noxObserverFactory_.create(Teuchos::sublist(piroParams, "NOX Observer"));
    const Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData =
      saveEigFactory_.create(Teuchos::sublist(piroParams, "Save Eigen Data Strategy"));
    const Teuchos::RCP<LOCA::StatusTest::Abstract> statusTest =
      statusTestFactory_.create(Teuchos::sublist(piroParams, "Status Test"));
    result = Teuchos::rcp(new Piro::Epetra::LOCASolver(piroParams, model, observer, saveEigData, statusTest));
  } else if (type == "LOCA Adaptive") {
    const Teuchos::RCP<Piro::Epetra::AdaptiveSolutionManager> adaptMgr =
      adaptSolMgrFactory_.create(Teuchos::sublist(piroParams, "Adaptive Solution Manager"));
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      noxObserverFactory_.create(Teuchos::sublist(piroParams, "NOX Observer"));
    const Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData =
      saveEigFactory_.create(Teuchos::sublist(piroParams, "Save Eigen Data Strategy"));
    result = Teuchos::rcp(new Piro::Epetra::LOCAAdaptiveSolver(piroParams, model, adaptMgr, observer, saveEigData));
  } else if (type == "Trapezoid Rule") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      noxObserverFactory_.create(Teuchos::sublist(piroParams, "NOX Observer"));
    result = Teuchos::rcp(new Piro::Epetra::TrapezoidRuleSolver(piroParams, model, observer));
  } else if (type == "Velocity Verlet") {
    const Teuchos::RCP<NOX::Epetra::Observer> observer =
      noxObserverFactory_.create(Teuchos::sublist(piroParams, "NOX Observer"));
    result = Teuchos::rcp(new Piro::Epetra::VelocityVerletSolver(piroParams, model, observer));
  } else
#endif

#ifdef Piro_ENABLE_Rythmos
  if (type == "Rythmos") {
    const Teuchos::RCP<Rythmos::IntegrationObserverBase<double> > observer =
      rythmosObserverFactory_.create(Teuchos::sublist(piroParams, "Rythmos Observer"));
    result = Teuchos::rcp(new Piro::Epetra::RythmosSolver(piroParams, model, observer));
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

#ifdef Piro_ENABLE_NOX
void
Piro::Epetra::SolverFactory::setNOXObserverProvider(
    const std::string &key,
    const Piro::Provider<NOX::Epetra::Observer> &p)
{
  noxObserverFactory_.setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setNOXInterfaceProvider(
    const std::string &key,
    const Piro::Provider<NOX::Epetra::ModelEvaluatorInterface> &p)
{
  noxInterfaceFactory_.setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setNOXLinearSystemProvider(
    const std::string &key,
    const Piro::Provider<NOX::Epetra::LinearSystem> &p)
{
  noxLinearSystemFactory_.setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setLOCASaveEigenDataProvider(
    const std::string &key,
    const Piro::Provider<LOCA::SaveEigenData::AbstractStrategy> &p)
{
  saveEigFactory_.setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setLOCAStatusTestProvider(
    const std::string &key,
    const Piro::Provider<LOCA::StatusTest::Abstract> &p)
{
  statusTestFactory_.setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setAdapativeSolutionManagerProvider(
    const std::string &key,
    const Provider<Piro::Epetra::AdaptiveSolutionManager> &p)
{
  adaptSolMgrFactory_.setProvider(key, p);
}
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
void
Piro::Epetra::SolverFactory::setRythmosObserverProvider(
    const std::string &key,
    const Piro::Provider<Rythmos::IntegrationObserverBase<double> > &p)
{
  rythmosObserverFactory_.setProvider(key, p);
}
#endif /* Piro_ENABLE_Rythmos */
