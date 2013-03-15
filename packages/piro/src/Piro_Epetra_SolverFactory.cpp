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

namespace Piro {

namespace Epetra {

#ifdef Piro_ENABLE_NOX
template <>
std::string
SolverFactory::getLabel<NOX::Epetra::Observer>()
{
  return "NOX Observer";
}

template <>
ExtensibleFactory<NOX::Epetra::Observer> &
SolverFactory::getFactory<NOX::Epetra::Observer>()
{
  return noxObserverFactory_;
}


template <>
std::string
SolverFactory::getLabel<NOX::Epetra::ModelEvaluatorInterface>()
{
  return "Interface";
}

template <>
ExtensibleFactory<NOX::Epetra::ModelEvaluatorInterface> &
SolverFactory::getFactory<NOX::Epetra::ModelEvaluatorInterface>()
{
  return noxInterfaceFactory_;
}


template <>
std::string
SolverFactory::getLabel<NOX::Epetra::LinearSystem>()
{
  return "Linear System";
}

template <>
ExtensibleFactory<NOX::Epetra::LinearSystem> &
SolverFactory::getFactory<NOX::Epetra::LinearSystem>()
{
  return noxLinearSystemFactory_;
}


template <>
std::string
SolverFactory::getLabel<LOCA::SaveEigenData::AbstractStrategy>()
{
  return "Save Eigen Data Strategy";
}

template <>
ExtensibleFactory<LOCA::SaveEigenData::AbstractStrategy> &
SolverFactory::getFactory<LOCA::SaveEigenData::AbstractStrategy>()
{
  return saveEigFactory_;
}


template <>
std::string
SolverFactory::getLabel<LOCA::StatusTest::Abstract>()
{
  return "Status Test";
}

template <>
ExtensibleFactory<LOCA::StatusTest::Abstract> &
SolverFactory::getFactory<LOCA::StatusTest::Abstract>()
{
  return statusTestFactory_;
}


template <>
std::string
SolverFactory::getLabel<Piro::Epetra::AdaptiveSolutionManager>()
{
  return "Adaptive Solution Manager";
}

template <>
ExtensibleFactory<Piro::Epetra::AdaptiveSolutionManager> &
SolverFactory::getFactory<Piro::Epetra::AdaptiveSolutionManager>()
{
  return adaptSolMgrFactory_;
}
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
template <>
std::string
SolverFactory::getLabel<Rythmos::IntegrationObserverBase<double> >()
{
  return "Rythmos Observer";
}

template <>
ExtensibleFactory<Rythmos::IntegrationObserverBase<double> > &
SolverFactory::getFactory<Rythmos::IntegrationObserverBase<double> >()
{
  return rythmosObserverFactory_;
}
#endif /* Piro_ENABLE_Rythmos */

} // namespace Epetra

} // namespace Piro

template <typename T>
Teuchos::RCP<T>
Piro::Epetra::SolverFactory::create(const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  const std::string &label = this->getLabel<T>();
  return this->getFactory<T>().create(Teuchos::sublist(params, label));
}


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
  } else
#endif

#ifdef Piro_ENABLE_Rythmos
  if (type == "Rythmos") {
    const Teuchos::RCP<Rythmos::IntegrationObserverBase<double> > observer =
      this->create<Rythmos::IntegrationObserverBase<double> >(piroParams);
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
  this->setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setNOXInterfaceProvider(
    const std::string &key,
    const Piro::Provider<NOX::Epetra::ModelEvaluatorInterface> &p)
{
  this->setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setNOXLinearSystemProvider(
    const std::string &key,
    const Piro::Provider<NOX::Epetra::LinearSystem> &p)
{
  this->setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setLOCASaveEigenDataProvider(
    const std::string &key,
    const Piro::Provider<LOCA::SaveEigenData::AbstractStrategy> &p)
{
  this->setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setLOCAStatusTestProvider(
    const std::string &key,
    const Piro::Provider<LOCA::StatusTest::Abstract> &p)
{
  this->setProvider(key, p);
}

void
Piro::Epetra::SolverFactory::setAdapativeSolutionManagerProvider(
    const std::string &key,
    const Provider<Piro::Epetra::AdaptiveSolutionManager> &p)
{
  this->setProvider(key, p);
}
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
void
Piro::Epetra::SolverFactory::setRythmosObserverProvider(
    const std::string &key,
    const Piro::Provider<Rythmos::IntegrationObserverBase<double> > &p)
{
  this->setProvider(key, p);
}
#endif /* Piro_ENABLE_Rythmos */
