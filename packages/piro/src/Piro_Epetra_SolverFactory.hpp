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

#ifndef PIRO_EPETRA_SOLVERFACTORY_H
#define PIRO_EPETRA_SOLVERFACTORY_H

#include "EpetraExt_ModelEvaluator.h"

#include "Piro_ExtensibleFactory.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "NOX_Epetra_Observer.H"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_LinearSystem.H"

#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "LOCA_StatusTest_Abstract.H"

#include "Piro_Epetra_AdaptiveSolutionManager.hpp"
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
#include "Rythmos_IntegrationObserverBase.hpp"
#endif /* Piro_ENABLE_Rythmos */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <string>
#include <map>

namespace Piro {

namespace Epetra {

//! Factory for creating Epetra-based Piro solvers
class SolverFactory {
public:
  //! Create solver
  Teuchos::RCP<EpetraExt::ModelEvaluator>
  createSolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<EpetraExt::ModelEvaluator> &model);

  //! Replace the default auxiliary object provider
  template <typename T>
  void setDefaultProvider(const Provider<T> &p);

  //! Register a new auxiliary object provider
  template <typename T>
  void setProvider(const std::string &key, const Provider<T> &p);

#ifdef Piro_ENABLE_NOX
  //! Register a new NOX::Epetra::Observer provider
  void setNOXObserverProvider(
      const std::string &key,
      const Provider<NOX::Epetra::Observer> &p);

  //! Register a new NOX::Epetra::ModelEvaluatorInterface provider
  void setNOXInterfaceProvider(
      const std::string &key,
      const Provider<NOX::Epetra::ModelEvaluatorInterface> &p);

  //! Register a new NOX::Epetra::LinearSystem provider
  void setNOXLinearSystemProvider(
      const std::string &key,
      const Provider<NOX::Epetra::LinearSystem> &p);

  //! Register a new LOCA::SaveEigenData::AbstractStrategy provider
  void setLOCASaveEigenDataProvider(
      const std::string &key,
      const Provider<LOCA::SaveEigenData::AbstractStrategy> &p);

  //! Register a new LOCA::StatusTest::Abstract provider
  void setLOCAStatusTestProvider(
      const std::string &key,
      const Provider<LOCA::StatusTest::Abstract> &p);

  //! Register a new Piro::Epetra::AdaptiveSolutionManager provider
  void setAdapativeSolutionManagerProvider(
      const std::string &key,
      const Provider<Piro::Epetra::AdaptiveSolutionManager> &p);
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
  //! Register a new Rythmos::IntegrationObserverBase<double> provider
  void setRythmosObserverProvider(
      const std::string &key,
      const Provider<Rythmos::IntegrationObserverBase<double> > &p);
#endif /* Piro_ENABLE_Rythmos */

private:
  template <typename T>
  static std::string getLabel();

  template <typename T>
  Teuchos::RCP<T> create(const Teuchos::RCP<Teuchos::ParameterList> &params);

  template <typename T>
  ExtensibleFactory<T> &getFactory();

#ifdef Piro_ENABLE_NOX
  ExtensibleFactory<NOX::Epetra::Observer> noxObserverFactory_;
  ExtensibleFactory<NOX::Epetra::ModelEvaluatorInterface> noxInterfaceFactory_;
  ExtensibleFactory<NOX::Epetra::LinearSystem> noxLinearSystemFactory_;
  ExtensibleFactory<LOCA::SaveEigenData::AbstractStrategy> saveEigFactory_;
  ExtensibleFactory<LOCA::StatusTest::Abstract> statusTestFactory_;
  ExtensibleFactory<Piro::Epetra::AdaptiveSolutionManager> adaptSolMgrFactory_;
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
  ExtensibleFactory<Rythmos::IntegrationObserverBase<double> > rythmosObserverFactory_;
#endif /* Piro_ENABLE_Rythmos */
};

template <typename T>
inline
void
SolverFactory::setDefaultProvider(const Piro::Provider<T> &p)
{
  this->getFactory<T>().setDefaultProvider(p);
}

template <typename T>
inline
void
SolverFactory::setProvider(const std::string &key, const Piro::Provider<T> &p)
{
  this->getFactory<T>().setProvider(key, p);
}

} // namespace Epetra

} // namespace Piro

#endif /* PIRO_EPETRA_SOLVERFACTORY_H */
