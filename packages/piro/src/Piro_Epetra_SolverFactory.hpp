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

#include "Piro_Provider.hpp"

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

  //! Set the auxiliary object provider of type T
  template <typename T>
  void setSource(const Provider<T> &p);

  //! Deprecated - Set the auxiliary object provider of type T
  template <typename T>
  void setDefaultProvider(const Provider<T> &p);

private:
  template <typename T>
  static std::string getLabel();

  template <typename T>
  Teuchos::RCP<T> create(const Teuchos::RCP<Teuchos::ParameterList> &params);

  template <typename T>
  Provider<T> &getSource();

#ifdef Piro_ENABLE_NOX
  Provider<NOX::Epetra::Observer> noxObserverSource_;
  Provider<NOX::Epetra::ModelEvaluatorInterface> noxInterfaceSource_;
  Provider<NOX::Epetra::LinearSystem> noxLinearSystemSource_;
  Provider<LOCA::SaveEigenData::AbstractStrategy> saveEigSource_;
  Provider<LOCA::StatusTest::Abstract> statusTestSource_;
  Provider<Piro::Epetra::AdaptiveSolutionManager> adaptSolMgrSource_;
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
  Provider<Rythmos::IntegrationObserverBase<double> > rythmosObserverSource_;
#endif /* Piro_ENABLE_Rythmos */
};

template <typename T>
inline
void
SolverFactory::setSource(const Piro::Provider<T> &p)
{
  this->getSource<T>() = p;
}

template <typename T>
inline
void
SolverFactory::setDefaultProvider(const Piro::Provider<T> &p)
{
  this->setSource(p);
}

} // namespace Epetra

} // namespace Piro

#endif /* PIRO_EPETRA_SOLVERFACTORY_H */
