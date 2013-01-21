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

#include "Piro_ConfigDefs.hpp"
#include "Piro_Epetra_Factory.hpp"
#include "Teuchos_Assert.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#include "Piro_Epetra_VelocityVerletSolver.hpp"
#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#endif

#ifdef Piro_ENABLE_Rythmos
#include "Piro_Epetra_RythmosSolver.hpp"
#endif



Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::Factory::
createSolver(Teuchos::RCP<Teuchos::ParameterList> piroParams,
	     Teuchos::RCP<EpetraExt::ModelEvaluator> model)
{
  std::string name = piroParams->get("Solver Type", "NOX");
  bool found = false;

#ifdef Piro_ENABLE_NOX
  if (name == "NOX") {
    found = true;
    Teuchos::RCP<NOX::Epetra::Observer> observer =
      piroParams->get< Teuchos::RCP<NOX::Epetra::Observer> >(
	"Observer", Teuchos::null);
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface =
      piroParams->get< Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> >(
	"Interface", Teuchos::null);
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys =
      piroParams->get< Teuchos::RCP<NOX::Epetra::LinearSystem> >(
	"Linear System", Teuchos::null);
    return Teuchos::rcp(new Piro::Epetra::NOXSolver(
			  piroParams, model, observer, interface, linsys));
  }

  if (name == "LOCA") {
    found = true;
    Teuchos::RCP<NOX::Epetra::Observer> observer =
      piroParams->get< Teuchos::RCP<NOX::Epetra::Observer> >(
	"Observer", Teuchos::null);
    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData =
      piroParams->get< Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> >(
	"Save Eigen Data Strategy", Teuchos::null);
    Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTest =
      piroParams->get< Teuchos::RCP<LOCA::StatusTest::Abstract> >(
	"Status Test", Teuchos::null);
    return Teuchos::rcp(new Piro::Epetra::LOCASolver(
			  piroParams, model, observer, saveEigData,
			  locaStatusTest));
  }

  if (name == "Trapezoid Rule") {
    found = true;
    Teuchos::RCP<NOX::Epetra::Observer> observer =
      piroParams->get< Teuchos::RCP<NOX::Epetra::Observer> >(
	"Observer", Teuchos::null);
    return Teuchos::rcp(new Piro::Epetra::TrapezoidRuleSolver(
			  piroParams, model, observer));
  }

  if (name == "Velocity Verlet") {
    found = true;
    Teuchos::RCP<NOX::Epetra::Observer> observer =
      piroParams->get< Teuchos::RCP<NOX::Epetra::Observer> >(
	"Observer", Teuchos::null);
    return Teuchos::rcp(new Piro::Epetra::VelocityVerletSolver(
			  piroParams, model, observer));
  }
#endif

#ifdef Piro_ENABLE_Rythmos
  if (name == "Rythmos") {
    found = true;
    typedef Piro::Epetra::RythmosSolver::Scalar Scalar;
    Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer =
      piroParams->get< Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > >(
	"Observer", Teuchos::null);
    return Teuchos::rcp(new Piro::Epetra::RythmosSolver(
			  piroParams, model, observer));
  }
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(!found, Teuchos::Exceptions::InvalidParameter,
		     std::endl << "Error!  Piro::Epetra::Factory():  " <<
		     "Invalid Solver Type " << name << std::endl);

  return Teuchos::null;
}
