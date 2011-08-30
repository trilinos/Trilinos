// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"
#include "Piro_Epetra_Factory.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#endif

#ifdef Piro_ENABLE_Rythmos
#include "Piro_Epetra_RythmosSolver.hpp"
#endif

#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#include "Piro_Epetra_VelocityVerletSolver.hpp"


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

  TEST_FOR_EXCEPTION(!found, Teuchos::Exceptions::InvalidParameter,
		     std::endl << "Error!  Piro::Epetra::Factory():  " <<
		     "Invalid Solver Type " << name << std::endl);

  return Teuchos::null;
}
