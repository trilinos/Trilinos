//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"

#include "../example/LorenzModel.hpp"
#include "Rythmos_GAASPErrorEstimator.hpp"
#include "Rythmos_ErrorEstimateBase.hpp"
#include "Rythmos_GAASPInterface.hpp"
#include "Rythmos_Types.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_GAASPErrorEstimator, defaultConstruct ) {
  GAASPErrorEstimator gaaspEE;
  Teuchos::RCP<const ErrorEstimateBase<double> > errorEstimate;
  TEST_THROW( errorEstimate = gaaspEE.getErrorEstimate(), std::logic_error ); 
}

TEUCHOS_UNIT_TEST( Rythmos_GAASPInterface, defaultConstruct ) {
  GAASPInterface gI;
  TEST_THROW( gI.forwardSolve(), std::logic_error );
}

RCP<Rythmos::GAASPErrorEstimator> createGAASPEELorenz(double dt) {
  RCP<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = Teuchos::rcp( new Epetra_SerialComm );
#endif // HAVE_MPI
  RCP<ParameterList>
    paramList = Teuchos::parameterList();
  RCP<LorenzModel>
    lorenzModel = Teuchos::rcp(new LorenzModel( epetra_comm, *paramList ));
  RCP<Thyra::ModelEvaluator<double> >
    thyraLorenzModel = Thyra::epetraModelEvaluator(lorenzModel,Teuchos::null);
  RCP<Rythmos::GAASPErrorEstimator> gaaspEE = Teuchos::rcp(new Rythmos::GAASPErrorEstimator);
  gaaspEE->setModel(thyraLorenzModel);
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->sublist("GAASP Interface Parameters").set<double>("eTime",1.0);
  pl->sublist("GAASP Interface Parameters").set<double>("timeStep",dt);
  RCP<ParameterList> voPL = Teuchos::sublist(pl,"VerboseObject");
  voPL->set("Verbosity Level","none");
  RCP<ParameterList> gaaspIPL = Teuchos::sublist(pl,"GAASP Interface Parameters");
  RCP<ParameterList> gaaspIPLvopl = Teuchos::sublist(gaaspIPL,"VerboseObject");
  gaaspIPLvopl->set("Verbosity Level", "none");
  gaaspEE->setParameterList(pl);
  return gaaspEE;
}

// This is a regression test.  The behavior has changed over time.... (1/20/09)
TEUCHOS_UNIT_TEST( Rythmos_GAASPInterface, lorenz ) {
  double dt = 0.1;
  RCP<Rythmos::GAASPErrorEstimator> gaaspEE = createGAASPEELorenz(dt);
  RCP<const Rythmos::ErrorEstimateBase<double> > errObj = gaaspEE->getErrorEstimate();
  double err = errObj->getTotalError();
  double tol = 1.0e-14;
  TEST_FLOATING_EQUALITY( err, -0.0296387721978618, tol );
}

// This is a regression test also.
//TEUCHOS_UNIT_TEST( Rythmos_GAASPInterface, lorenz_controlError ) {
//  double dt = 0.1;
//  RCP<Rythmos::GAASPErrorEstimator> gaaspEE = createGAASPEELorenz(dt);
//  double uTOL = 1.0e-2;
//  RCP<const Rythmos::ErrorEstimateBase<double> > errObj = gaaspEE->controlGlobalError(uTOL);
//  double err = errObj->getTotalError();
//  double tol = 1.0e-14;
//  TEST_FLOATING_EQUALITY( err, -0.00293816063746374, tol );
//}

} // namespace Rythmos


