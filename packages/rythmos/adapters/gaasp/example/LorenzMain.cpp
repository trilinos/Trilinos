//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER 


#include "Rythmos_GAASPErrorEstimator.hpp"
#include "LorenzModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

int main(int argc, char *argv[])
{

  using std::endl;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::CommandLineProcessor;
  
  bool result, success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = rcp( new Epetra_SerialComm );
#endif // HAVE_MPI

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read commandline options
    //

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    bool dumpFinalSolutions = false;
    clp.setOption(
      "dump-final-solutions", "no-dump-final-solutions", &dumpFinalSolutions,
      "Determine if the final solutions are dumpped or not." );
    
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    if ( Teuchos::VERB_DEFAULT == verbLevel )
      verbLevel = Teuchos::VERB_LOW;

    const Teuchos::EVerbosityLevel
      solnVerbLevel = ( dumpFinalSolutions ? Teuchos::VERB_EXTREME : verbLevel );

    //
    // Get the base parameter list that all other parameter lists will be read
    // from.
    //
    
    RCP<ParameterList>
      paramList = Teuchos::parameterList();

    //
    // Create the underlying EpetraExt::ModelEvaluator
    //

    RCP<LorenzModel>
      lorenzModel = rcp(new LorenzModel( epetra_comm, *paramList ));

    //
    // Create the Thyra-wrapped ModelEvaluator
    //
    
    RCP<Thyra::ModelEvaluator<double> >
      thyraLorenzModel = Thyra::epetraModelEvaluator(lorenzModel,Teuchos::null);

    //
    // Create the Rythmos GAASP ErrorEstimator
    //

    RCP<Rythmos::GAASPErrorEstimator> gaaspEE = rcp(new Rythmos::GAASPErrorEstimator);
    gaaspEE->setModel(thyraLorenzModel);
    //gaaspEE->setQuantityOfInterest( AVERAGE_ERROR_QTY );  // Not passed through yet.
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->sublist("GAASP Interface Parameters").set<double>("eTime",1.0);
    pl->sublist("GAASP Interface Parameters").set<double>("timeStep",0.1);
    //pl->sublist("GAASP Interface Parameters").set<double>("timeStep",0.5);
    gaaspEE->setParameterList(pl);

    //RCP<const Rythmos::ErrorEstimateBase<double> > error = gaaspEE->getErrorEstimate();

    //double uTOL = 1.0e-8;
    double uTOL = 1.0e-2;
    RCP<const Rythmos::ErrorEstimateBase<double> > error = gaaspEE->controlGlobalError(uTOL);


    double err = error->getTotalError();
    out->precision(15);
    *out << "err = " << err << std::endl;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

//  if(success)
//    *out << "\nEnd Result: TEST PASSED" << endl;
//  else
//    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

