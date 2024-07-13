// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_EpetraMVPreprocessorFactory.h"

#include "Epetra_MultiVector.h"

#include "Teuchos_Assert.hpp"

namespace RBGen {

  EpetraMVPreprocessorFactory::EpetraMVPreprocessorFactory()
  {
    preproc_types.push_back("None");
    preproc_types.push_back("Modified Snapshot");
  }

  Teuchos::RCP<Preprocessor<Epetra_MultiVector> >
  EpetraMVPreprocessorFactory::create( const Teuchos::ParameterList& params )
  {
    // See if the "Preprocessing" sublist exists
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist( "Preprocessing Method" ), std::invalid_argument, "Preprocessing Method sublist does not exist!");

    // Get the preprocessing method sublist.
    const Teuchos::ParameterList& preproc_params = params.sublist( "Preprocessing Method" );

    // Get the preprocessing parameter.
    std::string method = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(preproc_params),
                                                             "Method" );

    // Create the preprocessor.
    Teuchos::RCP<Preprocessor< Epetra_MultiVector > > RBPreprocessor;

    // Basic preprocessor does nothing
    if ( method == "None" ) {
      RBPreprocessor = Teuchos::rcp( new NoPreprocessor<Epetra_MultiVector>() );
    } else
    // Modified snapshot preprocessor
    if ( method == "Modified Snapshot" ) {
      RBPreprocessor = Teuchos::rcp( new MSPreprocessor() );
    } else
    {
       // Throw an exception because the preprocessing method was not recognized by this factory.
    }
    //
    // Return the preprocessor created
    //
    return RBPreprocessor;

  }

} // end of RBGen namespace

