// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

