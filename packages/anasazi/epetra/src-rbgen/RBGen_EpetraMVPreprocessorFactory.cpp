// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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

