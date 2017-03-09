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

#ifndef RBGEN_EPETRAMV_PREPROCESSOR_FACTORY_H
#define RBGEN_EPETRAMV_PREPROCESSOR_FACTORY_H

#include "RBGen_PreprocessorFactory.hpp"
#include "RBGen_NoPreprocessor.hpp"
#include "RBGen_MSPreprocessor.h"
#include "RBGen_ConfigDefs.h"

// Forward declaration of Epetra_MultiVector.
class Epetra_MultiVector;

namespace RBGen {

  //! Specialization of a PreprocessorFactor for Epetra_MultiVector datasets.
  class EpetraMVPreprocessorFactory : public virtual PreprocessorFactory<Epetra_MultiVector> {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    EpetraMVPreprocessorFactory();

    //! Destructor.
    virtual ~EpetraMVPreprocessorFactory() {};
    //@}

    //! @name Factory methods
    //@{

    Teuchos::RCP<Preprocessor<Epetra_MultiVector> > create( const Teuchos::ParameterList& params );

    //@}

  private:

    // Available preprocessing types
    std::vector<std::string> preproc_types;

  };
 
} // end of RBGen namespace

#endif
