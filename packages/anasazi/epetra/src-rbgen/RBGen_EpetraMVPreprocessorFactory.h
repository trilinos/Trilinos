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
