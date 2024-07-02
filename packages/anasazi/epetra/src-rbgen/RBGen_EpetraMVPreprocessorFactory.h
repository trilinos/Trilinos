// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
