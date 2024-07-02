// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_MS_PREPROCESSOR_H
#define RBGEN_MS_PREPROCESSOR_H

#include "RBGen_Preprocessor.hpp"
#include "RBGen_FileIOHandler.hpp"
#include "RBGen_ConfigDefs.h"

// Forward declaration of Epetra_MultiVector class
class Epetra_MultiVector;

namespace RBGen {
  
  //! Specialization for Preprocessor using Epetra_MultiVector.
  class MSPreprocessor : public Preprocessor< Epetra_MultiVector > {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    MSPreprocessor();

    //! Destructor.
    virtual ~MSPreprocessor() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize preprocessor
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
                     const Teuchos::RCP< FileIOHandler <Epetra_MultiVector> >& fileio = Teuchos::null );

    void Reset() { isInitialized_ = false; }
    //@}

    //! @name Preprocess Methods
    //@{

    //! Preprocess the snapshot set passed in
    void Preprocess( Teuchos::RCP<Epetra_MultiVector>& ss );
    //@}

    //! @name Return Methods
    //@{

    //! Return the multivector used to modify the snapshot set
    Teuchos::RCP< Epetra_MultiVector > getMSVector() const { return msVector_; }

    //! @name Status Methods
    //@{

    //! Return initialized status of the preprocessor
    bool isInitialized() const { return isInitialized_; }

    //@}
  private:

    //! Initialization flag.
    bool isInitialized_;

    //! Preprocessing type
    std::string preprocType_;

    //! Input filename.
    std::string input_file_;

    //! Scalar Scaling.
    double scale_;

    //! Scaling vector for the snapshots.
    std::vector< double > scalings_;
    
    //! Scaling indices for the snapshots.
    std::vector< std::pair<int,int> > scaling_idx_;
	
    //! Available preprocessing types
    std::vector< std::string > preproc_types_;

    //! Pointer to the multivector used to modify the snapshot set
    Teuchos::RCP< Epetra_MultiVector > msVector_;

    //! Pointer to the File I/O Handler object.
    Teuchos::RCP< FileIOHandler< Epetra_MultiVector > > fileio_;
  };
  
} // end of RBGen namespace

#endif // RBGEN_MS_PREPROCESSOR_H
