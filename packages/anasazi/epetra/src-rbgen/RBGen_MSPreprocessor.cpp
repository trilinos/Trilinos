// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "RBGen_MSPreprocessor.h"
#include "Epetra_MultiVector.h"

namespace RBGen {

  MSPreprocessor::MSPreprocessor(): isInitialized_(false), scale_(1.0)
  {
    // Insert the acceptable preprocessing types for this class
    // Subtract off vectors given through an input file.
    preproc_types_.push_back("Input File");
    // Subtract off the arithmetic mean of the snapshot set.
    preproc_types_.push_back("Arithmetic Mean");
  }

  void MSPreprocessor::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
                                   const Teuchos::RCP< FileIOHandler<Epetra_MultiVector> >& fileio )
  {
    // Save the file i/o handler.
    fileio_ = fileio;

    // Get the preprocessing method sublist.
    const Teuchos::ParameterList& preproc_params = params->sublist( "Preprocessing Method" );

    // Get the preprocessing method.
    std::string preprocMethod = Teuchos::getParameter< std::string >( preproc_params, "Method" );
    // if (preprocMethod != "Modified Snapshot") THROW AN ERROR!!!    

    // Get the preprocessing type.
    preprocType_ = Teuchos::getParameter< std::string >( preproc_params, "Type" );

    // Check that this type is valid.
    bool typeValid = false;
    std::vector<std::string>::iterator p = preproc_types_.begin();
    while ( p != preproc_types_.end() && !typeValid ) {
      if ( *p == preprocType_ )
        typeValid = true;
      ++p;
    }
    // if (!typeValid) THROW AN ERROR!!!

    if (preprocType_ == "Input File") {
      //
      // Try to get the input file from the parameter list
      //
      try {
        input_file_ = Teuchos::getParameter< std::string >( preproc_params, "Filename" );
      }
      catch (std::exception &e) {
        std::cout<<"The modified snapshot input file has not been specified"<<std::endl;
      }
    }      
    //
    // Check if there is a scaling factor on the steady state vector
    //
    scale_ = 1.0;
    if ( preproc_params.isParameter("Scaling Factor") ) {
      scale_ = Teuchos::getParameter< double >( preproc_params, "Scaling Factor" );
    }
    
    /*    try {
      scalings_ = Teuchos::getParameter< std::vector< double > >( *params, "Snapshot Scaling" );
      }
      catch (std::exception &e) {
      std::cout<<"The snapshot scaling vector has not been specified in the input!"<<std::endl;
      }
      //
      // Try to get the snapshot scaling indices vector
      //
      try {
      scaling_idx_ = Teuchos::getParameter< std::vector< std::pair<int,int> > >( *params, "Snapshot Scaling Indices" );
      }
      catch (std::exception &e) {
      std::cout<<"The snapshot scaling indices have not been specified in the input!"<<std::endl;
    }      
    */
    if ( scalings_.size() != scaling_idx_.size() ) {
      std::cout<<"The scaling vector is not the same size as the number of index pairs!"<< std::endl;
    }
    //
    // We're initialized now.
    // 
    isInitialized_ = true;
  }
  
  void MSPreprocessor::Preprocess( Teuchos::RCP<Epetra_MultiVector>& ss )
  {
    if (isInitialized_) {      
  
      if (preprocType_ == "Input File") {
        //
        // Get the steady state vector from the file.
        //
        std::vector< std::string > steady_file( 1, input_file_ );
        msVector_ = fileio_->Read( steady_file );
        Epetra_MultiVector* colMV = 0;
        // 
        // Remove the scaled steady state vector.
        //
        for (int i=0; i<ss->NumVectors(); i++) {
          colMV = new Epetra_MultiVector( View, *ss, i, 1 );
          colMV->Update( -1.0*scale_, *msVector_, 1.0 );
          delete colMV;
        }
      }
      else if (preprocType_ == "Arithmetic Mean") {
        //
        // Compute the arithmetic mean from the input snapshots
        //
        msVector_ = Teuchos::rcp( new Epetra_MultiVector( ss->Map(), 1 ) );
        Epetra_MultiVector* colMV = 0;
        // 
        // Sum up all the multivectors; compute the average.
        //
        for (int i=0; i<ss->NumVectors(); i++) {
          colMV = new Epetra_MultiVector( View, *ss, i, 1 );
          msVector_->Update( 1.0, *colMV, 1.0 );
          delete colMV;
        }
        msVector_->Scale( 1.0 / ss->NumVectors() );
        // 
        // Remove the arithmetic mean from the snapshot set.
        //
        for (int i=0; i<ss->NumVectors(); i++) {
          colMV = new Epetra_MultiVector( View, *ss, i, 1 );
          colMV->Update( -1.0*scale_, *msVector_, 1.0 );
          delete colMV;
        }      
      }
    }
    else {
      // Throw error here!
    }
  }
  
} // end of RBGen namespace

