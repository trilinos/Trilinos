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
