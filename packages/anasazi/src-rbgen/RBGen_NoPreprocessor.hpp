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

#ifndef RBGEN_NO_PREPROCESSOR_HPP
#define RBGEN_NO_PREPROCESSOR_HPP

#include "RBGen_Preprocessor.hpp"

namespace RBGen {

  //! Null preprocessor concrete class.
  template< class DataSetType > 
  class NoPreprocessor : public Preprocessor< DataSetType > {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    NoPreprocessor() {};

    //! Destructor.
    virtual ~NoPreprocessor() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize preprocessor
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
                     const Teuchos::RCP< FileIOHandler<DataSetType> >& fileio ) {};

    //! Reset preprocessor
    void Reset() {};
    //@}

    //! @name Preprocess Methods
    //@{

    //! Preprocess the snapshot set passed in
    void Preprocess( Teuchos::RCP<DataSetType>& ss ) {};
    //@}

    //! @name Status Methods
    //@{

    //! Return initialized status of the preprocessor
    bool isInitialized() const { return true; };
    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_NO_PREPROCESSOR_HPP
