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

#ifndef RBGEN_FILE_IO_HANDLER_HPP
#define RBGEN_FILE_IO_HANDLER_HPP

#include "RBGen_ConfigDefs.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  //! Abstract base class for reading datasets from files.
  /*!
   */
  template<class DataSetType>  
  class FileIOHandler {  
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOHandler() {};

    //! Destructor.
    virtual ~FileIOHandler() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize file reader using 
    virtual void Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params ) = 0;
    
    void Reset() {};
    //@}
    
    //! @name File Reading Methods
    //@{
    
    //! Method for reading multiple files and putting them into an data set.
    virtual Teuchos::RCP< DataSetType > Read( const std::vector<std::string>& filenames ) = 0;

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one data set into a file.
    virtual void Write( const Teuchos::RCP< const DataSetType >& MV, const std::string& filename ) = 0;

    //@}

    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    virtual bool isInitialized() const = 0;

    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_FILE_IO_HANDLER_HPP
