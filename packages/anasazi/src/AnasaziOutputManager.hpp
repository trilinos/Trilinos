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

#ifndef ANASAZI_OUTPUT_MANAGER_HPP
#define ANASAZI_OUTPUT_MANAGER_HPP

/*!     \file AnasaziOutputManager.hpp
        \brief Abstract class definition for Anasazi Output Managers.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

/*!  \class Anasazi::OutputManager

  \brief Output managers remove the need for the eigensolver to know any information
  about the required output.  Calling isVerbosity( MsgType type ) informs the solver if
  it is supposed to output the information corresponding to the message type.

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

namespace Anasazi {

template <class ScalarType>
class OutputManager {

  public:

    //!@name Constructors/Destructor 
  //@{ 

  //! Default constructor
  OutputManager( int vb = Anasazi::Errors ) : vb_(vb) {};

  //! Destructor.
  virtual ~OutputManager() {};
  //@}
  
  //! @name Set/Get methods
  //@{ 

  //! Set the message output types for this manager.
  virtual void setVerbosity( int vb ) { vb_ = vb; }

  //! Get the message output types for this manager.
  virtual int getVerbosity( ) const { return vb_; }

  //@}

  //! @name Output methods
  //@{ 

  //! Find out whether we need to print out information for this message type.
  /*! This method is used by the solver to determine whether computations are
      necessary for this message type.
  */
  virtual bool isVerbosity( MsgType type ) const = 0;

  //! Send output to the output manager.
  virtual void print( MsgType type, const std::string output ) = 0;

  //! Create a stream for outputting to.
  virtual std::ostream &stream( MsgType type ) = 0;

  //@}

  private:

  //! @name Undefined methods
  //@{ 

  //! Copy constructor.
  OutputManager( const OutputManager<ScalarType>& OM );

  //! Assignment operator.
  OutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );

  //@}

  protected:
  int vb_;
};

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
