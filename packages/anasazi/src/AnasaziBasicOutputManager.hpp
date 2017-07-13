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

#ifndef ANASAZI_BASIC_OUTPUT_MANAGER_HPP
#define ANASAZI_BASIC_OUTPUT_MANAGER_HPP

/*!     \file AnasaziBasicOutputManager.hpp
        \brief Basic output manager for sending information of select verbosity levels to the appropriate output stream
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_oblackholestream.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*!  \class Anasazi::BasicOutputManager

  \brief Anasazi's basic output manager for sending information of select verbosity levels
  to the appropriate output stream.

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

namespace Anasazi {

  using std::ostream;

  template <class ScalarType>
  class BasicOutputManager : public OutputManager<ScalarType> {

    public:

      //! @name Constructors/Destructor
      //@{ 

      //! Default constructor
      BasicOutputManager(int vb = Anasazi::Errors, 
                         Teuchos::RCP<ostream> os = Teuchos::rcpFromRef(std::cout),
                         int printingRank = 0);

      //! Destructor.
      virtual ~BasicOutputManager() {};
      //@}

      //! @name Set/Get methods
      //@{ 

      //! Set the output stream for this manager.
      void setOStream( Teuchos::RCP<ostream> os );

      //! Get the output stream for this manager.
      Teuchos::RCP<ostream> getOStream();

      //@}

      //! @name Output methods
      //@{ 

      //! Find out whether we need to print out information for this message type.
      /*! This method is used by the solver to determine whether computations are
        necessary for this message type.
        */
      bool isVerbosity( MsgType type ) const;

      //! Send some output to this output stream.
      void print( MsgType type, const std::string output );

      //! Return a stream for outputting to.
      ostream &stream( MsgType type );

      //@}

    private:

      //! @name Undefined methods
      //@{ 

      //! Copy constructor.
      BasicOutputManager( const OutputManager<ScalarType>& OM );

      //! Assignment operator.
      BasicOutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );

      //@}

      Teuchos::RCP<ostream> myOS_;
      Teuchos::oblackholestream myBHS_;
      bool iPrint_;
  };

  template<class ScalarType>
  BasicOutputManager<ScalarType>::BasicOutputManager(int vb, Teuchos::RCP<ostream> os, int printingRank)
  : OutputManager<ScalarType>(vb), myOS_(os) {
    // print only on proc 0
    int MyPID;
#ifdef HAVE_MPI
    // Initialize MPI
    int mpiStarted = 0;
    MPI_Initialized(&mpiStarted);
    if (mpiStarted) MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
    else MyPID=0;
#else 
    MyPID = 0;
#endif
    iPrint_ = (MyPID == printingRank);
  }

  template<class ScalarType>
  void BasicOutputManager<ScalarType>::setOStream( Teuchos::RCP<ostream> os ) { 
    myOS_ = os; 
  }

  template<class ScalarType>
  Teuchos::RCP<ostream> BasicOutputManager<ScalarType>::getOStream() { 
    return myOS_; 
  }

  template<class ScalarType>
  bool BasicOutputManager<ScalarType>::isVerbosity( MsgType type ) const {
    if ( (type & this->vb_) == type ) {
      return true;
    }
    return false;
  }

  template<class ScalarType>
  void BasicOutputManager<ScalarType>::print( MsgType type, const std::string output ) {
    if ( (type & this->vb_) == type && iPrint_ ) {
      *myOS_ << output;
    }
  }

  template<class ScalarType>
  ostream & BasicOutputManager<ScalarType>::stream( MsgType type ) {
    if ( (type & this->vb_) == type && iPrint_ ) {
      return *myOS_;
    }
    return myBHS_;
  }

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
