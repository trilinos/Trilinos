//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER

#ifndef BELOS_OUTPUT_MANAGER_HPP
#define BELOS_OUTPUT_MANAGER_HPP

/*!     \file BelosOutputManager.hpp
        \brief Class which manages the output and verbosity of the Belos solvers.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_oblackholestream.hpp"	
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*!	\class Belos::OutputManager

	\brief Belos's basic output manager for sending information of select verbosity levels
	to the appropriate output stream.

	This output manager will remove the need for the solver or linear problem to know any information
	about the required output.  Calling <tt>isVerbosity( MsgType vb )</tt> will inform the solver if
	it is supposed to output the information corresponding to the verbosity type (\c vb ).

	\author Michael Heroux and Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType>
  class OutputManager {
    
  public:
    
    //! @name Constructors/Destructor
    //@{ 
    
    //! Basic constructor.
    OutputManager( int vb = Belos::Errors, const Teuchos::RCP< std::ostream > &os = Teuchos::rcp(&std::cout,false) );
    
    //! Destructor.
    virtual ~OutputManager() {};
    //@}
    
    //! @name Set methods
    //@{ 
    
    //! Set the output stream for this manager.
    void setOStream( const Teuchos::RCP<std::ostream> &os ) { myOS_ = os; };
    
    //! Set the verbosity level for this manager.
    void setVerbosity( int vb ) { vb_ = vb; }; 
    
    //@}
    
    //! @name Get methods
    //@{ 

    //! Get an output stream for outputting the input message type.
    std::ostream& stream( MsgType type ) 
    {
      if ( (type & vb_) && iPrint_ ) {
	return *myOS_;
      }
      return myBHS_;
    }
 
    //! Get the output stream for this manager.
    Teuchos::RCP<std::ostream> getOStream() { return myOS_; };
    
    //@}
    
    //! @name Query methods
    //@{ 
    
    //! Find out whether we need to print out information for this message type.
    /*! This method is used by the solver to determine whether computations are
      necessary for this message type.
    */
    bool isVerbosity( MsgType type ) const { return (( type == Belos::Errors ) || ( vb_ & type )); }; 
    
    //@}

    //! @ name Print methods
    //@{
    
    //! Send some output of a specified message type to the output stream.
    void print( MsgType type, const std::string output );

    //@}

  private:
    
    //! @name Undefined methods
    //@{ 
    
    //! Copy constructor.
    OutputManager( const OutputManager<ScalarType>& OM );
    
    //! Assignment operator.
    OutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );
    
    //@}
    
    int vb_;
    Teuchos::RCP<std::ostream> myOS_;	
    Teuchos::oblackholestream myBHS_;  
    bool iPrint_;
  };
  
  template<class ScalarType>
  OutputManager<ScalarType>::OutputManager( int vb, const Teuchos::RCP<std::ostream> &os ) :
    vb_(vb),
    myOS_(os)
  {
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
    iPrint_ = (MyPID == 0);
  }
 
  template<class ScalarType>
  void OutputManager<ScalarType>::print( MsgType type, const std::string output ) {
  if ( (type & vb_) && iPrint_ ) {
    *myOS_ << output;
  }
}
 
} // end Belos namespace

#endif

// end of file BelosOutputManager.hpp
