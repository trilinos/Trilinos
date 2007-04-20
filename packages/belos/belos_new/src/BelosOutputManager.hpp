// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_OUTPUT_MANAGER_HPP
#define BELOS_OUTPUT_MANAGER_HPP

/*!     \file BelosOutputManager.hpp
        \brief Class which manages the output and verbosity of the Belos solvers.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_oblackholestream.hpp"	

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
    OutputManager( int vb = Belos::Errors, const Teuchos::RefCountPtr<ostream> &os = Teuchos::rcp(&std::cout,false) );
    
    //! Destructor.
    virtual ~OutputManager() {};
    //@}
    
    //! @name Set methods
    //@{ 
    
    //! Set the output stream for this manager.
    void setOStream( const Teuchos::RefCountPtr<ostream> &os ) { myOS_ = os; };
    
    //! Set the verbosity level for this manager.
    void setVerbosity( int vb ) { vb_ = vb; }; 
    
    //@}
    
    //! @name Get methods
    //@{ 

    //! Get an output stream for outputting the input message type.
    ostream& stream( MsgType type ) 
    {
      if ( (type & vb_) && iPrint_ ) {
	return *myOS_;
      }
      return myBHS_;
    }
 
    //! Get the output stream for this manager.
    Teuchos::RefCountPtr<ostream> getOStream() { return myOS_; };
    
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
    void print( MsgType type, const string output );

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
    Teuchos::RefCountPtr<ostream> myOS_;	
    Teuchos::oblackholestream myBHS_;  
    bool iPrint_;
  };
  
  template<class ScalarType>
  OutputManager<ScalarType>::OutputManager( int vb, const Teuchos::RefCountPtr<ostream> &os ) :
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
  
} // end Belos namespace

#endif

// end of file BelosOutputManager.hpp
