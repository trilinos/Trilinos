// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_OUTPUT_MANAGER_HPP
#define TPETRA_OUTPUT_MANAGER_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MsgType.hpp"

namespace Tpetra {

//! Tpetra::OutputManager: Primary class for managing Tpetra object output content.
/*!  Tpetra's basic output manager for sending information of select verbosity levels
  to the appropriate output stream.

  This output manager will remove the need for the derived class to know any information
  about the required output.  Calling <tt>isVerbosity( MsgType type )</tt> will inform the solver if
  it is supposed to output the information corresponding to the message type (\c type ).
  
  \author Mike Heroux
*/

class OutputManager {

public:

  //@{ \name Constructors/Destructor.

  //! Default constructor
  OutputManager();

  //! Basic constructor.
  OutputManager( int myID, int vb = Tpetra::Signature, int printID = 0, ostream& os = std::cout );

  //! Destructor.
  virtual ~OutputManager() {}
  //@}
	
  //@{ \name Set methods.

  //! Set the output stream for this manager.
  //void setOstream ( const ostream & os ) { myOS_ = os; }

  //! Set the message output types for this manager.
  void setVerbosity( int vb ) { vb_ = vb; }

  //@}

  //@{ \name Get methods.

  //! Get the output stream for this manager.
  ostream& getOstream() { return myOS_; }

  //@}

  //@{ \name Query methods.

  //! Find out whether we need to print out information for this message type.
  /*! This method is used by the solver to determine whether computations are
    necessary for this message type.
  */
  bool isVerbosity( MsgType type ) { return (( type == Tpetra::Signature ) || ( vb_ & type )); }

  //! Find out whether this processor needs to print out information for this message type.
  /*! This method is used by the solver to determine whether this output stream has been
    selected to output the information for this message type.
  */		
  bool isVerbosityAndPrint( MsgType type ) { return ( iPrint_ && isVerbosity( type ));}

  //! Find out whether information can be outputted through this output stream.
  bool doPrint( void ) { return (iPrint_);}	

  //@}

private:

  //@{ \name Undefined methods.

  //! Copy constructor.
  OutputManager( const OutputManager & om );

  //! Assignment operator.
  OutputManager& operator=( const OutputManager & om );

  //@}

  int myID_, printID_;
  int vb_;
  bool iPrint_;
  ostream& myOS_;	
};

OutputManager::OutputManager() :
  myID_(0),
  printID_(0),
  vb_(Tpetra::Signature),
  iPrint_(true),
  myOS_(std::cout)
{
}

OutputManager::OutputManager( int myID, int vb, int printID, ostream& os ) :
  myID_(myID),
  printID_(printID),
  vb_(vb),
  iPrint_(myID == printID),
  myOS_(os)
{
}

} // end Tpetra namespace

#endif

// end of file TpetraOutputManager.hpp
