
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef _NEWP_HELLO_H_
#define _NEWP_HELLO_H_
#include "New_Package_config.h"
#include "Epetra_ConfigDefs.h"
#include "Newp_Hello.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
//! Newp_Hello: A sample class 

/*! The Newp_Hello class prints out a "Hello World" type message

<b>A typical heading</b>
<ul>
  <li> A typical first list entry
  <li> A typical second list entry
</ul>

<b>Another typical heading</b>

*/

//=========================================================================
class Newp_Hello {

  public:

  //@{ \name Constructors/destructors.
  //! Basic Newp_Hello constuctor.
  /*! Creates a Newp_Hello object and fills with default values.  

    \warning Newp_Hello is English language only.  In Africa use Newp_Jambo.

    \param Comm In
           An Epetra Communicator 

    \return  Newp_Hello object

  */
  Newp_Hello(const Epetra_Comm& Comm);

  //! Newp_Hello copy constructor.
  
  Newp_Hello(const Newp_Hello& Source);
  //@}
  
  //@{ \name Print methods

  //! Print method
  virtual void Print(ostream & os) const;
  //@}


 private:

  const Epetra_Comm& Comm_ ; 

};

#endif /* _NEWP_HELLO_H_ */
