
//@HEADER
// ***********************************************************************
// 
//                     Rythmos Package
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
//@HEADER

#ifndef _NEWP_HELLO_H_
#define _NEWP_HELLO_H_
#include "Rythmos_ConfigDefs.h"
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
