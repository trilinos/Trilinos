
//@HEADER
// ***********************************************************************
// 
//                     New_Package Example Package
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

#include "New_Package_ConfigDefs.h"
#include "Newp_Hello.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"

/** \brief A sample class.
 *
 * This class prints out a "Hello World" type message.
 *
 * Here is what output from this class looks like for an example program:
 *
 * \verbinclude hello_test.out
 *
 * <b>A typical heading</b>
 * <ul>
 * <li> A typical first list entry
 * <li> A typical second list entry
 * </ul>
 *
 * <b>Another typical heading</b>
 */
class Newp_Hello {
public:
  
  /** \name Constructors/destructors. */
  //@{
  
  /** \brief Creates a Newp_Hello object and fills with default values.
   * 
   * \param  Comm [in] An Epetra Communicator.
   *
   * \warning Newp_Hello is English language only.  In Africa use Newp_Jambo.
   * 
   * \return  Newp_Hello object
   */
  Newp_Hello(const Epetra_Comm& Comm);

  /** \brief . */
  Newp_Hello(const Newp_Hello& Source);

  //@}
  
  /** \name Print functions */
  //@{

  /** \brief . */
  virtual void Print(ostream & os) const;

  //@}

 private:

  const Epetra_Comm Comm_ ; // Must be stored by value since this is a handle!

};

#endif /* _NEWP_HELLO_H_ */
