//@HEADER
/*
************************************************************************

Epetra: Linear Algebra Services Package 
Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef TPETRA_MSGTYPE_H
#define TPETRA_MSGTYPE_H

namespace Tpetra {
	/*! \file Tpetra_MsgType.hpp
	    \brief Tpetra::MsgType Enum
	*/

	/*! \enum MsgType
	    Enumerated list of available message types recognized by Tpetra output object.
	*/
	enum MsgType {	
		Signature = 0x1,     /*!< Signature: brief, unique identification information */
		Summary = 0x2,	     /*!< Summary information */
		Watch = 0x4, 	     /*!< Print notification anytime the object is modified */
	};
  
}

#endif // TPETRA_MSGTYPE_H
