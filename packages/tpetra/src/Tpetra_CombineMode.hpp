// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef _TPETRA_COMBINEMODE_HPP_
#define _TPETRA_COMBINEMODE_HPP_

/*! \file Tpetra_CombineMode.hpp 
    \brief Tpetra::Combine Mode enumerable type
*/

namespace Tpetra {
	
	/*! \enum CombineMode
      If set to Add, existing values will be summed with new values.
      If set to Replace, existing values will be replaced with new values.
	*/
	
	enum CombineMode {
		Add, /*!< Existing values will be summed with new values. */
		Replace /*!< Existing values will be replaced with new values. */
	};

} // namespace Tpetra

#endif /* _TPETRA_COMBINEMODE_HPP_ */
