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

// ////////////////////////////////////////////////////////////////////////
// Teuchos_RawMPITraits.hpp

#ifndef TEUCHOS_RAW_MPI_TRAITS_H
#define TEUCHOS_RAW_MPI_TRAITS_H

#include "Teuchos_ConfigDefs.hpp"

/** \file Teuchos_RawMPITraits.hpp
 *  \brief Declaration of a templated traits class that returns raw MPI data types.
 */

namespace Teuchos {

///
/** \brief Templated class that returns raw MPI data types.
 *
 * \note 
 * <ul> <li> This class should not compile if it is instantiated by accident.  
 * 	<li> It should only be included if MPI is available and the MPI header 
 *		<b>must</b> be included before this header file.
 *	<li> Template specializations exist for datatypes: <tt>char</tt>, <tt>int</tt>,
 *		<tt>float</tt>, and <tt>double</tt>.
 * </ul>
 */
template <class T> class RawMPITraits {
public:
	/** \brief Return the raw MPI data type of the template argument */
	static MPI_Datatype type() { bool *junk1; T *junk2 = &junk1; return MPI_DATATYPE_NULL; }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
///
/** Specialization of <tt>RawMPITraits</tt> for <tt>char</tt>
 */
template <> class RawMPITraits<char> {
public:
	///
	static MPI_Datatype type() { return MPI_CHAR; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>int</tt>
 */
template <> class RawMPITraits<int> {
public:
	///
	static MPI_Datatype type() { return MPI_INT; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>float</tt>
 */
template <> class RawMPITraits<float> {
public:
	///
	static MPI_Datatype type() { return MPI_FLOAT; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>double</tt>
 */
template <> class RawMPITraits<double> {
public:
	///
	static MPI_Datatype type() { return MPI_DOUBLE; }
};

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Teuchos

#endif // TEUCHOS_RAW_MPI_TRAITS_H
