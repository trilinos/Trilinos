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


#ifndef TEUCHOS_EXCEPTIONS_HPP
#define TEUCHOS_EXCEPTIONS_HPP


#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


/** \brief Base exception class for Teuchos
 *
 * \ingroup teuchos_mem_mng_grp
 */
class ExceptionBase : public std::logic_error
{public:ExceptionBase(const std::string& what_arg) : std::logic_error(what_arg) {}};
// 2007/11/07: rabartl: Above, I had to change the name from Exception to
// ExceptionBase because Marzio did a 'using namespace Teuchos' and then he
// declared his own Exception class.  The file Laplacian3D.cpp failed to
// compile.  STOP DOING USING NAMESPACE BLAH!!!!!!


/** \brief Thrown if a duplicate owning RCP is creatd the the same object.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class DuplicateOwningRCPError : public ExceptionBase
{public:DuplicateOwningRCPError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Null reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class NullReferenceError : public ExceptionBase
{public:NullReferenceError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Null reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class NonconstAccessError : public ExceptionBase
{public:NonconstAccessError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Range error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class RangeError : public ExceptionBase
{public:RangeError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Dangling reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class DanglingReferenceError : public ExceptionBase
{public:DanglingReferenceError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Incompatiable iterators error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class IncompatibleIteratorsError : public ExceptionBase
{public:IncompatibleIteratorsError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


} // end namespace Teuchos


#endif	// TEUCHOS_EXCEPTIONS_HPP
