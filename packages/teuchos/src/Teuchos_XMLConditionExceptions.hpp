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


#ifndef TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_
#define TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_

/*! \file Teuchos_XMLConditionExceptions.hpp
 * \brief A collection of Exceptions thrown
 * when converting Conditions to and from
 * XML.
 */
#include <stdexcept>

namespace Teuchos {

/** \brief Thrown when a StringConditon is missing it's Value tag.
 */
class MissingValuesTagException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingValuesTagException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValuesTagException(const std::string& what_arg):
    std::logic_error(what_arg){}

};

/** \brief Thrown when an appropriate Condition Converter can't be found */
class CantFindConditionConverterException : public std::logic_error{

public:

  /**
   * \brief Constructs an CantFindConditionConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindConditionConverterException(const std::string& what_arg):
    std::logic_error(what_arg){}

};



} // namespace Teuchos
#endif //TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_

