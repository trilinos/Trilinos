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



#ifndef TEUCHOS_CONDITION_HPP_
#define TEUCHOS_CONDITION_HPP_

/*! \file Teuchos_Condition.hpp
    \brief An object to determin if a particular set of conditions
    are occuring.
*/

#include "Teuchos_Dependency.hpp"
#include "Teuchos_Describable.hpp"

namespace Teuchos {

/**
 * \brief A Condition determines whether or not 
 * a particular set of conditions are currently
 * occuring.
 */
class TEUCHOS_LIB_DLL_EXPORT Condition : public Describable{
public:

  /** \name Public types */
  //@{
  
  /**
   * \brief Convenience typedef for defining a list of conditions.
   */
  typedef Teuchos::Array<Teuchos::RCP<Condition> > ConditionList;

  typedef Teuchos::Array<Teuchos::RCP<const Condition> > ConstConditionList;

    //@}

  /** \name Constructors/Destructor */
  //@{
  
  /**
   * \brief Constructs a Condition
   */
  Condition(){}

  /**
   * \brief Destructs a condition.
   */
  virtual ~Condition(){}

  //@}


  /** Attribute/Query Functions */
  //@{
  
  /** \brief Determins whether or not a condition is true. */
  virtual bool isConditionTrue() const = 0;

  /**
   * \brief Determines whether or not the evaluation of a parameter
   * occurs somewhere in this condition.
   *
   * @return Whether or not the evaluation of a parameter
   * occurs somewhere in this condition.
   */
  virtual bool containsAtLeasteOneParameter() const = 0;

  /**
   * \brief Gets all of the parameters that are evaluated in this
   * condition.
   */
  virtual Dependency::ConstParameterEntryList getAllParameters() const = 0;

  /** \brief Get the value that should be used for the condition type
   * attribute when converting a condition to XML.
   */
  virtual std::string getTypeAttributeValue() const = 0;

  /** \brief . */
  static const std::string& getXMLTagName(){
    static const std::string xmlTagName = "Condition";
    return xmlTagName;
  }

  //@}
  
};

}
#endif //TEUCHOS_CONDITION_HPP_
