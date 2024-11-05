// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_FUNCTION_OBJECT_H
#define Teuchos_FUNCTION_OBJECT_H

#include "Teuchos_Describable.hpp"


/*! \file Teuchos_FunctionObject.hpp
    \brief An object representation of a function
*/


namespace Teuchos{


/**
 * \brief A function object represents an arbitrary function.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT FunctionObject: public Describable {

public:

  /** \name XMLSerialiation Functions */
  //@{

  /** \brief Returns the string to be used for the value of the
   * type attribute when converting the function to XML.
   */
  virtual std::string getTypeAttributeValue() const=0;

  /** \brief Returns the name of the XML tag used to indicate
   * a funciton object.
   */
  static const std::string& getXMLTagName(){
    static const std::string funcitonTagName = "Function";
    return funcitonTagName;
  }

  //@}

};


} // namespace Teuchos


#endif
