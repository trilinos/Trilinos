// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_FUNCTIONOBJECTXMLCONVERTER_HPP
#define TEUCHOS_FUNCTIONOBJECTXMLCONVERTER_HPP

/*! \file Teuchos_FunctionObjectXMLConverter.hpp
 * \brief Converts back and forth between XML
 * and FunctionObjects.
*/

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_FunctionObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


namespace Teuchos {


/** \brief An abstract base class for converting FunctionObjects to
 * and from XML.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT FunctionObjectXMLConverter : public Describable {

public:

  /** \name Converter Functions */
  //@{

  /** \brief Converts a given XMLObject to a FunctionObject.
   *
   * @param xmlObj The XMLObject to convert to a FunctionObject.
   * @return The converted FunctionObject.
   */
  RCP<FunctionObject> fromXMLtoFunctionObject(const XMLObject& xmlObj) const;

  /** \brief Preforms any and all special xml conversion that
   * is specific to a
   * particular FunctionObject.
   *
   * @param xmlObj The xml to be converted.
   * @return The converted FunctionObject.
   */
  virtual RCP<FunctionObject> convertXML(const XMLObject& xmlObj) const=0;

  /** \brief Converters a given FunctionObject to XML.
   *
   * @param function The FunctionObject to be converted to XML.
   * @return An XML representation of the given FunctionObject.
   */
  XMLObject fromFunctionObjecttoXML(
    const RCP<const FunctionObject> function) const;

  /** \brief Preforms any and all special FunctionObject conversion that is
   * specific to a particlar FunctionObject.
   *
   * @param function The FunctionObject to be converted.
   * @param xmlObj The XMLObject to store all serialization in.
   */
  virtual void convertFunctionObject(
    const RCP<const FunctionObject> function,
    XMLObject& xmlObj) const = 0;

  //@}

  //! \name Attribute/Query Functions
  //@{

  /** \brief . */
  static const std::string& getTypeAttributeName(){
    static const std::string typeAttributeName = "type";
    return typeAttributeName;
  }

  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_FUNCTIONOBJECTXMLCONVERTER_HPP
