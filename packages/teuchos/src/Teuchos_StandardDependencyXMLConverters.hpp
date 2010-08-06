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

#ifndef TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardDependencyXMLConverters.hpp
 * \brief A collection of standard DependencyXMLConverters.
*/


#include "Teuchos_DependencyXMLConverter.hpp"


namespace Teuchos {


/** \brief An xml converter for VisualDepenencies
 */
class VisualDependencyConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific visual dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @return An XMLObject representing the VisualDepenency.
   */
  virtual void convertSpecialVisualAttributes(
    RPC<const VisualDepenency> dependency,
    XMLObject& xmlObj) const = 0;

  /** \brief Converts any special aspects of a
   * specific visual dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependees The dependees of the visual dependency.
   * @param dependents The dependents of the visua dependency.
   * @param showIf The showIf attribute of the visual dependency.
   * @return The converted VisualDependency.
   */
  RCP<VisualDepenency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const ParameterParentMap& dependees,
    const ParameterParentMap& dependets,
    bool showIf) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const ParameterParentMap& dependees, 
    const ParameterParentMap& dependents) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj) const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  static const std::string& getShowIfAttributeName(){
    static const std::string showIfAttributeName = "showIf";
    return showIfAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for StringVisualDepenencies
 */
class StringVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDepenency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDepenency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const ParameterParentMap& dependees,
    const ParameterParentMap& dependets,
    bool showIf) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{
  
  inline
  const std::string& getTypeAttributeValue() const{
    return TypeNameTraits<StringVisualDependency>::name();
  }
  
  //@}


private:

  /** \name Private Members */
  //@{
  
  static const std::string& getStringValuesTagName(){
    static const std::string stringValuesTagName = "StringValues";
    return stringValuesTagName;
  }

  static const std::string& getStringTagName(){
    static const std::string stringTagName = "string";
    return stringTagName;
  }

  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for StringVisualDepenencies
 */
class StringVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDepenency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDepenency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const ParameterParentMap& dependees,
    const ParameterParentMap& dependets,
    bool showIf) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{
  
  inline
  const std::string& getTypeAttributeValue() const{
    return TypeNameTraits<BoolVisualDependency>::name();
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  static const std::string& getStringValuesTagName(){
    static const std::string stringValuesTagName = "StringValues";
    return stringValuesTagName;
  }

  static const std::string& getStringTagName(){
    static const std::string stringTagName = "string";
    return stringTagName;
  }

  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for BoolVisualDepenencies
 */
class BoolVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  XMLObject convertSpecialVisualAttributes(
    RCP<const VisualDepenency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDepenency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const ParameterParentMap& dependees,
    const ParameterParentMap& dependets,
    bool showIf) const;
  
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

