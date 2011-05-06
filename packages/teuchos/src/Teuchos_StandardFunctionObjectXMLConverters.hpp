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

#ifndef TEUCHOS_STANDARDFUNCTIONOBJECTXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDFUNCTIONOBJECTXMLCONVERTERS_HPP

/*! \file Teuchos_StandardFunctionObjectXMLConverters.hpp
 * \brief A collection of standard FunctionObjectXMLConverters.
*/


#include "Teuchos_FunctionObjectXMLConverter.hpp"
#include "Teuchos_StandardFunctionObjects.hpp"


namespace Teuchos {


/** \brief An xml converter for SimpleFunctionObjects */
template<class OperandType>
class TEUCHOS_LIB_DLL_EXPORT SimpleFunctionXMLConverter : public FunctionObjectXMLConverter{

public:

  /** \name Special Converter Functions */
  //@{

  /**
   * \brief Gets the specific SimpleFunction to be returned
   *  by this converter when converting from XML.
   *
   *  @parameter operand The operand to be used with the SimpleFunction.
   */
  virtual RCP<SimpleFunctionObject<OperandType> > 
    getSpecificSimpleFunction(OperandType operand) const = 0;

  /**
   * \brief Add and extra XML traits that are specific to a certain
   * FuncitonOjbect when converting that function object to xml.
   *
   * @param functionObject The function object being convertered.
   * @param xmlObj The XMLObject to which any special traits should be added.
   */
  virtual void getSpecialSimpleFunctionXMLTraits(
    const RCP<const SimpleFunctionObject<OperandType> > functionObject,
    XMLObject& xmlObj) const{}

  //@}

  /** \name Overridden from FunctionObjectXMLConverter */
  //@{

  /** \brief . */
  RCP<FunctionObject> 
    convertXML(const XMLObject& xmlObj) const;

  /** \brief . */
  void convertFunctionObject(
    const RCP<const FunctionObject> functionObject, 
    XMLObject& xmlObj) const;
  
  //@}

  //! \name Constant Retrieving Functions
  //@{

  static std::string getOperandAttributeName(){
    static const std::string operandAttributeName = "operand";
    return operandAttributeName;
  }

  //@}
  
};

template<class OperandType>
RCP<FunctionObject> 
SimpleFunctionXMLConverter<OperandType>::convertXML(
  const XMLObject& xmlObj) const
{
  OperandType operand = 
    xmlObj.getRequired<OperandType>(getOperandAttributeName());
  return getSpecificSimpleFunction(operand);
}

template<class OperandType>
void SimpleFunctionXMLConverter<OperandType>::convertFunctionObject(
  const RCP<const FunctionObject> functionObject, 
  XMLObject& xmlObj) const
{
  RCP<const SimpleFunctionObject<OperandType> > castedFunction =
    rcp_dynamic_cast<const SimpleFunctionObject<OperandType> >(
      functionObject, true); 
  OperandType operand = castedFunction->getModifiyingOperand();
  xmlObj.addAttribute(getOperandAttributeName(),operand);
  getSpecialSimpleFunctionXMLTraits(castedFunction, xmlObj);
}

/**
 * \brief Class for converting SubtractionFunction objects to and from XML.
 *
 * The valid XML represntation of a SubtractionFunction is:
 * \code
   <Function operand="operand value" type="SubtractionFunction(operand_type)"/>
 \endcode
 */
template<class OperandType>
class TEUCHOS_LIB_DLL_EXPORT SubtractionFunctionXMLConverter :
  public SimpleFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SimpleFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SimpleFunctionObject<OperandType> > 
    getSpecificSimpleFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SimpleFunctionObject<OperandType> >
SubtractionFunctionXMLConverter<OperandType>::getSpecificSimpleFunction(
  OperandType operand) const
{
  return rcp(new SubtractionFunction<OperandType>(operand));
}

/**
 * \brief Class for converting AdditionFunction objects to and from XML.
 *
 * The valid XML represntation of a AdditionFunction is:
 * \code
   <Function operand="operand value" type="AdditionFunction(operand_type)"/>
 \endcode
 */
template<class OperandType>
class TEUCHOS_LIB_DLL_EXPORT AdditionFunctionXMLConverter :
  public SimpleFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SimpleFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SimpleFunctionObject<OperandType> > 
    getSpecificSimpleFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SimpleFunctionObject<OperandType> >
AdditionFunctionXMLConverter<OperandType>::getSpecificSimpleFunction(
  OperandType operand) const
{
  return rcp(new AdditionFunction<OperandType>(operand));
}

/**
 * \brief Class for converting MultiplicationFunction objects to and from XML.
 *
 * The valid XML represntation of a MultiplicationFunction is:
 * \code
   <Function operand="operand value" type="MultiplicationFunction(operand_type)"/>
 \endcode
 */
template<class OperandType>
class TEUCHOS_LIB_DLL_EXPORT MultiplicationFunctionXMLConverter :
  public SimpleFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SimpleFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SimpleFunctionObject<OperandType> > 
    getSpecificSimpleFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SimpleFunctionObject<OperandType> >
MultiplicationFunctionXMLConverter<OperandType>::getSpecificSimpleFunction(
  OperandType operand) const
{
  return rcp(new MultiplicationFunction<OperandType>(operand));
}

/**
 * \brief Class for converting DivisionFunction objects to and from XML.
 *
 * The valid XML represntation of a DivisionFunction is:
 * \code
   <Function operand="operand value" type="DivisionFunction(operand_type)"/>
 \endcode
 */
template<class OperandType>
class TEUCHOS_LIB_DLL_EXPORT DivisionFunctionXMLConverter :
  public SimpleFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SimpleFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SimpleFunctionObject<OperandType> > 
    getSpecificSimpleFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SimpleFunctionObject<OperandType> >
DivisionFunctionXMLConverter<OperandType>::getSpecificSimpleFunction(
  OperandType operand) const
{
  return rcp(new DivisionFunction<OperandType>(operand));
}


} // namespace Teuchos


#endif // TEUCHOS_STANDARDFUNCTIONOBJECTXMLCONVERTERS_HPP

