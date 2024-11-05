// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
class SimpleFunctionXMLConverter : public FunctionObjectXMLConverter{

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
    const RCP<const SimpleFunctionObject<OperandType> > /* functionObject */,
    XMLObject& /* xmlObj */) const{}

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
class SubtractionFunctionXMLConverter :
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
class AdditionFunctionXMLConverter :
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
class MultiplicationFunctionXMLConverter :
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
class DivisionFunctionXMLConverter :
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

