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


/** \brief An xml converter for SingleOperatorFunctions */
template<class OperandType>
class SingleOperatorFunctionXMLConverter : 
  public FunctionObjectXMLConverter<OperandType, OperandType>
{

public:

  /** \name Special Converter Functions */
  //@{

  /**
   * \brief Gets the specific SingleOperatorFunction to be returned
   *  by this converter when converting from XML.
   *
   *  @parameter operand The operand to be used with the SingleOperatorFunction.
   */
  virtual RCP<SingleOperatorFunction<OperandType> > 
    getSpecificSingleOperatorFunction(OperandType operand) const = 0;

  /**
   * \brief Add and extra XML traits that are specific to a certain
   * FuncitonOjbect when converting that function object to xml.
   *
   * @param functionObject The function object being convertered.
   * @param xmlObj The XMLObject to which any special traits should be added.
   */
  virtual void getSpecialSingleOperatorFunctionXMLTraits(
    const RCP<const SingleOperatorFunction<OperandType> > functionObject,
    XMLObject& xmlObj) const{}

  //@}

  /** \name Overridden from FunctionObjectXMLConverter */
  //@{

  /** \brief . */
  RCP<FunctionObject<OperandType, OperandType> > 
    convertXML(const XMLObject& xmlObj) const;

  /** \brief . */
  void convertFunction(
    const RCP<const FunctionObject<OperandType, OperandType> > functionObject, 
    XMLObject& xmlObj) const;
  
  //@}
 
private:

  //! \name Private mebers.
  //@{

  static std::string getOperandAttributeName(){
    static const std::string operandAttributeName = "operand";
    return operandAttributeName;
  }

  //@}
  
};

RCP<FunctionObject<OperandType, OperandType> > 
  SingleOperatorFunctionXMLConverter::convertXML(const XMLObject& xmlObj) const
{
  OperandType operand = 
    xmlObj.getRequired<OperandType>(getOperandAttributeName());
  return getSpecificSingleOperatorFunction(operand);
}

void convertFunction(
  const RCP<const FunctionObject<OperandType, OperandType> > functionObject, 
  XMLObject& xmlObj) const
{
  RCP<const SingleOperatorFunctions<OperandType> > castedFunction =
    rcp_dynamic_cast<const SingleOperatorFunction<OperandType> >(
      functionObject, true); 
  OperandType operand = castedFunction->getModifiyingOperand();
  xmlObj.addAttribute(getOperandAttributeName(),operand);
  getSpecialSingleOperatorFunctionXMLTraits(xmlObj, castedFunction);
}

template<class OperandType>
class SubtractionFunctionXMLConverter :
  public SingleOperatorFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SingleOperatorFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SingleOperatorFunction<OperandType> > 
    getSpecificSingleOperatorFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SingleOperatorFunction<OperandType> >
  getSpecificSingleOperatorFunction(OperandType operand) const
{
  return rcp(new SubtractionFunction<OperandType>(operand));
}

template<class OperandType>
class AdditionFunctionXMLConverter :
  public SingleOperatorFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SingleOperatorFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SingleOperatorFunction<OperandType> > 
    getSpecificSingleOperatorFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SingleOperatorFunction<OperandType> >
  getSpecificSingleOperatorFunction(OperandType operand) const
{
  return rcp(new AdditionFunction<OperandType>(operand));
}

template<class OperandType>
class MultiplicationFunctionXMLConverter :
  public SingleOperatorFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SingleOperatorFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SingleOperatorFunction<OperandType> > 
    getSpecificSingleOperatorFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SingleOperatorFunction<OperandType> >
  getSpecificSingleOperatorFunction(OperandType operand) const
{
  return rcp(new MultiplicationFunction<OperandType>(operand));
}

template<class OperandType>
class DivisionFunctionXMLConverter :
  public SingleOperatorFunctionXMLConverter<OperandType>
{
public:

  /** \name Overridden from SingleOperatorFunctionXMLConverter */
  //@{
  
  /** \brief. */
  RCP<SingleOperatorFunction<OperandType> > 
    getSpecificSingleOperatorFunction(OperandType operand) const;

  //@}
};

template<class OperandType>
RCP<SingleOperatorFunction<OperandType> >
  getSpecificSingleOperatorFunction(OperandType operand) const
{
  return rcp(new DivisionFunction<OperandType>(operand));
}


} // namespace Teuchos


#endif // TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP

