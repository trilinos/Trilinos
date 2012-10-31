// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

