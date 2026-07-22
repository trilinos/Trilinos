// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT Condition : public Describable{
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
