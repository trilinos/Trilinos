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
