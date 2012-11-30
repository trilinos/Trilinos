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


#ifndef TEUCHOS_DEPENDENCY_HPP_
#define TEUCHOS_DEPENDENCY_HPP_
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_InvalidDependencyException.hpp"
/*! \file Dependency.hpp
    \brief DataStructure for expressing dependencies between elements
    in ParameterLists.
*/


namespace Teuchos{


/**
 * \brief This class represents a depndency between elements in a Parameter List.
 * 
 * \reference DependencySheet
 * \reference ParameterList
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT Dependency : public Describable{

public:

  /** \name Public types */
  //@{

  /**
   * \brief A list of Dependees.
   */
  typedef std::set<RCP<ParameterEntry>, RCPComp > ParameterEntryList;

  /**
   * \brief A list of dependents.
   */
  typedef std::set<RCP<const ParameterEntry>, RCPConstComp > ConstParameterEntryList;

  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Dependency
   *
   * @param dependees A list of all the dependees.
   * @param dependents A list of all the dependents.
   */
  Dependency(
    ConstParameterEntryList dependees, 
    ParameterEntryList dependents);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependees A list of all the dependees.
   * @param dependent The dependent parameter.
   */
  Dependency(
    ConstParameterEntryList dependees, 
    RCP<ParameterEntry> dependent);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependee The dependee parameter.
   * @param dependents A List of all the dependents.
   */
  Dependency(
    RCP<const ParameterEntry> dependee, 
    ParameterEntryList dependents);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   */
  Dependency(
    RCP<const ParameterEntry> dependee, 
    RCP<ParameterEntry> dependent);

  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * \brief Gets the dependees of the dependency.
   *
   *  @return The dependees of the dependency.
   */
  inline const ConstParameterEntryList& getDependees() const{
    return dependees_;
  }

  /**
   * \brief Gets the dependents of the dependency.
   *
   * @return The dependents of the dependency.
   */
  inline ParameterEntryList& getDependents(){
    return dependents_;
  }

  /**
   * \brief Gets the dependents of the dependency.
   *
   * @return The dependents of the dependency.
   */
  inline const ConstParameterEntryList& getDependents() const{
    return constDependents_;
  }


  /** \brief Gets the first dependee in the dependees list.
   * This is a convience function.
   */
  inline RCP<const ParameterEntry> getFirstDependee() const{
    return *(dependees_.begin());
  }

  /**
   * \brief Convienence function. 
   * Returns the first dependee in the list of dependees.
   *
   * @return The first dependee in the list of dependees.
   */
  template<class S>
  inline S getFirstDependeeValue() const{
    return getValue<S>(*(*(dependees_.begin())));
  }

  /** 
   * \brief Returns the string to be used for the value of the 
   * type attribute when converting the dependency to XML.
   */
  virtual std::string getTypeAttributeValue() const = 0;
 

  /**
   * \brief Returns the XML tag to use when serializing Dependencies.
   */
  static const std::string& getXMLTagName(){
    static const std::string xmlTagName = "Dependency";
    return xmlTagName;
  }

  //@}

  /** \name Evalutation Functions */
  //@{

  /**
   * \brief Evaluates the dependency and makes any appropriate changes to the
   * dependee based on the dependent.
   */
  virtual void evaluate() = 0;
  
  //@}
  
  /** \name I/O Functions */
  //@{

  /** \brief prints out information about the dependency. */
  virtual void print(std::ostream& out) const;

  //@}

protected:

  /** \name Validation Functions */
  //@{

  /**
   * \brief Validates the dependency to make sure it's valid/has been setup 
   * properly. If subclassing, this fucntion should
   * be called in the new subclasses constructor.
   */
  virtual void validateDep() const = 0;

  //@}

 private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The parameters being depended upon.
   */
  ConstParameterEntryList dependees_;

  /**
   * \brief The dependent paramters.
   */
  ParameterEntryList dependents_;

  /**
   * \brief A const version dependent paramters.
   */
  ConstParameterEntryList constDependents_;

  /**
   * \brief Declaring and defining the default constructor as private.
   */
  Dependency(){}

  /**
   * \brief creates a const version of the dependent parameters.
   */
  void createConstDependents();

  /**
   * \brief makes sure the none of the dependess and dependets are null.
   */
  void checkDependeesAndDependents();

  //@}

};


} //namespace Teuchos
#endif //TEUCHOS_DEPENDENCY_HPP_
