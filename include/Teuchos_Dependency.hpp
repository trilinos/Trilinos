// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DEPENDENCY_HPP_
#define TEUCHOS_DEPENDENCY_HPP_
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_InvalidDependencyException.hpp"
#include <set>
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
