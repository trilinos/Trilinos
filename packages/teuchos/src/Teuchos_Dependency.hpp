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



#ifndef TEUCHOS_DEPENDENCY_HPP_
#define TEUCHOS_DEPENDENCY_HPP_
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_InvalidDependencyException.hpp"
namespace Teuchos{


/**
 * This class represents a depndency between elements in a Parameter List.
 * 
 * @see DependencyList 
 */
class Dependency{
public:
  /**
   * Allows two dependecies to be compared.
   */
  class DepComp{
  public:
    bool operator () (const RCP<Dependency> dep1, const RCP<Dependency> dep2) const{
      return dep1->getDependents().begin()->first >= dep2->getDependents().begin()->first;
    }
  };

  /**
   * Enum classifying various types of dependencies.
   */
  enum Type{VisualDep, ValidatorDep, NumberValidatorAspectDep, NumberArrayLengthDep};

  /**
   * Maps parameters to their associated parent ParametersList.
   */
  typedef std::map<const std::string, const RCP<ParameterList> > ParameterParentMap;

  /**
   * Constructs a Dependency
   *
   * @param dependees A map of all the dependees and their associated parent lists.
   * @param dependents A map of all the dependents and their associated parent lists.
   * @param depType The type of dependency.
   */
  Dependency(ParameterParentMap& dependees, ParameterParentMap& dependents, Type depType);

  /**
   * Constructs a Dependency
   *
   * @param dependees A map of all the dependees and their associated parent lists.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param depType The type of dependency.
   */
  Dependency(ParameterParentMap& dependees, std::string dependentName, RCP<ParameterList> dependentParentList, Type depType);

  /**
   * Constructs a Dependency
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map of all the dependents and their associated parent lists.
   * @param depType The type of dependency.
   */
  Dependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap& dependents, Type depType);

  /**
   * Constructs a Dependency
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param depType The type of dependency.
   */
  Dependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, Type depType);

  /**
   * Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~Dependency(){}

  /**
   * Determines whether or not a ParameterList or any of it's children lists contain a specific
   * ParameterList.
   *
   * @param parentList The ParameterList to search.
   * @param listToFind The ParameterList to for which we are searching.
   * @return True if the parentList or and or any of it's children ParameterLists contains the list
   * specified by the listToFind parameter.
   */
  static bool doesListContainList(RCP<ParameterList> parentList, RCP<ParameterList> listToFind);

  /**
   * Gets the dependees of the dependency.
   *
   *  @return The dependees of the dependency.
   */
  const ParameterParentMap& getDependees() const;

  /**
   * Gets the dependents of the dependency.
   *
   * @return The dependents of the dependency.
   */
  const ParameterParentMap& getDependents() const;

  /**
   * Gets the names of the dependee parameters.
   *
   * @return The names of the dependee parameters.
   */
  const std::set<std::string>& getDependeeNames() const;

  /**
   * Gets a set containing the names of the dependent parameters.
   *
   * @return A set containing the names of the dependent parameters.
   */
  const std::set<std::string>& getDependentNames() const;

  /**
   * Gets a string containing all the names of the dependee parameters.
   *
   * @return A string containing all the names of the dependee parameters.
   */
  std::string getDependeeNamesString() const;

  /**
   * Gets a string containing all the names of the dependent parameters.
   *
   * @return A string containing all the names of the dependent parameters.
   */
  std::string getDependentNamesString() const;

  /**
   * Gets the name of a dependee given a pointer to the dependee parameter.
   *
   * @param Pointer to the dependee parameter whose name is desired.
   * @return The name of the dependee parameter associated with the pointer specified in the arguments.
   */
  std::string getDependeeName(const ParameterEntry* dependee) const;

  /**
   * Gets the type of the dependency.
   *
   * @return The type of dependency.
   */
  const Type& getType() const;

  /**
   * Evaluates the dependency and makes any appropriate changes to the
   * dependee based on the dependent.
   */
  virtual void evaluate() = 0;

  /**
   * Convienence function. Returns the first dependee in the list of dependees.
   *
   * @return The first dependee in the list of dependees.
   */
  inline const ParameterEntry* getFirstDependee() const{
    return dependees_.begin()->second->getEntryPtr(dependees_.begin()->first);
  }

  /**
   * Convienence function. Returns the first dependee in the list of dependees.
   *
   * @return The first dependee in the list of dependees.
   */
  template<class S>
  inline const S getFirstDependeeValue() const{
    return dependees_.begin()->second->get<S>(dependees_.begin()->first);
  }

  /**
   * Gets the name of the first dependee in the dependees map.
   *
   * @return the name of the first dependee in the dependees map.
   */
  inline const std::string& getFirstDependeeName() const{
    return dependees_.begin()->first;
  }

private:

  /**
   * The dependee is the parameter being depended upon.
   * This is a map of all the dependees and their associated parent ParameterLists.
   */
  ParameterParentMap dependees_;

  /**
   * The dependent is the parameter that dependes on another parameter.
   * This is a map of all the dependents and their associated parent ParametersLists.
   */
  ParameterParentMap dependents_;

  /**
   * The names of all the dependees
   */
  std::set<std::string> dependeeNames_;

  /**
   * The names of all the dependents
   */
  std::set<std::string> dependentNames_;

  /**
   * The type of dependency.
   */
  Type type_;

  /**
   * Initializes all the dependnees and dependents along with checking to make sure
   * that their parents lists are actually valid.
   *
   * @param dependees The dependees to be initialized.
   * @param dependents The dependents to be initialized.
   */
  void intitializeDependeesAndDependents(ParameterParentMap& dependees, ParameterParentMap& dependents);

  /**
   * Validates the dependency to make sure it's valid/has been setup properly. If subclassing, this fucntion should
   * be called in the new subclasses constructor.
   */
  virtual void validateDep() const = 0;

};

}
#endif //TEUCHOS_DEPENDENCY_HPP_
