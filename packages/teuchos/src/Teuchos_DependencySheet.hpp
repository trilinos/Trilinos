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



#ifndef TEUCHOS_DEPENDENCYSHEET_HPP_
#define TEUCHOS_DEPENDENCYSHEET_HPP_

/*! \file Teuchos_DependencySheet.hpp
    \brief DataStructure keeping track of dependencies.
*/

#include "Teuchos_Dependency.hpp"

namespace Teuchos{

/**
 * \brief A Dependency sheet keeps track of dependencies between various elements located
 * somewhere within a "Root List". 
 *
 * All dependencies added to a DependencySheet
 * must have dependents and dependees who are either in the Root List or one of
 * its sublists.
 * Note that a DependencySheet never acts on these dependencies. It mearly keeps
 * track of them.
 */
class DependencySheet{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convience typedef representing a set of dependencies.
   */
  typedef std::set<RCP<Dependency>, Dependency::DepComp > DepSet;

  /**
   * \brief Convience typedef. Maps dependee parameter entries to a set of their corresponding
   * dependencies.
   */
  typedef map<const ParameterEntry*, DepSet > DepMap;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs an empty DependencySheet with the name DEP_ANONYMOUS.
   * 
   * @param rootList The Parameter List containing all parameters and sublists
   * for which this Dependency will keep track of dependencies.
   */
  DependencySheet(RCP<ParameterList> rootList);

  /**
   * \brief Constructs a DependencySheet.
   *
   * @param rootList The Parameter List containing all parameters and sublists
   * for which this Dependency will keep track of dependencies.
   * @param name Name of the Dependency Sheet.
   */
  DependencySheet(RCP<ParameterList> rootList, const std::string &name);
  
  //@}

  /** \name Add/Remove Functions */
  //@{

  /**
   * \brief Adds a dependency to the sheet.
   * 
   * @param dependency The dependency to be added.
   * @return True if the addition was sucessful, false otherwise.
   */
  bool addDependency(const RCP<Dependency> dependency);

  /**
   * \brief Removes a particular dependency between two parameters.
   *
   * @param dependency The dependency to be removed.
   * @return True if the removal was sucessfull, false otherwise.
   */
  bool removeDependency(RCP<Dependency> dependency);
  
  //@}

  //! \name Attribute/Query Functions
  //@{
  
  /**
   * \brief Determines whether or not a parameter is depended upon by any another
   * parameters or parameter lists.
   *
   * @param name The paramteter to be checked for dependents.
   * @return True if the parameter you're checking has other dependents, false otherwise.
   */
  bool hasDependents(const ParameterEntry *dependee) const;

  /**
   * \brief Returns a set of all the dependencies associated with a particular dependee.
   *
   * @param dependee The parameter whose dependencies are sought. 
   * @return A set of all dependencies associated with the dependee parameter.
   * */
  const DepSet& getDependenciesForParameter(const ParameterEntry *dependee) const;
  
  //@}

  /** \name Iterator Functions */
  //@{

  /**
   * \brief Returns an iterator to the beginning of all the dependees in the sheet.
   */
  DepMap::iterator depBegin();

  /**
   * Returns an iterator to the end of all of the dependees in the sheet.
   *
   * @return An iterator to the end of all of the dependees in the sheet.
   */
  DepMap::iterator depEnd();

  /**
   * \brief Returns a const iterator to the beginning of all the dependees in the sheet.
   */
  DepMap::const_iterator depBegin() const;

  /**
   * \brief Returns a const iterator to the end of all of the dependees in the sheet.
   */
  DepMap::const_iterator depEnd() const;
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief Prints out a list of the dependencies in the DependencySheet
   */
  void printDeps();

  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map containing all the depenecies for a list. 
   *
   * Dependencies with multiple dependees will be found in multiple
   * places within the map. Essentially, for each dependee, there will 
   * be a pointer to for dependency of which it is a part.
   */
  DepMap dependencies_;

  /**
   * \brief The Name of the dependency sheet.
   */
  std::string name_;

  /**
   * \brief The root parameterlist that this dependency sheet is associated with.
   */
  RCP<ParameterList> rootList_;

  /**
   * \brief Validates whether or not the dependees and dependents of a dependency exist
   * within the root ParameterList.
   */
  void validateExistanceInRoot(RCP<Dependency> dependency);
  
  //@}

};


}
#endif //TEUCHOS_DEPENDENCYSHEET_HPP_
