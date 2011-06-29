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
 * \brief A Dependency sheet keeps track of dependencies between various
 * ParameterEntries
 */
class TEUCHOS_LIB_DLL_EXPORT DependencySheet{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convience typedef representing a set of dependencies.
   */
  typedef std::set<RCP<Dependency>, RCPComp > DepSet;

  /**
   * \brief Convience typedef. Maps dependee parameter entries to a set of their corresponding
   * dependencies.
   */
  typedef map<RCP<const ParameterEntry>, DepSet, RCPConstComp > DepMap;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs an empty DependencySheet with the name DEP_ANONYMOUS.
   */
  DependencySheet();

  /**
   * \brief Constructs a DependencySheet.
   *
   * @param name Name of the Dependency Sheet.
   */
  DependencySheet(const std::string &name);
  
  //@}

  /** \name Add/Remove Functions */
  //@{

  /**
   * \brief Adds a dependency to the sheet.
   * 
   * @param dependency The dependency to be added.
   */
  void addDependency(RCP<Dependency> dependency);

  /**
   * \brief Adds a dependencies from another she 
   * to this sheet.
   * 
   * @param otherSheet The other sheet from which
   * to add dependencies.
   */
  void addDependencies(RCP<DependencySheet> otherSheet);

  /**
   * \brief Removes a particular dependency between two parameters.
   *
   * @param dependency The dependency to be removed.
   * @return True if the removal was sucessfull, false otherwise.
   */
  void removeDependency(RCP<Dependency> dependency);

  /**
   * \brief sets the name of the dependency sheet
   */
  void setName(const std::string newName){
    name_ = newName;
  }
  
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
  inline bool hasDependents(RCP<const ParameterEntry> dependee) const{
    return (dependenciesMap_.find(dependee) != dependenciesMap_.end() 
      && dependenciesMap_.find(dependee)->second.size() > 0);
  }

  /**
   * \brief Returns a set of all the dependencies associated with a 
   * particular dependee. If no dependencies with the given dependee have
   * ever been enetered into the sheet then a null reference is returned.
   *
   * @param dependee The parameter whose dependencies are sought. 
   * @return A set of all dependencies associated with the dependee parameter.
   * */
  RCP<const DepSet> getDependenciesForParameter(
    RCP<const ParameterEntry> dependee) const;

  /**
   * \brief Gets the name of the dependency sheet.
   */
  inline const std::string& getName() const{
    return name_;
  }

  /**
   * \brief Determines whether or not this dependency sheet has any dependencies.
   */
  inline bool empty() const{
    return dependencies_.empty();
  }
  
  //@}

  /** \name Iterator Functions */
  //@{

  /**
   * \brief Returns an iterator to the beginning of all 
   * the dependees in the sheet.
   */
  inline DepSet::iterator depBegin(){
    return dependencies_.begin();
  }

  /**
   * Returns an iterator to the end of all of the dependees in the sheet.
   *
   * @return An iterator to the end of all of the dependees in the sheet.
   */
  inline DepSet::iterator depEnd(){
    return dependencies_.end();
  }

  /**
   * \brief Returns a const iterator to the beginning of all the dependees in the sheet.
   */
  inline DepSet::const_iterator depBegin() const{
    return dependencies_.begin();
  }
    
  /**
   * \brief Returns a const iterator to the end of all of the dependees in the sheet.
   */
  inline DepSet::const_iterator depEnd() const{
    return dependencies_.end();
  }

  /**
   * \brief Returns the number of Dependencies in this
   * DependencySheet.
   *
   * @return The number of Depenedencies in this
   * DependencySheet.
   */
  inline DepSet::size_type size(){
    return dependencies_.size();
  }
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief Prints out a list of the dependencies in the DependencySheet
   */
  void printDeps(std::ostream& out) const;


  /** \brief When serializing to XML, this string should be used as the name
   * of the name attribute */
  static const std::string& getNameAttributeName(){
    static const std::string NAME_ATTRIBUTE_NAME="name";
    return NAME_ATTRIBUTE_NAME;
  }
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map containing all the depenecies for a list. 
   *
   * Dependencies with multiple dependees will be found in multiple
   * places within the map. Essentially, for each dependee, there will 
   * be a pointer to each dependency of which it is a part.
   */
  DepMap dependenciesMap_;

  /**
   *\brief A simple set of all the dependencies in this Dependency Sheet.
   */
  DepSet dependencies_;

  /**
   * \brief The Name of the dependency sheet.
   */
  std::string name_;

  //@}

};


}
#endif //TEUCHOS_DEPENDENCYSHEET_HPP_
