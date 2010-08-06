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
class Dependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Allows two dependecies to be compared.
   */
  class DepComp{
  public:
    bool operator () (const RCP<Dependency> dep1,
      const RCP<Dependency> dep2) const
    {
      return dep1->getDependents().begin().get()
        >= dep2->getDependents().begin().get();
    }
  };

  /**
   * \brief A list of Dependees.
   */
  typedef std::set<RCP<ParameterEntry> > DependeeList;

  /**
   * \brief A list of dependents.
   */
  typedef std::set<RCP<ParameterEntry> > DependentList;

  /**
   * \brief A list of Dependees.
   */
  typedef std::set<RCP<const ParameterEntry> > ConstDependeeList;

  /**
   * \brief A list of dependents.
   */
  typedef std::set<RCP<const ParameterEntry> > ConstDependentList;

  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a Dependency
   *
   * @param dependees A list of all the dependees.
   * @param dependents A list of all the dependents.
   */
  Dependency(DependeeList& dependees, DependentList& dependents);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependees A list of all the dependees.
   * @param dependent The dependent parameter.
   */
  Dependency(DependeeList& dependees, RCP<ParameterEntry> dependent);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependee The dependee parameter.
   * @param dependents A List of all the dependents.
   */
  Dependency(RCP<ParameterEntry> dependeeName, DependentList& dependents);

  /**
   * \brief Constructs a Dependency
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   */
  Dependency(RCP<ParameterEntry> dependee, RCP<ParameterEntry> dependent);

  /**
   * \brief Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~Dependency(){}
  
  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * \brief Gets the dependees of the dependency.
   *
   *  @return The dependees of the dependency.
   */
  inline const ConstDependeeList& getDependees() const{
    return dependees_;
  }

  /**
   * \brief Gets the dependents of the dependency.
   *
   * @return The dependents of the dependency.
   */
  inline const ConstDependentList& getDependents() const{
    return dependents_;
  }

  /** \brief Gets the first dependee in the dependees list.
   * This is a convience function.
   */
  inline RCP<const ParameterEntry> getFirstDependee() const{
    return dependees_.begin();
  }

  /**
   * \brief Convienence function. 
   * Returns the first dependee in the list of dependees.
   *
   * @return The first dependee in the list of dependees.
   */
  template<class S>
  inline const S getFirstDependeeValue() const{
    return getValue<S>(dependees_.begin());
  }

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
  DependeeList dependees_;

  /**
   * \brief The dependent paramters.
   */
  DependentList dependents_;

  /**
   * \brief Declaring and defining the default constructor as private.
   */
  Dependency(){}

  //@}

};


} //namespace Teuchos
#endif //TEUCHOS_DEPENDENCY_HPP_
