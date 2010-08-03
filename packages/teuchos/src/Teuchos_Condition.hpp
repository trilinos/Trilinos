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



#ifndef TEUCHOS_CONDITION_HPP_
#define TEUCHOS_CONDITION_HPP_

#include "Teuchos_Dependency.hpp"

namespace Teuchos{

/**
 * A Condition determines whether or not 
 * a particular set of conditions are currently
 * occuring.
 */
class Condition{
public:
	/**
	 * Conviencne typedef for defining a list of conditions.
	 */
	typedef Teuchos::Array<Teuchos::RCP<Condition> > ConditionList;

	/**
	 * Enum to determine the conditoin type.
	 */
	enum Type{BinLogicCon, NotCon, ParamCon};

	/**
	 * Construcst a Condition
	 */
	Condition(Type type):type_(type){}

	/**
	 * Deconstructs a condition.
	 */
	virtual ~Condition(){}

	/**
	 * Determins whether or not a condition is true.
	 */
	virtual bool isConditionTrue() const = 0;

	/**
	 * Determines whether or not the evaluation of a parameter
	 * occurs somewhere in this condition.
	 *
	 * @return Whether or not the evaluation of a parameter
	 * occurs somewhere in this condition.
	 */
	virtual bool containsAtLeasteOneParameter() const = 0;

	/**
	 * Gets all of the parameters that are evaluated in this
	 * condition.
	 *
	 * @return A map of all of the parameters that are evaluated in this
	 * condition.
	 */

	virtual Dependency::ParameterParentMap getAllParameters() const = 0;

	/**
	 * Gets the type of the condition.
	 *
	 * @return The type of the condition.
	 */
	const Type& getType() const{
		return type_;
	}

private:
	/**
	 * The Conditions's type.
	 */
	Type type_;

};

}
#endif //TEUCHOS_CONDITION_HPP_
