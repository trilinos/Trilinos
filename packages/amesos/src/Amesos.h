/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
*/

#ifndef _AMESOS_FACTORY_H_
#define _AMESOS_FACTORY_H_
#include "Amesos_BaseSolver.h"

//! Amesos:  A method for binding a third party direct solver to an Epetra_LinearProblem.
/*!  Amesos creates an instance of a solver, binding a third party direct 
solver to an Epetra_LinearProblem, allowing access to the specified third 
party solver through the Amesos interface (i.e. Numeric Factorization 
SymbolicFactrozation(), Solve() and support functions.)
*/
//
class Amesos { 

  //@{ \name Creation method
  //! Amesos Create method
  /*! Creates an instance of the Amesos_BaseSolver class specified by 
    ClassType.

      <br \>Preconditions:<ul>
      <li>ClassType must be one of the recognized class types.
Return 0 on failure.
      </li>
      <li>ClassType must specify a third party solver that has been 
linked with this particular implementation.  Return 0 on failure.
      </li>
      <li>Epetra_LinearProblem may be empty.  Although the linear problem 
is not checked at the time of construction, the operator must be an Epetra_RowMatrix, or derived from an Epetra_RowMatrix.

      </li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>If Create() returns a non-null pointer, that pointer points to an 
Amesos solver. 
      </li>
      </ul>
  */
public: 
  Amesos_BaseSolver *Create( char *ClassType, 
			     const Epetra_LinearProblem& LinearProblem );

  // @}
  
  // @{ Query methods
  /*! Queries whether a given interface is avaiable or not. */

  bool Query( char * ClassType );
  
};  // End of  class Amesos  
#endif /* _AMESOS_FACTORY_H_ */
