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
#include "AmesosClassType.h"

//! Amesos_Factory:  A method for creating Amesos classes
/*!  Amesos_Factory allows a code to delay the decision about which
concrete class to use to implement the Amesos_BaseSolver interface.  
*/
//
class Amesos_Factory { 

  //@{ \name Creation method
  //! Amesos_Factory Create method
  /*! Creates an instance of the Amesos_BaseSolver class specified by 
    ClassType 
  */
public: 
  Amesos_BaseSolver *Create( AmesosClassType ClassType, 
			     const Epetra_LinearProblem& LinearProblem, 
			     Teuchos::ParameterList &ParameterList );
};  // End of  class Amesos_Factory  
#endif /* _AMESOS_FACTORY_H_ */
