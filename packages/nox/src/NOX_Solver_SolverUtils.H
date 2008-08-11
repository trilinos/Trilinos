//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_SOLVERUTILS_H
#define NOX_SOLVER_SOLVERUTILS_H

#include "Teuchos_RCP.hpp"
#include "NOX_StatusTest_Generic.H"

namespace Teuchos {
  class ParameterList;
}

namespace NOX {
  
  namespace Solver{
    

    /*! \brief Nonmember method that returns the status test check type.

        This object parses the "Solver Options" parameter list for a
        key "Status Test Check Type" of type <std::string> with
        possible values:

	- "Complete"
	- "Minimal" (default)
	- "None"

	These options correspond to the NOX::StatusTest::CheckType.
	Please follow the link for this object for more information.
    */
    NOX::StatusTest::CheckType 
    parseStatusTestCheckType(Teuchos::ParameterList& solver_options_list);
    
  }  // namespace Solver
}  // namespace NOX
#endif
