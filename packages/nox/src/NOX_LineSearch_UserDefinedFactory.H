// $Id$ 
// $Source$ 

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

#ifndef NOX_LINESEARCH_USERDEFINED_FACTORY_H
#define NOX_LINESEARCH_USERDEFINED_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Common.H"

namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
}

namespace NOX {
namespace LineSearch {

class Generic;

/*!  \brief Pure virtual interface for users to supply their own line search objects.
  
*/  
class UserDefinedFactory {
    
public:
  
  //! Constructor
  UserDefinedFactory() {}
  
  //! Destructor
  virtual ~UserDefinedFactory() {}
  
  /*! \brief Builds a user defined line search object.
    
    @param gd A global data pointer that contains the top level
    parameter list.  Without storing this inside the line
    searchobject, there is no guarantee that the second parameter \c
    params will still exist.  It can be deleted by the top level RCP.
    
    @param params Sublist with line search construction parameters.
      
  */
  virtual Teuchos::RCP<NOX::LineSearch::Generic> 
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
		  Teuchos::ParameterList& params) const = 0;
  
};
} // namespace LineSearch
} // namespace NOX
#endif
