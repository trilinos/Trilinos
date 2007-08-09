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

#ifndef NOX_DIRECTION_FACTORY_H
#define NOX_DIRECTION_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Common.H"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}

namespace NOX {

class GlobalData;

namespace Direction {

class Generic;

/*! 

  \brief Factory to build direction objects derived from NOX::Direction::Generic.  

  <b>Parameters</b>
  
   <ul>
   <li> "Method" <std::string> Name of the direction. Valid choices are:
     <ul> 
     <li> "Newton" NOX::Direction::Newton
     <li> "Steepest Descent" NOX::Direction::SteepestDescent
     <li> "NonlinearCG" NOX::Direction::NonlinearCG
     <li> "Broyden" NOX::Direction::Broyden 
     <li> "Tensor" PRERELEASE ONLY! NOX::Direction::Tensor
     <li> "Modified-Newton" PRERELEASE ONLY! NOX::Direction::ModifiedNewton
     <li> "Quasi-Newton" PRERELEASE ONLY! NOX::Direction::QuasiNewton
     <li> "User Defined" - see below
     </ul>

   <li> "User Defined Constructor" - see below

   <li> "Newton" <sublist> Parameters to build a NOX::Direction::Newton object.
   <li> "Steepest Descent" <sublist> Parameters to build a NOX::Direction::SteepestDescent object.
   <li> "NonlinearCG" <sublist> Parameters to build a NOX::Direction::NonlinearCG object.
   <li> "Broyden" <sublist> Parameters to build a NOX::Direction::Broyden object.
   <li> "Tensor" <sublist> Parameters to build a NOX::Direction::Tensor object.
   <li> "Modified-Newton" <sublist> Parameters to build a NOX::Direction::ModifiedNewton object.
   <li> "Quasi-Newton" <sublist> Parameters to build a NOX::Direction::QuasiNewton object.
   <li> "User Defined Direction Factory" < RCP<NOX::Direction::UserDefinedFactory> > RCP to a NOX::Direction::UserDefinedFactory derived object.  This factory object is used to build user defined direction objects.
   </ul>

<b>Using a User-Defined Direction</b>

The user has the option of passing in a user-defined direction.  First, they must implement their own direction, deriving from the base class interface NOX::Direction::Generic:

\code
   class MyDirection : public NOX::Direction::Generic {
   
     // Ctor that takes the standard direction arguments.
     MyDirection(const Teuchos::RCP<NOX::GlobalData>& gd, 
                 Teuchos::ParameterList& params);

     .
     .
     .

   };
\endcode

Next they must write a factory to build their object, deriving from the NOX::Direction::UserDefinedFactory base class interface:

\code
class MyFactory {
  MyDirFactory();
  ~MyDirFactory();
   Teuchos::RCP<NOX::Direction::Generic> buildDirection(const Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
   .
   .
   .
};
\endcode

Then under the "Direction" parameter sublist, they need to set the
method to "User Defined" and register the factory:

\code
using namespace Teuchos; // for RCP and ParameterList
ParameterList& dir_params = p.sublist("Direction");
RCP<NOX::Direction::UserDefinedFactory> dir_facotry = rcp(new MyDirectionFactory);
dir_params.set("Method", "User Defined");
dir_params.set("User Defined Direction Factory", dir_factory);
\endcode

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::Direction::UserDefinedFactory.

 */  
class Factory {

public:

  //! Constructor
  Factory();

  //! Destructor
  ~Factory();

  /*! \brief Factory to build a direction object.

      @param gd A global data pointer that contains the top level parameter list.  Without storing this inside the direction object, there is no guarantee that the second parameter \c params will still exist.  It can be deleted by the top level RCP.
      @param params Sublist with direction construction parameters.

  */
  Teuchos::RCP<NOX::Direction::Generic> 
  buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
		 Teuchos::ParameterList& params);
  
};

  /*! Nonmember function to build a direction object.

  \relates NOX::Direction::Factory

  */
  Teuchos::RCP<NOX::Direction::Generic> 
  buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
		 Teuchos::ParameterList& params);
  

} // namespace Direction
} // namespace NOX
#endif
