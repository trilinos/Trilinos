// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_MANAGER_H
#define LOCA_CONTINUATION_MANAGER_H

#include "NOX_Abstract_Group.H" // For ReturnType

namespace LOCA {
  namespace Abstract {
    class Group;
  }
  namespace Continuation {
    class Group;
  }
}

namespace LOCA {

  namespace Continuation {

    /*! 
      \brief A manager class for creating different continuation groups.
    */
    /*!
      The parameters passed to the constructor or reset should specify
      the "Method", as described below, as well as any additional
      parameters for that particular method.
 
      <ul>
      <li> "Method" - Name of the continuation method. Valid choices are
      <ul> 
      <li> "%Natural" (LOCA::Continuation::NaturalGroup) [Default]
      <li> "%Arc Length" (LOCA::Continuation::ArcLengthGroup) 
      </ul>
      </ul>

      Call the createContinuationGroup to instantiate a continuation group
      of the appropriate type.
    */  

    class Manager {

    public:

      //! Constructor. 
      Manager(NOX::Parameter::List& p);

      //! Destructor
      virtual ~Manager();

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& p);
  
      //! Creates continuation group for given method
      virtual LOCA::Continuation::Group* 
      createContinuationGroup(const LOCA::Abstract::Group& grp, 
			      NOX::Parameter::List& linSolverParams);

      //! Returns method name
      const string& getMethod() const;

      //! Returns continuation parameter id (string)
      const string& getConParamID() const;

    protected:

      //! Stores continuation method name
      string method;
      
      //! Stores continuation parameter id (string)
      string conParamID;

      //! Stores parameter list
      NOX::Parameter::List* paramsPtr;
    };
  } // namespace Predictor
} // namespace LOCA

#endif
