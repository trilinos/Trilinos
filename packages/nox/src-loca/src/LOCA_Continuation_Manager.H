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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_MANAGER_H
#define LOCA_CONTINUATION_MANAGER_H

#include "NOX_Abstract_Group.H" // For ReturnType

namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
    class ExtendedGroup;
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
      <li> "%Householder Arc Length" (LOCA::Continuation::HouseholderGroup) 
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
      virtual LOCA::Continuation::ExtendedGroup* 
      createContinuationGroup(LOCA::Continuation::AbstractGroup& grp);

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
