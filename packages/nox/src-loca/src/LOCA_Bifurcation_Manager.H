// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_BIFURCATION_MANAGER_H
#define LOCA_BIFURCATION_MANAGER_H

#include "NOX_Abstract_Group.H" // For ReturnType

namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
  }
}

namespace LOCA {

  namespace Bifurcation {

    /*! 
      \brief A manager class for creating different bifurcation groups.
    */
    /*!
      The parameters passed to the constructor or reset should specify
      the "Bifurcation Method" , as described below, as well as any additional
      parameters for that particular method.  In particular these parameters
      should specify the initial starting data for bifurcation.
 
      <ul>
      <li> "Method" - Name of the bifurcation method. Valid choices are
      <ul> 
      <li> "None" - No bifurcation [Default]
      <li> "Turning Point" (LOCA::Bifurcation::TPBord::ExtendedGroup)
      <li> "Modified Turning Point" (LOCA::Bifurcation::TPBord::ModifiedBorderingGroup)
      <li> "Pitchfork" (LOCA::Bifurcation::PitchforkBord::ExtendedGroup) 
      <li> "Hopf" (LOCA::Bifurcation::HopfBord::ExtendedGroup)
      </ul>
      </ul>

      Call the createBifurcationGroup to instantiate a bifurcation group
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
  
      //! Creates bifurcation group for given method
      virtual LOCA::Continuation::AbstractGroup* 
      createBifurcationGroup(LOCA::Continuation::AbstractGroup& grp);

      //! Returns method name
      const string& getMethod() const;

    protected:

      //! Stores continuation method name
      string method;

      //! Parameters
      NOX::Parameter::List* paramsPtr;
    };
  } // namespace Predictor
} // namespace LOCA

#endif
