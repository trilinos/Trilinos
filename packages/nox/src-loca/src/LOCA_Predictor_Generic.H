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

#ifndef LOCA_PREDICTOR_GENERIC_H
#define LOCA_PREDICTOR_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List

// Forward declarations
namespace LOCA {
  namespace Continuation {
    class Vector;
    class Group;
  }
}

namespace LOCA {

  namespace Predictor {

    class Generic {

    public:

      //! Constructor. 
      Generic() {};

      //! Destructor
      virtual ~Generic() {};

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params) = 0;
  
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::Group& prevGroup,
	      LOCA::Continuation::Group& curGroup,
	      LOCA::Continuation::Vector& result) = 0;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
