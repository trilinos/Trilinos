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

#ifndef LOCA_PREDICTOR_TANGENT_H
#define LOCA_PREDICTOR_TANGENT_H

#include "LOCA_Predictor_Generic.H"  // base class
#include "NOX_Parameter_List.H"      // class data element

namespace LOCA {

  namespace Predictor {

    //! %Tangent predictor computation class
    /*! 
      This class computes the predictor direction given by the tangent
      vector to the continuation curve.  The tangent vector calculation
      is implemented by the computeTangent method
      of the LOCA::Continuation::Group class.
    */
    class Tangent : public LOCA::Predictor::Generic {

    public:

      //! Constructor. 
      Tangent(NOX::Parameter::List& params);

      //! Destructor
      virtual ~Tangent();

      //! Reset parameters.
      /*!
	The parameters are eventually passed to the 
	LOCA::Continuation::Group::computeTangent method and therefore should
	be the parameters for the linear solver.
      */
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      //! Computes the tangent predictor by calling curGroup.computeTangent().
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::Group& prevGroup,
	      LOCA::Continuation::Group& curGroup,
	      LOCA::Continuation::Vector& result);

    protected:
      //! Stores parameter list for computeTangent call.
      NOX::Parameter::List parameters;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
