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

#ifndef LOCA_PREDICTOR_RANDOM_H
#define LOCA_PREDICTOR_RANDOM_H

#include "LOCA_Predictor_Generic.H"  // base class

namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace LOCA {

  namespace Predictor {

    //! %Random predictor computation class
    /*! 
      This class computes the predictor direction where the solution 
      component is filled with random values and the parameter component
      equal to 1.  Each componenet of the solution vector \f$v_i\f$ of the 
      predictor is given by \f$v_i = \epsilon r_i x_i\f$ where 
      \f$r_i\f$ is a random value between 0 and 1, \f$x_i\f$ is the
      corresponding componenet of the solution vector, and \f$\epsilon\f$
      is a parameter.

      The parameters used by this class supplied in the constructor or reset
      method are:
      <ul>
      <li> "Epsilon" - \f$\epsilon\f$ as defined above (Default 1.0e-3)
      </ul>
    */
    class Random : public LOCA::Predictor::Generic {

    public:

      //! Constructor. 
      Random(NOX::Parameter::List& params);

      //! Destructor
      virtual ~Random();

      //! Reset parameters.
      /*!
	The only parameter for the random predictor is \f$epsilon\f$.
      */
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      //! Computes the random predictor.
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::Group& prevGroup,
	      LOCA::Continuation::Group& curGroup,
	      LOCA::Continuation::Vector& result);

    protected:
      double epsilon;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
