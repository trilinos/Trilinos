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

#ifndef LOCA_PREDICTOR_MANAGER_H
#define LOCA_PREDICTOR_MANAGER_H

#include "LOCA_Predictor_Generic.H"  // base class

namespace LOCA {

  namespace Predictor {

    //! %Manager for all predictor direction computations
    /*! 
      The parameters passed to the constructor or reset should specify
      the "Method", as described below, as well as any additional
      parameters for that particular method.
 
      <ul>
      <li> "Method" - Name of the predictor direction method. Valid choices are
      <ul> 
      <li> "%Constant" (LOCA::Predictor::Constant) [Default]
      <li> "%Tangent" (LOCA::Predictor::Tangent) 
      <li> "%Secant" (LOCA::Predictor::Secant) 
      <li> "%Random" (LOCA::Predictor::Random)
      </ul>
      </ul>
    */  
    class Manager : public LOCA::Predictor::Generic {

    public:

      //! Constructor. 
      Manager(NOX::Parameter::List& params);

      //! Destructor
      virtual ~Manager();

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      //! Computes predictor based on predictor method parameter
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, double stepSize,
	      LOCA::Continuation::ExtendedGroup& prevGroup,
	      LOCA::Continuation::ExtendedGroup& curGroup,
	      LOCA::Continuation::ExtendedVector& result);

      //! Computes predictor based on predictor method parameter
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result);

      //! Returns current predictor method
      const string& getMethod() const;

    protected:
      //! Stores predictor method name
      string method;

      //! Stores pointer to current predictor object
      LOCA::Predictor::Generic *predictorPtr;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
