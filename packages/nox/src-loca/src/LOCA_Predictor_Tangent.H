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
  
      //! Computes the tangent predictor.
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, double stepSize,
	      LOCA::Continuation::ExtendedGroup& prevGroup,
	      LOCA::Continuation::ExtendedGroup& curGroup,
	      LOCA::Continuation::ExtendedVector& result);

      //! Computes the tangent predictor.
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result);

    protected:

      //! Stores derivative of residual w.r.t. parameter
      NOX::Abstract::Vector* dfdpVecPtr;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
