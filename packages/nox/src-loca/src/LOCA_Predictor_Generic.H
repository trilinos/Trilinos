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

#ifndef LOCA_PREDICTOR_GENERIC_H
#define LOCA_PREDICTOR_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List

// Forward declarations
namespace LOCA {
  namespace Continuation {
    class ExtendedVector;
    class ExtendedGroup;
  }
  namespace MultiContinuation {
    class ExtendedGroup;
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  //! %Predictor direction strategies
  namespace Predictor {
    
    //! %Generic predictor direction interface.
    /*! 
     * %Generic interface for computing a predictor for continuation and 
     * serves as a base class for all predictor classes.
     *
     * Also provides a generic method for determining the orientation of
     * the predictor vector by comparing to the secant vector.
     */
    class Generic {

    public:

      //! Constructor. 
      /*! 
       * Constructors of derived objects should look like reset.
       */
      Generic();

      //! Destructor
      virtual ~Generic();

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      /*! 
       * \brief Compute the predictor direction given the current and
       * previous solution groups.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, double stepSize,
	      LOCA::Continuation::ExtendedGroup& prevGroup,
	      LOCA::Continuation::ExtendedGroup& curGroup,
	      LOCA::Continuation::ExtendedVector& result) = 0;

      /*! 
       * \brief Compute the predictor direction given the current and
       * previous solution groups.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result) = 0;

    protected:
      /*! 
       * \brief Sets orientation of predictor based on parameter change 
       * from previous steps
       */
      /*!
       * The implementation here looks at the sign of the scaled dot product
       * between the secant vector and the predictor and changes the sign
       * of the predictor if this scaled dot product is a different sign
       * than the step size.  If \em baseOnSecant is \em false, then the
       * sign is chosen so the parameter component of the predictor is 
       * positive for cases when a secant vector is not available (first
       * step in a continuation run) or may give incorrect information
       * (last step of a continuation run).
       */
      virtual void
      setPredictorOrientation(bool baseOnSecant, double stepSize,
			      LOCA::Continuation::ExtendedGroup& prevGroup,
			      LOCA::Continuation::ExtendedGroup& curGroup,
			      LOCA::Continuation::ExtendedVector& result);

      /*! 
       * \brief Sets orientation of predictor based on parameter change 
       * from previous steps
       */
      /*!
       * The implementation here looks at the sign of the scaled dot product
       * between the secant vector and the predictor and changes the sign
       * of the predictor if this scaled dot product is a different sign
       * than the step size.  If \em baseOnSecant is \em false, then the
       * sign is chosen so the parameter component of the predictor is 
       * positive for cases when a secant vector is not available (first
       * step in a continuation run) or may give incorrect information
       * (last step of a continuation run).
       */
      virtual void
      setPredictorOrientation(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result);

    protected:
      /*!
       * \brief Holds secant vector for setting predictor orientation based on
       * the secant vector orientation.
       */
      LOCA::Continuation::ExtendedVector *secantVecPtr;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
