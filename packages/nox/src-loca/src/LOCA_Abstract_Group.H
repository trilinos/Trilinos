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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "NOX_Abstract_Group.H"    // base class

#include "NOX_Common.H"            // for string 
#include "LOCA_DerivUtils.H"

namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace LOCA {
  class ParameterVector;
}

//! Library of continuation algorithms namespace
namespace LOCA { 

  //! %LOCA abstract interface namespace
  namespace Abstract { 

    /*!  \brief %LOCA abstract interface to a group, derived from the 
      NOX::Abstract::Group.  %LOCA %abstract groups provide the interface 
      necessary to perform continuation, i.e., compute families of solutions to
      \f$ F(x,p) = 0 \f$.
    */
    /*! 
      Concrete implemenations of this interface must provide implementations of
      all of the methods in the NOX::Abstract::Group interface as well as the
      parameter and scale vector setting/retrieval methods.
    */

    class Group : public virtual NOX::Abstract::Group {

    public:
  
      //! Default constructor.
      Group(const DerivUtils& deriv = DerivUtils());
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();
  
      //! Assignment operator
      virtual Group& operator=(const Group& source);

      //! Generic routine for applying Nic method to solve near-singular system
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseNic(NOX::Parameter::List& params,
			      const NOX::Abstract::Vector& input,
			      const NOX::Abstract::Vector& approxNullVec,
			      const NOX::Abstract::Vector& JacApproxNullVec,
			      NOX::Abstract::Vector& result) const;

      /*! Generic routine for applying Nic method to solve near-singular system
	with multiple right hand sides */
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseNicMulti(NOX::Parameter::List& params,
			      const NOX::Abstract::Vector*const* inputs,
			      const NOX::Abstract::Vector& approxNullVec,
			      const NOX::Abstract::Vector& JacApproxNullVec,
			      NOX::Abstract::Vector** results, 
			      int nVecs) const;

      //! Generic routine for applying Nic method to solve near-singular system
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseNicDay(NOX::Parameter::List& params,
				 const NOX::Abstract::Vector& input,
				 const NOX::Abstract::Vector& approxNullVec,
				 const NOX::Abstract::Vector& JacApproxNullVec,
				 NOX::Abstract::Vector& result) const;

      /*! Generic routine for applying Nic method to solve near-singular system
	with multiple right hand sides */
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseNicDayMulti(NOX::Parameter::List& params,
				 const NOX::Abstract::Vector*const* inputs,
				 const NOX::Abstract::Vector& approxNullVec,
				 const NOX::Abstract::Vector& JacApproxNullVec,
				 NOX::Abstract::Vector** results,
				 int nVecs) const;

      //! Generic routine for applying Iterative Refinement for solving near-singular system
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseItRef(NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input,
				NOX::Abstract::Vector& result) const;


      //! Generic routine for applying Iterative Refinement for solving near-singular system with multiple right hand sides
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseItRefMulti(NOX::Parameter::List& params,
				const NOX::Abstract::Vector*const* inputs,
				NOX::Abstract::Vector** results, 
				int nVecs) const;

      //! applyJacobianInverse for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseMulti(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const;

      //! Set the parameter vector in the group to p (pVector = p).  
      virtual void setParams(const ParameterVector& p) = 0;

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const ParameterVector& getParams() const = 0;

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val) = 0;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const = 0;

      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val) = 0;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const = 0;

      //! Compute \f$dF/dp\f$ using DerivUtils derivative object
      virtual NOX::Abstract::Group::ReturnType
      computeDfDp(int paramID, NOX::Abstract::Vector& result);

      //! For debugging purposes only
      virtual void print() const = 0;

      //! Set scale vector
      virtual void setScaleVec(const NOX::Abstract::Vector& s) = 0;

      //! Get scale vector
      virtual const NOX::Abstract::Vector& getScaleVec() const = 0;
  
    protected:
      //! Pointer to current DerivUtils derivative computation object
      DerivUtils* derivPtr;
    };
  } // namespace Abstract
} // namespace LOCA

#endif
