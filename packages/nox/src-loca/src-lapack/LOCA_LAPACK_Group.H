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

#ifndef LOCA_LAPACK_GROUP_H
#define LOCA_LAPACK_GROUP_H

#include "LOCA_Abstract_Group.H"	// base class
#include "NOX_LAPACK_Group.H"           // base class
#include "LOCA_Parameter_Vector.H"      // class data element
#include "LOCA_LAPACK_Interface.H"      // class data element

namespace LOCA {

  //! %LOCA BLAS/%LAPACK support
  namespace LAPACK {

    //! Extension of the NOX::LAPACK::Group to LOCA.  
    /*!
      This class is derived both from the NOX::LAPACK::Group and
      LOCA::Abstract::Group classes and therefore inherits the implementation
      of the NOX::Abstract::Group interface provided by NOX::LAPACK::Group.

      This class additionally stores a ParameterVector params and a scaling
      vector (NOX::LAPACK::Vector) scaleVec and implements the pure virtual
      methods of the LOCA::Abstract::Group class.  

      The computeF() and computeJacobian() methods of the NOX::LAPACK::Group
      parent class are overloaded here.  They both set the entire contents
      of the parameter vector in the problem interface before calling the
      NOX::LAPACK::Group computeF() and computeJacobian().
    */
    
    class Group : public NOX::LAPACK::Group, public LOCA::Abstract::Group {

    public:

      //! Constructor with default scaling vector set to all one's
      Group(LOCA::LAPACK::Interface& i);

      //! Constructor with specified scaling vector
      Group(LOCA::LAPACK::Interface& i, const NOX::LAPACK::Vector& s);

      //! Copy constructor
      Group(const LOCA::LAPACK::Group& source, 
	    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~Group();

      //! Assignment operator
      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      //! Assignment operator
      LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

      //! Assignment operator
      NOX::LAPACK::Group& operator=(const NOX::LAPACK::Group& source);

      //! Assignment operator
      LOCA::LAPACK::Group& operator=(const LOCA::LAPACK::Group& source);

      //! Cloning function
      NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      /*!
	Calls LOCA::LAPACK::Interface::setParams before evalulating F.
      */
      NOX::Abstract::Group::ReturnType computeF();

      //! Overloaded computeJacobian()
      /*!
	Calls LOCA::LAPACK::Interface::setParams before evalulating J.
      */
      NOX::Abstract::Group::ReturnType computeJacobian();

      //! applyJacobianInverse for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseMulti(NOX::Parameter::List& params,
			   const NOX::Abstract::Vector* const* inputs,
			   NOX::Abstract::Vector** outputs, int nVecs) const;

      //! Set the parameter vector
      void setParams(const LOCA::ParameterVector& p);

      //! Return a const reference to the parameter vector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //! Print out the solution vector and continuation parameter
      void printSolution(const double conParam) const;

      //! Print out a vector and a parameter
      void printSolution(const NOX::LAPACK::Vector& x_, const double conParam) const;

      //! Print out a vector and a parameter
      void printSolution(const NOX::Abstract::Vector& x_, const double conParam) const;

      //! Set scale vector
      virtual void setScaleVec(const NOX::Abstract::Vector& s);

      //! Set scale vector
      void setScaleVec(const NOX::LAPACK::Vector& s);

      //! Get scale vector
      virtual const NOX::Abstract::Vector& getScaleVec() const;

      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

    protected:

      //! Compute a default scale vector (vector of all one's)
      void computeScaleVec();

    protected:
      //! Referece to current problem
      LOCA::LAPACK::Interface& locaProblemInterface;

      //! vector of parameters
      ParameterVector params;

      //! Scaling vector
      NOX::LAPACK::Vector scaleVec;
    };

  } // namespace LAPACK
} // namespace LOCA


#endif
