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
  namespace LAPACK {

    class Group : public NOX::LAPACK::Group, public LOCA::Abstract::Group {

    public:

      //! Constructor
      Group(Interface& i);

      //! Copy constructor
      Group(const LOCA::LAPACK::Group& source, 
	    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~Group();

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);
      NOX::LAPACK::Group& operator=(const NOX::LAPACK::Group& source);
      LOCA::LAPACK::Group& operator=(const LOCA::LAPACK::Group& source);

      NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeTangent(NOX::Parameter::List& params, int paramId);

      //! Set the bifurcation/continuation parameter list in the group.
      void setParams(const LOCA::ParameterVector& p);

      //! Compute a new parameter vector where this.pVector = grp.pVector + step*d.
      void computeParams(const LOCA::ParameterVector& oldParams, 
			 const LOCA::ParameterVector& direction, double step);

      //! Return a const reference to the ParameterVector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return Tangent vector (for predictor).
      const NOX::Abstract::Vector& getTangent() const;

      //! Print out the group
      void print() const;

    protected:
      //! resets the isValid flags to false
      void resetIsValid();

    protected:
      Interface& locaProblemInterface; 
      ParameterVector params;
      NOX::LAPACK::Vector tangentVec;
      bool isValidTangent;
    };

  } // namespace LAPACK
} // namespace LOCA


#endif
