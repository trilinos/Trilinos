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

#ifndef LOCA_BLAS_GROUP_H
#define LOCA_BLAS_GROUP_H

#include "LOCA_Abstract_Group.H"	// base class
#include "NOX_BLAS_Group.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_BLAS_Interface.H"

namespace LOCA {
  namespace BLAS {

    class Group : public NOX::BLAS::Group, public LOCA::Abstract::Group {

    public:

      //! Constructor
      Group(Interface& i);

      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~Group();

      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      virtual Abstract::Group& operator=(const Abstract::Group& source);
      virtual Group& operator=(const Group& source);

      virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      virtual bool computeF();

      virtual bool computeJacobian();

      virtual bool computeTangent(NOX::Parameter::List& params, int paramId);

      //! Set the bifurcation/continuation parameter list in the group.
      virtual bool setParams(const ParameterVector& p);

      //! Compute a new parameter vector where this.pVector = grp.pVector + step*d.
      virtual bool computeParams(const ParameterVector& oldParams, 
				 const ParameterVector& direction, 
				 double step);

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const ParameterVector& getParams() const;

      //! Return Tangent vector (for predictor).
      virtual const NOX::Abstract::Vector& getTangent() const;

      //! Print out the group
      virtual bool print() const;

    protected:
      Interface& locaProblemInterface; 
      ParameterVector params;
      NOX::BLAS::Vector tangentVec;
      bool isValidTangent;
    };

  } // namespace BLAS
} // namespace LOCA


#endif
