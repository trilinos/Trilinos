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

#ifndef LOCA_CONTINUATION_GROUP_H
#define LOCA_CONTINUATION_GROUP_H

#include "NOX_Parameter_List.H"
#include "NOX_Abstract_Group.H"        // base class
#include "LOCA_Abstract_Group.H"       // class data element
#include "LOCA_Continuation_Vector.H"  // class data element
#include "NOX_Common.H"                // for string 

namespace LOCA { 
  namespace Continuation { 

    class Group : public virtual NOX::Abstract::Group {

    public:
  
      //! Constructor.
      Group(const LOCA::Abstract::Group& grp, int paramID,
	    const NOX::Parameter::List& linSolverParams);
      Group(const LOCA::Abstract::Group& grp, string paramID,
	    const NOX::Parameter::List& linSolverParams);

      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      virtual Group& operator=(const Group& source);

      //! Set the previous solution vector, x_old, to y.
      virtual void setPrevX(const NOX::Abstract::Vector& y) = 0;
      virtual void setPrevX(const LOCA::Continuation::Vector& y) = 0;

      //! Gets the previous solution vector
      virtual const LOCA::Continuation::Vector& getPrevX() const = 0;

      virtual bool isPrevXVec() const = 0;

      //! Set step size for continuation constraint equation
      virtual void setStepSize(double deltaS) = 0;

      //! Compute tangent direction
      virtual NOX::Abstract::Group::ReturnType
      computeTangent(NOX::Parameter::List& params);

      //! Get tangent direction vector, needed for predictor
      virtual const LOCA::Continuation::Vector& getTangent() const;

      virtual bool isTangent() const;

      virtual void setContinuationParameter(double val);

      virtual double getContinuationParameter() const;

      virtual const LOCA::Abstract::Group& getGroup() const;

    protected:

      void resetIsValid();

    protected:

      LOCA::Abstract::Group* grpPtr;
      int conParamID;
      LOCA::Continuation::Vector tangentVec;
      bool isValidTangent;
      NOX::Parameter::List linearSolverParams;
  
    };
  } // namespace Abstract
} // namespace LOCA

#endif
