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

#ifndef LOCA_CONTINUATION_FINITEDIFFERENCEGROUP_H
#define LOCA_CONTINUATION_FINITEDIFFERENCEGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class
#include "LOCA_DerivUtils.H"                    // class data element

namespace LOCA { 

  namespace Continuation { 

    /*!  \brief %LOCA abstract interface to a group, derived from the 
      NOX::Abstract::Group.  %LOCA %continuation underlying groups provide 
      the interface necessary to perform continuation, i.e., compute families 
      of solutions to \f$ F(x,p) = 0 \f$.
    */
    /*! 
      Concrete implemenations of this interface must provide implementations of
      all of the methods in the NOX::Abstract::Group interface as well as the
      parameter and scale vector setting/retrieval methods.
    */

    class FiniteDifferenceGroup : 
      public virtual LOCA::Continuation::AbstractGroup {

    public:
  
      //! Default constructor.
      FiniteDifferenceGroup(const DerivUtils& deriv);
 
      //! Copy constructor
      FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~FiniteDifferenceGroup();

      //! Assignment operator
      virtual FiniteDifferenceGroup& 
      operator=(const FiniteDifferenceGroup& source);

      //! Compute \f$dF/dp\f$ using DerivUtils derivative object
      virtual NOX::Abstract::Group::ReturnType
      computeDfDp(int paramID, NOX::Abstract::Vector& result);

    protected:

      //! Pointer to current DerivUtils derivative computation object
      DerivUtils* derivPtr;

    };
  } // namespace Continuation
} // namespace LOCA

#endif
