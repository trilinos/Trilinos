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

#ifndef LOCA_MULTICONTINUATION_FINITEDIFFERENCEGROUP_H
#define LOCA_MULTICONTINUATION_FINITEDIFFERENCEGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class
#include "LOCA_DerivUtils.H"                    // class data element

namespace LOCA { 

  namespace MultiContinuation { 

    /*!  
     * \brief Concrete class that provides a concrete implementation of the
     * computeDfDp() method of the LOCA::Continuation::AbstractGroup using
     * first-order finite differencing.
    */
    /*! 
     * The finite-differencing calculations are actually implemented by
     * the LOCA::DerivUtils class, and a custom DerivUtils object can be
     * passed through the constructor of this class.  However, in the future
     * the calculations encapsulated in the DerivUtils class may be 
     * incorporated directly into this class and other finite-differencing
     * child classes.
    */

    class FiniteDifferenceGroup : 
      public virtual LOCA::MultiContinuation::AbstractGroup {

    public:
  
      //! Constructor.
      FiniteDifferenceGroup(const DerivUtils& deriv);
 
      //! Copy constructor
      FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~FiniteDifferenceGroup();

      //! Assignment operator
      virtual FiniteDifferenceGroup& 
      operator=(const FiniteDifferenceGroup& source);

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      /*! 
       * The calculation is implemented by calling the corresponding
       * LOCA::DerivUtils::computeDfDp() method of the passed 
       * LOCA::DerivUtils object.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const vector<int>& paramIDs, 
		       NOX::Abstract::MultiVector& dfdp, 
		       bool isValidF);

    protected:

      //! Pointer to current DerivUtils derivative computation object
      DerivUtils* derivPtr;

    };
  } // namespace MultiContinuation
} // namespace LOCA

#endif
