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

#ifndef LOCA_MULTICONTINUATION_ABSTRACTGROUP_H
#define LOCA_MULTICONTINUATION_ABSTRACTGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class
#include "NOX_Abstract_MultiVector.H"

namespace LOCA { 

  //! Groups and vectors for multi-parameter continuation
  namespace MultiContinuation { 

    /*! 
     * \brief %LOCA abstract interface for continuation, derived from the 
     * NOX::Abstract::Group.  This abstract class provides 
     * the interface necessary to perform continuation, i.e., compute 
     * families of solutions to \f$ F(x,p) = 0 \f$.
    */
    /*! 
     * Concrete implemenations of this interface must provide implementations 
     * of all of the methods in the NOX::Abstract::Group interface as well as 
     * the additional interface defined here.
    */

    class AbstractGroup : public virtual LOCA::Continuation::AbstractGroup {

    public:
  
      //! Default constructor.
      AbstractGroup() {}
 
      //! Destructor.
      virtual ~AbstractGroup() {}

      /*! 
       * @name Pure virtual methods 
       * These methods must be defined by any concrete implementation
       */
      //@{

      //! Assignment operator
      virtual AbstractGroup& operator=(const AbstractGroup& source) = 0;

       //! Set parameters indexed by (integer) paramIDs
      virtual void setParams(
		     const vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals) = 0;

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDp(const vector<int>& paramIDs, 
		  NOX::Abstract::MultiVector& dfdp, 
		  bool isValidF) = 0;

      //@}

    };
  } // namespace MultiContinuation
} // namespace LOCA

#endif
