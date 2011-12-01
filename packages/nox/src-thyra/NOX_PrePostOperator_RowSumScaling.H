// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H
#define NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_PrePostOperator.H"
#include "Thyra_VectorBase.hpp"

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

  /** \brief Updates an inverse row sum scaling vector at the beginning of a solve.

  */
  class RowSumScaling : public NOX::Abstract::PrePostOperator {
    
  public:
    
    enum ENOX_WhenToUpdateScaling {
      UpdateInvRowSumVectorAtBeginningOfSolve,
      UpdateInvRowSumVectorAtBeginningOfIteration
    };

    /** <tt>sv</tt> is the scaling vector */
    RowSumScaling(const Teuchos::RCP< ::Thyra::VectorBase<double> >& inv_row_sum_vec,
		  ENOX_WhenToUpdateScaling s);

    void runPreIterate(const NOX::Solver::Generic& solver);
    
    void runPreSolve(const NOX::Solver::Generic& solver);
    
    Teuchos::RCP<const ::Thyra::VectorBase<double> > 
    getInvRowSumScalingVector() const;

  private:

    void computeScaling(const NOX::Solver::Generic& solver);

    //! Inverse row sum vector
    Teuchos::RCP< ::Thyra::VectorBase<double> > inv_row_sum_vec_;
    
    ENOX_WhenToUpdateScaling when_to_update;
  };
}

#endif
