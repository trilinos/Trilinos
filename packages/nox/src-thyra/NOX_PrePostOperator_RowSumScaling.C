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

#include "NOX_PrePostOperator_RowSumScaling.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Thyra_Group.H"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Teuchos_Assert.hpp"


NOX::RowSumScaling::
RowSumScaling(const Teuchos::RCP< ::Thyra::VectorBase<double> >& inv_row_sum_vec,
	      ENOX_WhenToUpdateScaling s) :
  inv_row_sum_vec_(inv_row_sum_vec),
  when_to_update(s)
{

}

void NOX::RowSumScaling::
runPreIterate(const NOX::Solver::Generic& solver)
{  
  if (when_to_update == UpdateInvRowSumVectorAtBeginningOfIteration)
    computeScaling(solver);
}

void NOX::RowSumScaling::
runPreSolve(const NOX::Solver::Generic& solver)
{  
  if (when_to_update == UpdateInvRowSumVectorAtBeginningOfSolve)
    computeScaling(solver);
}

Teuchos::RCP<const ::Thyra::VectorBase<double> > 
NOX::RowSumScaling::getInvRowSumScalingVector() const
{
  return inv_row_sum_vec_;
}


void NOX::RowSumScaling::
computeScaling(const NOX::Solver::Generic& solver)
{  
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<const NOX::Abstract::Group> group = solver.getSolutionGroupPtr();
  
  RCP<const NOX::Thyra::Group> thyra_group = 
    rcp_dynamic_cast<const NOX::Thyra::Group>(group);
  
  if (!thyra_group->isJacobian()) {
    RCP<NOX::Thyra::Group> tmp_nox_thyra_group = 
      Teuchos::rcp_const_cast<NOX::Thyra::Group>(thyra_group);
    TEUCHOS_ASSERT( !tmp_nox_thyra_group.is_null() );
    tmp_nox_thyra_group->computeJacobian();
  }
  
  RCP< const ::Thyra::LinearOpBase< double > > jac = 
    thyra_group->getJacobianOperator(); 	

  RCP< const ::Thyra::RowStatLinearOpBase< double > > row_stat_jac = 
    Teuchos::rcp_dynamic_cast< const ::Thyra::RowStatLinearOpBase< double > >(jac);

  TEUCHOS_ASSERT( !row_stat_jac.is_null() );

  if (inv_row_sum_vec_.is_null())
    inv_row_sum_vec_ = ::Thyra::createMember(jac->range());

  row_stat_jac->getRowStat( ::Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
			    inv_row_sum_vec_.ptr());

}
