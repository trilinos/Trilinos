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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
