/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "EpetraExt_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"


namespace EpetraExt {


DiagonalQuadraticResponseOnlyModelEvaluator
::DiagonalQuadraticResponseOnlyModelEvaluator(
  const Teuchos::RCP<Epetra_Comm> &comm,
  const int localDim, const double &pt, const double &p0, const double &scale
  )
	:epetra_comm_(comm), scale_(scale)
{

  using Teuchos::rcp;

  const int ng = 1;

  map_p_ = rcp(new Epetra_Map(-1, localDim, 0, *epetra_comm_));
  map_g_ = rcp(new Epetra_Map(ng, ng, 0, *epetra_comm_));

  pt_ = rcp(new Epetra_Vector(*map_p_));
  pt_->PutScalar(pt);

  p0_ = rcp(new Epetra_Vector(*map_p_));
  p0_->PutScalar(p0);

}


// Overridden from EpetraExt::ModelEvaluator


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalQuadraticResponseOnlyModelEvaluator::get_x_map() const
{
  return Teuchos::null;
}


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalQuadraticResponseOnlyModelEvaluator::get_f_map() const
{
  return Teuchos::null;
}


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalQuadraticResponseOnlyModelEvaluator::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return map_p_;
}


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalQuadraticResponseOnlyModelEvaluator::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=0);
  return map_g_;
}


Teuchos::RefCountPtr<const Epetra_Vector>
DiagonalQuadraticResponseOnlyModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return p0_;
}


EpetraExt::ModelEvaluator::InArgs
DiagonalQuadraticResponseOnlyModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  return inArgs;
}


EpetraExt::ModelEvaluator::OutArgs
DiagonalQuadraticResponseOnlyModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    0, 0, DerivativeProperties(
      DERIV_LINEARITY_NONCONST,
      DERIV_RANK_DEFICIENT,
      true // supportsAdjoint
      )
    );
  return outArgs;
}


void DiagonalQuadraticResponseOnlyModelEvaluator::evalModel(
  const InArgs& inArgs, const OutArgs& outArgs
  ) const
{

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

  //
  // Get the input arguments
  //

  const Epetra_Vector &p = *inArgs.get_p(0);

  //
  // Get the output arguments
  //

  const RCP<Epetra_Vector> g_out = outArgs.get_g(0);

  const RCP<Epetra_MultiVector> DgDp_trans_out =
    get_DgDp_mv(0, 0,outArgs,DERIV_TRANS_MV_BY_ROW);

  //
  // Compute the functions
  //

  if (nonnull(g_out) || nonnull(DgDp_trans_out)) {

    Epetra_Vector p_minus_pt(*map_p_);

    p_minus_pt = p;
    p_minus_pt.Update(-1.0, *pt_, 1.0);

    if (nonnull(g_out)) {
      double dot[1];
      p_minus_pt.Dot(p_minus_pt, dot);
      (*g_out)[0] = scale_ * 0.5 * dot[0];
    }
    
    if (nonnull(DgDp_trans_out)) {
      (*DgDp_trans_out) = p_minus_pt;
      DgDp_trans_out->Scale(scale_);
    }

  }

}


} // namespace EpetraExt
