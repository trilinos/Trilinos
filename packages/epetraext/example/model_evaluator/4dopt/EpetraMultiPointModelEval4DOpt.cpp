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

#include "EpetraMultiPointModelEval4DOpt.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

#include "Epetra_MpiComm.h"

namespace {

inline double sqr( const double& s ) { return s*s; }

} // namespace

EpetraMultiPointModelEval4DOpt::EpetraMultiPointModelEval4DOpt(
  Teuchos::RefCountPtr<Epetra_Comm> epetra_comm
  ,const double         xt0
  ,const double        xt1
  ,const double        pt0
  ,const double        pt1
  ,const double        d
  ,const double        x00
  ,const double        x01
  ,const double        p00
  ,const double        p01
  ,const double        q0
  )
  :isInitialized_(false), epetra_comm_(epetra_comm),
   xt0_(xt0),xt1_(xt1),pt0_(pt0),pt1_(pt1),d_(d)
{
  using Teuchos::rcp;

  const int nx = 2, np = 2, ng = 1, nq = 1;

  map_x_ = rcp(new Epetra_Map(nx,0,*epetra_comm_));
  map_p_ = rcp(new Epetra_Map(np,0,*epetra_comm_));
  map_q_ = rcp(new Epetra_Map(nq,0,*epetra_comm_));
  map_g_ = rcp(new Epetra_Map(ng,0,*epetra_comm_));

  //const double inf = infiniteBound();
  const double inf = 1e+50;

  xL_ = rcp(new Epetra_Vector(*map_x_));  xL_->PutScalar(-inf);
  xU_ = rcp(new Epetra_Vector(*map_x_));  xU_->PutScalar(+inf);
  x0_ = rcp(new Epetra_Vector(*map_x_));  (*x0_)[0] = x00; (*x0_)[1] = x01;
  pL_ = rcp(new Epetra_Vector(*map_p_));  pL_->PutScalar(-inf);
  pU_ = rcp(new Epetra_Vector(*map_p_));  pU_->PutScalar(+inf);
  p0_ = rcp(new Epetra_Vector(*map_p_));  (*p0_)[0] = p00; (*p0_)[1] = p01;
  gL_ = rcp(new Epetra_Vector(*map_g_));  gL_->PutScalar(-inf);
  gU_ = rcp(new Epetra_Vector(*map_g_));  gU_->PutScalar(+inf);
  q_  = rcp(new Epetra_Vector(*map_q_));  (*q_)[0] = q0;
  qL_ = rcp(new Epetra_Vector(*map_q_));  (*qL_)[0] = -inf;
  qU_ = rcp(new Epetra_Vector(*map_q_));  (*qU_)[0] =  inf;

  // Initialize the graph for W CrsMatrix object
  W_graph_ = rcp(new Epetra_CrsGraph(::Copy,*map_x_,nx));
  {
    int indices[nx] = { 0, 1 };
    for( int i = 0; i < nx; ++i )
      W_graph_->InsertGlobalIndices(i,nx,indices);
  }
  W_graph_->FillComplete();

  isInitialized_ = true;

}

void EpetraMultiPointModelEval4DOpt::set_p_bounds(
  double pL0, double pL1, double pU0, double pU1
  )
{
  (*pL_)[0] = pL0;
  (*pL_)[1] = pL1;
  (*pU_)[0] = pU0;
  (*pU_)[1] = pU1;
}

void EpetraMultiPointModelEval4DOpt::set_x_bounds(
  double xL0, double xL1, double xU0, double xU1
  )
{
  (*xL_)[0] = xL0;
  (*xL_)[1] = xL1;
  (*xU_)[0] = xU0;
  (*xU_)[1] = xU1;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
EpetraMultiPointModelEval4DOpt::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraMultiPointModelEval4DOpt::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraMultiPointModelEval4DOpt::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(l>1);
  if (l==0) return map_p_;
  else return map_q_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraMultiPointModelEval4DOpt::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=0);
  return map_g_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(l>1);
  if (l==0) return p0_;
  else return q_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_x_lower_bounds() const
{
  return xL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_x_upper_bounds() const
{
  return xU_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_p_lower_bounds(int l) const
{
  TEST_FOR_EXCEPT(l>1);
  if (l==0) return pL_;
  else      return qL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraMultiPointModelEval4DOpt::get_p_upper_bounds(int l) const
{
  TEST_FOR_EXCEPT(l>1);
  if (l==0) return pU_;
  else      return qU_;
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraMultiPointModelEval4DOpt::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

EpetraExt::ModelEvaluator::InArgs
EpetraMultiPointModelEval4DOpt::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(2);
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
EpetraMultiPointModelEval4DOpt::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(2,1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DfDp,0,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_CONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDx,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,0,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    0,0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

void EpetraMultiPointModelEval4DOpt::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //
  Teuchos::RefCountPtr<const Epetra_Vector> p_in = inArgs.get_p(0);
  Teuchos::RefCountPtr<const Epetra_Vector> q_in = inArgs.get_p(1);
  const Epetra_Vector &p = (p_in.get() ? *p_in : *p0_);
  const Epetra_Vector &q = (q_in.get() ? *q_in : *q_);
  const Epetra_Vector &x = *inArgs.get_x();
  //
  // Get the output arguments
  //
  Epetra_Vector       *f_out = outArgs.get_f().get();
  Epetra_Vector       *g_out = outArgs.get_g(0).get();
  Epetra_Operator     *W_out = outArgs.get_W().get();
  Epetra_MultiVector  *DfDp_out = get_DfDp_mv(0,outArgs).get();
  Epetra_MultiVector  *DgDx_trans_out = get_DgDx_mv(0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
  Epetra_MultiVector  *DgDp_trans_out = get_DgDp_mv(0,0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
  //
  // Compute the functions
  //
  if(f_out) {
    Epetra_Vector &f = *f_out;
    // zero-based indexing for Epetra!
    //q[0] added for multipoint problem!
    f[0] =        x[0]      + x[1]*x[1] - p[0] + q[0];
    f[1] = d_ * ( x[0]*x[0] - x[1]      - p[1] );
  }
  if(g_out) {
    Epetra_Vector &g = *g_out;
    g[0] = 0.5 * ( sqr(x[0]-xt0_) + sqr(x[1]-xt1_) + sqr(p[0]-pt0_) + sqr(p[1]-pt1_) );
  }
  if(W_out) {
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_out);
    DfDx.PutScalar(0.0);
    //
    // Fill W = DfDx
    //
    // W = DfDx = [      1.0,  2*x[1] ]
    //            [ 2*d*x[0],     -d  ]
    //
    double values[2];
    int indexes[2];
    // Row [0]
    values[0] = 1.0;           indexes[0] = 0;
    values[1] = 2.0*x[1];      indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
    // Row [1]
    values[0] = 2.0*d_*x[0];   indexes[0] = 0;
    values[1] = -d_;           indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
  }
  if(DfDp_out) {
    Epetra_MultiVector &DfDp = *DfDp_out;
    DfDp.PutScalar(0.0);
    DfDp[0][0] = -1.0;
    DfDp[1][1] = -d_;
  }
  if(DgDx_trans_out) {
    Epetra_Vector &DgDx_trans = *(*DgDx_trans_out)(0);
    DgDx_trans[0] = x[0]-xt0_;
    DgDx_trans[1] = x[1]-xt1_;
  }
  if(DgDp_trans_out) {
    Epetra_Vector &DgDp_trans = *(*DgDp_trans_out)(0);
    DgDp_trans[0] = p[0]-pt0_;
    DgDp_trans[1] = p[1]-pt1_;
  }
}
