#include "EpetraModelEval4DOpt.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

namespace {

inline double sqr( const double& s ) { return s*s; }

} // namespace

EpetraModelEval4DOpt::EpetraModelEval4DOpt(
  const double         xt0
  ,const double        xt1
  ,const double        pt0
  ,const double        pt1
  ,const double        d
  ,const double        x00
  ,const double        x01
  ,const double        p00
  ,const double        p01
  )
  :xt0_(xt0),xt1_(xt1),pt0_(pt0),pt1_(pt1),d_(d)
{
  using Teuchos::rcp;

	epetra_comm_ = rcp(new Epetra_SerialComm());

  const int nx = 2, np = 2, ng = 1;

  map_x_ = rcp(new Epetra_Map(nx,0,*epetra_comm_));
  map_p_ = rcp(new Epetra_Map(np,0,*epetra_comm_));
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

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval4DOpt::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval4DOpt::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval4DOpt::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return map_p_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval4DOpt::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=1);
  return map_g_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return p0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_x_lower_bounds() const
{
  return xL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_x_upper_bounds() const
{
  return xU_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_p_lower_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return pL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval4DOpt::get_p_upper_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return pU_;
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraModelEval4DOpt::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

EpetraExt::ModelEvaluator::InArgs
EpetraModelEval4DOpt::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_beta,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
EpetraModelEval4DOpt::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1,1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DfDp,1,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    1,DerivativeProperties(
      DERIV_LINEARITY_CONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDx,1,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    1,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,1,1,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    1,1,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

void EpetraModelEval4DOpt::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //
  Teuchos::RefCountPtr<const Epetra_Vector> p_in = inArgs.get_p(1);
  const Epetra_Vector &p = (p_in.get() ? *p_in : *p0_);
  const Epetra_Vector &x = *inArgs.get_x();
  //
  // Get the output arguments
  //
  Teuchos::RefCountPtr<Epetra_Vector>       f_inout = outArgs.get_f();
  Teuchos::RefCountPtr<Epetra_Vector>       g_inout = outArgs.get_g(1);
  Teuchos::RefCountPtr<Epetra_Operator>     W_inout = outArgs.get_W();
  Teuchos::RefCountPtr<Epetra_MultiVector>  DfDp_inout = get_DfDp_mv(1,outArgs);
  Teuchos::RefCountPtr<Epetra_MultiVector>  DgDx_trans_inout = get_DgDx_mv(1,outArgs,DERIV_TRANS_MV_BY_ROW);
  Teuchos::RefCountPtr<Epetra_MultiVector>  DgDp_trans_inout = get_DgDp_mv(1,1,outArgs,DERIV_TRANS_MV_BY_ROW);
  //
  // Compute the functions
  //
  if(f_inout.get()) {
    Epetra_Vector &f = *f_inout;
    // zero-based indexing for Epetra!
    f[0] =        x[0]      + x[1]*x[1] - p[0];
    f[1] = d_ * ( x[0]*x[0] - x[1]      - p[1] );
  }
  if(g_inout.get()) {
    Epetra_Vector &g = *g_inout;
    // zero-based indexing for Epetra, 1-based indexing for Thyra!
    g[0] = 0.5 * ( sqr(x[0]-xt0_) + sqr(x[1]-xt1_) + sqr(p[0]-pt0_) + sqr(p[1]-pt1_) );
  }
  if(W_inout.get()) {
    const double beta = inArgs.get_beta();
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_inout);
    DfDx.PutScalar(0.0);
    //
    // Fill W = beta*DfDx
    //
    // W = beta*DfDx = beta * [      1.0,  2*x[1] ]
    //                        [ 2*d*x[0],     -d  ]
    //
    double values[2];
    int indexes[2];
    // Row [0]
    values[0] = beta;               indexes[0] = 0;
    values[1] = 2.0*beta*x[1];      indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
    // Row [1]
    values[0] = 2.0*beta*d_*x[0];   indexes[0] = 0;
    values[1] = -beta*d_;           indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
  }
  if(DfDp_inout.get()) {
    Epetra_MultiVector &DfDp = *DfDp_inout;
    DfDp.PutScalar(0.0);
    DfDp[0][0] = -1.0;
    DfDp[1][1] = -d_;
  }
  if(DgDx_trans_inout.get()) {
    Epetra_Vector &DgDx_trans = *(*DgDx_trans_inout)(0);
    DgDx_trans[0] = x[0]-xt0_;
    DgDx_trans[1] = x[1]-xt1_;
  }
  if(DgDp_trans_inout.get()) {
    Epetra_Vector &DgDp_trans = *(*DgDp_trans_inout)(0);
    DgDp_trans[0] = p[0]-pt0_;
    DgDp_trans[1] = p[1]-pt1_;
  }
}
