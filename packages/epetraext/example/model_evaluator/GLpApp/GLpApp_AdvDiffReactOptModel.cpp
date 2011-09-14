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

#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_TimeMonitor.hpp"

// 2007/06/14: rabartl: The dependence on Thyra is non-optional right now
// since I use it to perform a mat-vec that I could not figure out how to do
// with just epetra.  If this becomes a problem then we have to change this
// later!  Just grep for 'Thyra' outside of the ifdef for
// HAVE_THYRA_EPETRAEXT.
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"

#ifdef HAVE_THYRA_EPETRAEXT
// For orthogonalization of the basis B_bar
#include "sillyModifiedGramSchmidt.hpp" // This is just an example!
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#endif // HAVE_THYRA_EPETRAEXT

//#define GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF

namespace {

Teuchos::RefCountPtr<Teuchos::Time>
  initalizationTimer    = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Init"),
  evalTimer             = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval"),
  p_bar_Timer           = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:p_bar"),
  R_p_bar_Timer         = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:R_p_bar"),
  f_Timer               = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:f"),
  g_Timer               = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:g"),
  W_Timer               = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:W"),
  DfDp_Timer            = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:DfDp"),
  DgDx_Timer            = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:DgDx"),
  DgDp_Timer            = Teuchos::TimeMonitor::getNewTimer("AdvDiffReact:Eval:DgDp");

inline double sqr( const double& s ) { return s*s; }

inline double dot( const Epetra_Vector &x, const Epetra_Vector &y )
{
  double dot[1];
  x.Dot(y,dot);
  return dot[0];
}

} // namespace

namespace GLpApp {

AdvDiffReactOptModel::AdvDiffReactOptModel(
  const Teuchos::RefCountPtr<const Epetra_Comm>  &comm
  ,const double                                  beta
  ,const double                                  len_x     // Ignored if meshFile is *not* empty
  ,const double                                  len_y     // Ignored if meshFile is *not* empty
  ,const int                                     local_nx  // Ignored if meshFile is *not* empty
  ,const int                                     local_ny  // Ignored if meshFile is *not* empty
  ,const char                                    meshFile[]
  ,const int                                     np
  ,const double                                  x0
  ,const double                                  p0
  ,const double                                  reactionRate
  ,const bool                                    normalizeBasis
  ,const bool                                    supportDerivatives
  )
  : supportDerivatives_(supportDerivatives), np_(np)
{
  Teuchos::TimeMonitor initalizationTimerMonitor(*initalizationTimer);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out << "\nEntering AdvDiffReactOptModel::AdvDiffReactOptModel(...) ...\n";
#endif
  //
  // Initalize the data pool object
  //
  dat_ = Teuchos::rcp(new GLpApp::GLpYUEpetraDataPool(comm,beta,len_x,len_y,local_nx,local_ny,meshFile,false));
  //
  // Get the maps
  //
  map_x_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorDomainMap()));
  map_p_.resize(Np_);
  map_p_[p_rx_idx] = Teuchos::rcp(new Epetra_Map(1,1,0,*comm));
  map_p_bar_ = Teuchos::rcp(new Epetra_Map(dat_->getB()->OperatorDomainMap()));
  map_f_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorRangeMap()));
  map_g_ = Teuchos::rcp(new Epetra_Map(1,1,0,*comm));
  //
  // Initialize the basis matrix for p such that p_bar = B_bar * p
  //
  if( np_ > 0 ) {
    //
    // Create a sine series basis for y (the odd part of a Fourier series basis)
    //
    const Epetra_SerialDenseMatrix &ipcoords = *dat_->getipcoords();
    const Epetra_IntSerialDenseVector &pindx = *dat_->getpindx();
    Epetra_Map overlapmap(-1,pindx.M(),const_cast<int*>(pindx.A()),1,*comm);
    const double pi = 2.0 * std::asin(1.0);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *out << "\npi="<<pi<<"\n";
#endif
    map_p_[p_bndy_idx] = Teuchos::rcp(new Epetra_Map(np_,np_,0,*comm));
    B_bar_ = Teuchos::rcp(new Epetra_MultiVector(*map_p_bar_,np_));
    (*B_bar_)(0)->PutScalar(1.0); // First column is all ones!
    if( np_ > 1 ) {
      const int numBndyNodes        = map_p_bar_->NumMyElements();
      const int *bndyNodeGlobalIDs  = map_p_bar_->MyGlobalElements();
      for( int i = 0; i < numBndyNodes; ++i ) {
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
        const int global_id = bndyNodeGlobalIDs[i];
#endif
        const int local_id = overlapmap.LID(bndyNodeGlobalIDs[i]);
        const double x = ipcoords(local_id,0), y = ipcoords(local_id,1);
        double z = -1.0, L = -1.0;
        if( x < 1e-10 || len_x - 1e-10 < x ) {
          z = y;
          L = len_y;
        }
        else {
          z = x;
          L = len_x;
        }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
        *out << "\ni="<<i<<",global_id="<<global_id<<",local_id="<<local_id<<",x="<<x<<",y="<<y<<",z="<<z<<"\n";
#endif
        for( int j = 1; j < np_; ++j ) {
          (*B_bar_)[j][i] = std::sin(j*pi*z/L);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
          *out << "\nB("<<i<<","<<j<<")="<<(*B_bar_)[j][i]<<"\n";
#endif
        }
      }
      if(normalizeBasis) {
#ifdef HAVE_THYRA_EPETRAEXT
        //
        // Use modified Gram-Schmidt to create an orthonormal version of B_bar!
        //
        Teuchos::RefCountPtr<Thyra::MultiVectorBase<double> >
          thyra_B_bar = Thyra::create_MultiVector(
            B_bar_
            ,Thyra::create_VectorSpace(Teuchos::rcp(new Epetra_Map(*map_p_bar_)))
            ,Thyra::create_VectorSpace(Teuchos::rcp(new Epetra_Map(*map_p_[p_bndy_idx])))
            ),
          thyra_fact_R;
        Thyra::sillyModifiedGramSchmidt(thyra_B_bar.ptr(), Teuchos::outArg(thyra_fact_R));
        Thyra::scale(double(numBndyNodes)/double(np_),  thyra_B_bar.ptr()); // Each row should sum to around one!
        // We just discard the "R" factory thyra_fact_R ...
#else // HAVE_THYRA_EPETRAEXT
        TEST_FOR_EXCEPTION(
          true,std::logic_error
          ,"Error, can not normalize basis since we do not have Thyra support enabled!"
          );
#endif // HAVE_EPETRAEXT_THYRA
      }
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *out << "\nB_bar =\n\n";
    B_bar_->Print(Teuchos::OSTab(out).o());
#endif
  }
  else {
    // B_bar = I implicitly!
    map_p_[p_bndy_idx] = map_p_bar_;
  }
  //
  // Create vectors
  //
  x0_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xL_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xU_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  p0_.resize(Np_);
  p0_[p_bndy_idx] = Teuchos::rcp(new Epetra_Vector(*map_p_[p_bndy_idx]));
  p0_[p_rx_idx] = Teuchos::rcp(new Epetra_Vector(*map_p_[p_rx_idx]));
  pL_.resize(Np_);
  pU_.resize(Np_);
  //pL_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //pU_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //
  // Initialize the vectors
  //
  x0_->PutScalar(x0);
  p0_[p_bndy_idx]->PutScalar(p0);
  p0_[p_rx_idx]->PutScalar(reactionRate);
  //
  // Initialize the graph for W
  //
  dat_->computeNpy(x0_);
  //if(dumpAll) { *out << "\nA =\n"; { Teuchos::OSTab tab(out); dat_->getA()->Print(*out); } }
  //if(dumpAll) { *out << "\nNpy =\n"; {  Teuchos::OSTab tab(out); dat_->getNpy()->Print(*out); } }
  W_graph_ = Teuchos::rcp(new Epetra_CrsGraph(dat_->getA()->Graph())); // Assume A and Npy have same graph!
  //
  // Get default objective matching vector q
  //
  q_ = Teuchos::rcp(new Epetra_Vector(*(*dat_->getq())(0))); // From Epetra_FEVector to Epetra_Vector!
  //
  isInitialized_ = true;
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  *out << "\nLeaving AdvDiffReactOptModel::AdvDiffReactOptModel(...) ...\n";
#endif
}

void AdvDiffReactOptModel::set_q( Teuchos::RefCountPtr<const Epetra_Vector> const& q )
{
  q_ = q;
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    dout = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab dtab(dout);
  *dout << "\nIn AdvDiffReactOptModel::set_q(q):  q =\n\n";
  q_->Print(Teuchos::OSTab(dout).o());
#endif
}

Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>
AdvDiffReactOptModel::getDataPool()
{
  return dat_;
}

Teuchos::RefCountPtr<const Epetra_MultiVector>
AdvDiffReactOptModel::get_B_bar() const
{
  return B_bar_;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(!(0<=l<=Np_));
  return map_p_[l];
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=0);
  return map_g_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(!(0<=l<=Np_));
  return p0_[l];
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_lower_bounds() const
{
  return xL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_upper_bounds() const
{
  return xU_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_lower_bounds(int l) const
{
  TEST_FOR_EXCEPT(!(0<=l<=Np_));
  return pL_[l];
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_upper_bounds(int l) const
{
  TEST_FOR_EXCEPT(!(0<=l<=Np_));
  return pU_[l];
}

Teuchos::RefCountPtr<Epetra_Operator>
AdvDiffReactOptModel::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

Teuchos::RefCountPtr<Epetra_Operator>
AdvDiffReactOptModel::create_DfDp_op(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,dat_->getB()->Graph()));
  // See DfDp in evalModel(...) below for details
}

EpetraExt::ModelEvaluator::InArgs
AdvDiffReactOptModel::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(2);
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
AdvDiffReactOptModel::createOutArgs() const
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
  if(supportDerivatives_) {
    outArgs.setSupports(
      OUT_ARG_DfDp,0
      ,( np_ > 0
         ? DerivativeSupport(DERIV_MV_BY_COL)
         : DerivativeSupport(DERIV_LINEAR_OP,DERIV_MV_BY_COL)
        )
      );
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
  }
  return outArgs;
}

void AdvDiffReactOptModel::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  Teuchos::TimeMonitor evalTimerMonitor(*evalTimer);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    dout = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab dtab(dout);
  *dout << "\nEntering AdvDiffReactOptModel::evalModel(...) ...\n";
#endif
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  const bool trace = ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_LOW) );
  const bool dumpAll = ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_EXTREME) );
  //
  Teuchos::OSTab tab(out);
  if(out.get() && trace) *out << "\n*** Entering AdvDiffReactOptModel::evalModel(...) ...\n"; 
  //
  // Get the input arguments
  //
  const Epetra_Vector *p_in = inArgs.get_p(p_bndy_idx).get();
  const Epetra_Vector &p = (p_in ? *p_in : *p0_[p_bndy_idx]);
  const Epetra_Vector *rx_vec_in = inArgs.get_p(p_rx_idx).get();
  const double        reactionRate = ( rx_vec_in ? (*rx_vec_in)[0] : (*p0_[p_rx_idx])[0] );
  const Epetra_Vector &x = *inArgs.get_x();
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *dout << "\nx =\n\n";
    x.Print(Teuchos::OSTab(dout).o());
    *dout << "\np =\n\n";
    p.Print(Teuchos::OSTab(dout).o());
#endif
  //
  // Get the output arguments
  //
  Epetra_Vector       *f_out = outArgs.get_f().get();
  Epetra_Vector       *g_out = outArgs.get_g(0).get();
  Epetra_Operator     *W_out = outArgs.get_W().get();
  Derivative          DfDp_out;
  Epetra_MultiVector  *DgDx_trans_out = 0;
  Epetra_MultiVector  *DgDp_trans_out = 0;
  if(supportDerivatives_) {
    DfDp_out = outArgs.get_DfDp(0);
    DgDx_trans_out = get_DgDx_mv(0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
    DgDp_trans_out = get_DgDp_mv(0,0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
  }
  //
  // Precompute some shared quantities
  //
  // p_bar = B_bar*p
  Teuchos::RefCountPtr<const Epetra_Vector> p_bar;
  if(np_ > 0) {
    Teuchos::TimeMonitor p_bar_TimerMonitor(*p_bar_Timer);
    Teuchos::RefCountPtr<Epetra_Vector> _p_bar;
    _p_bar = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
    _p_bar->Multiply('N','N',1.0,*B_bar_,p,0.0);
    p_bar = _p_bar;
  }
  else {
    p_bar = Teuchos::rcp(&p,false);
  }
  // R_p_bar = R*p_bar = R*(B_bar*p)
  Teuchos::RefCountPtr<const Epetra_Vector> R_p_bar;
  if( g_out || DgDp_trans_out ) {
    Teuchos::TimeMonitor R_p_bar_TimerMonitor(*R_p_bar_Timer);
      Teuchos::RefCountPtr<Epetra_Vector>
      _R_p_bar = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
    dat_->getR()->Multiply(false,*p_bar,*_R_p_bar);
    R_p_bar = _R_p_bar;
  }
  //
  // Compute the functions
  //
  if(f_out) {
    //
    // f = A*x + reationRate*Ny(x) + B*(B_bar*p)
    //
    Teuchos::TimeMonitor f_TimerMonitor(*f_Timer);
    Epetra_Vector &f = *f_out;
    Epetra_Vector Ax(*map_f_);
    dat_->getA()->Multiply(false,x,Ax);
    f.Update(1.0,Ax,0.0);
    if(reactionRate!=0.0) {
      dat_->computeNy(Teuchos::rcp(&x,false));
      f.Update(reactionRate,*dat_->getNy(),1.0);
    }
    Epetra_Vector Bp(*map_f_);
    dat_->getB()->Multiply(false,*p_bar,Bp);
    f.Update(1.0,Bp,-1.0,*dat_->getb(),1.0);
  }
  if(g_out) {
    //
    // g = 0.5 * (x-q)'*H*(x-q) + 0.5*regBeta*(B_bar*p)'*R*(B_bar*p)
    //
    Teuchos::TimeMonitor g_TimerMonitor(*g_Timer);
    Epetra_Vector &g = *g_out;
    Epetra_Vector xq(x);
    xq.Update(-1.0, *q_, 1.0);
    Epetra_Vector Hxq(x);
    dat_->getH()->Multiply(false,xq,Hxq);
    g[0] = 0.5*dot(xq,Hxq) + 0.5*dat_->getbeta()*dot(*p_bar,*R_p_bar);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *dout << "q =\n\n";
    q_->Print(Teuchos::OSTab(dout).o());
    *dout << "x =\n\n";
    x.Print(Teuchos::OSTab(dout).o());
    *dout << "x-q =\n\n";
    xq.Print(Teuchos::OSTab(dout).o());
    *dout << "H*(x-q) =\n\n";
    Hxq.Print(Teuchos::OSTab(dout).o());
    *dout << "0.5*dot(xq,Hxq) = " << (0.5*dot(xq,Hxq)) << "\n";
    *dout << "0.5*regBeta*(B_bar*p)'*R*(B_bar*p) = " << (0.5*dat_->getbeta()*dot(*p_bar,*R_p_bar)) << "\n";
#endif
  }
  if(W_out) {
    //
    // W = A + reationRate*Npy(x)
    //
    Teuchos::TimeMonitor W_TimerMonitor(*W_Timer);
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_out);
    if(reactionRate!=0.0)
      dat_->computeNpy(Teuchos::rcp(&x,false));
    Teuchos::RefCountPtr<Epetra_CrsMatrix>
      dat_A = dat_->getA(),
      dat_Npy = dat_->getNpy();
    const int numMyRows = dat_A->NumMyRows();
    for( int i = 0; i < numMyRows; ++i ) {
      int dat_A_num_row_entries=0; double *dat_A_row_vals=0; int *dat_A_row_inds=0;
      dat_A->ExtractMyRowView(i,dat_A_num_row_entries,dat_A_row_vals,dat_A_row_inds);
      int DfDx_num_row_entries=0; double *DfDx_row_vals=0; int *DfDx_row_inds=0;
      DfDx.ExtractMyRowView(i,DfDx_num_row_entries,DfDx_row_vals,DfDx_row_inds);
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPT(DfDx_num_row_entries!=dat_A_num_row_entries);
#endif
      if(reactionRate!=0.0) {
        int dat_Npy_num_row_entries=0; double *dat_Npy_row_vals=0; int *dat_Npy_row_inds=0;
        dat_Npy->ExtractMyRowView(i,dat_Npy_num_row_entries,dat_Npy_row_vals,dat_Npy_row_inds);
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPT(dat_A_num_row_entries!=dat_Npy_num_row_entries);
#endif
        for(int k = 0; k < DfDx_num_row_entries; ++k) {
#ifdef TEUCHOS_DEBUG
          TEST_FOR_EXCEPT(dat_A_row_inds[k]!=dat_Npy_row_inds[k]||dat_A_row_inds[k]!=DfDx_row_inds[k]);
#endif
          DfDx_row_vals[k] = dat_A_row_vals[k] + reactionRate * dat_Npy_row_vals[k];
        }
      }
      else {
        for(int k = 0; k < DfDx_num_row_entries; ++k) {
#ifdef TEUCHOS_DEBUG
          TEST_FOR_EXCEPT(dat_A_row_inds[k]!=DfDx_row_inds[k]);
#endif
          DfDx_row_vals[k] = dat_A_row_vals[k];
        }
      }
    }
  }
  if(!DfDp_out.isEmpty()) {
    if(out.get() && trace) *out << "\nComputing DfDp ...\n"; 
    //
    // DfDp = B*B_bar
    //
    Teuchos::TimeMonitor DfDp_TimerMonitor(*DfDp_Timer);
    Epetra_CrsMatrix   *DfDp_op = NULL;
    Epetra_MultiVector *DfDp_mv = NULL;
    if(out.get() && dumpAll)
    { *out << "\nB =\n"; { Teuchos::OSTab tab(out); dat_->getB()->Print(*out); } }
    if(DfDp_out.getLinearOp().get()) {
      DfDp_op = &dyn_cast<Epetra_CrsMatrix>(*DfDp_out.getLinearOp());
    }
    else {
      DfDp_mv = &*DfDp_out.getDerivativeMultiVector().getMultiVector();
      DfDp_mv->PutScalar(0.0);
    }
    Teuchos::RefCountPtr<Epetra_CrsMatrix>
      dat_B = dat_->getB();
    if(np_ > 0) {
      //
      // We only support a Multi-vector form when we have a non-idenity basis
      // matrix B_bar for p!
      //
      TEST_FOR_EXCEPT(DfDp_mv==NULL);
      dat_B->Multiply(false,*B_bar_,*DfDp_mv);
    }
    else {
      //
      // Note: We copy from B every time in order to be safe.  Note that since
      // the client knows that B is constant (sense we told them so in
      // createOutArgs()) then it should only compute this matrix once and keep
      // it if it is smart.
      //
      // Note: We support both the CrsMatrix and MultiVector form of this object
      // to make things easier for the client.
      //
      if(DfDp_op) {
        const int numMyRows = dat_B->NumMyRows();
        for( int i = 0; i < numMyRows; ++i ) {
          int dat_B_num_row_entries=0; double *dat_B_row_vals=0; int *dat_B_row_inds=0;
          dat_B->ExtractMyRowView(i,dat_B_num_row_entries,dat_B_row_vals,dat_B_row_inds);
          int DfDp_num_row_entries=0; double *DfDp_row_vals=0; int *DfDp_row_inds=0;
          DfDp_op->ExtractMyRowView(i,DfDp_num_row_entries,DfDp_row_vals,DfDp_row_inds);
#ifdef TEUCHOS_DEBUG
          TEST_FOR_EXCEPT(DfDp_num_row_entries!=dat_B_num_row_entries);
#endif
          for(int k = 0; k < DfDp_num_row_entries; ++k) {
#ifdef TEUCHOS_DEBUG
            TEST_FOR_EXCEPT(dat_B_row_inds[k]!=DfDp_row_inds[k]);
#endif
            DfDp_row_vals[k] = dat_B_row_vals[k];
          }
          // ToDo: The above code should be put in a utility function called copyValues(...)!
        }
      }
      else if(DfDp_mv) {
        // We must do a mat-vec to get this since we can't just copy out the
        // matrix entries since the domain map may be different from the
        // column map!  I learned this the very very hard way!  I am using
        // Thyra wrappers here since I can't figure out for the life of me how
        // to do this cleanly with Epetra alone!
        Teuchos::RefCountPtr<Epetra_Vector>
          etaVec = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
        Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> >
          space_p_bar = Thyra::create_VectorSpace(Teuchos::rcp(new Epetra_Map(*map_p_bar_)));
        Teuchos::RefCountPtr<Thyra::VectorBase<double> >
          thyra_etaVec = Thyra::create_Vector(etaVec,space_p_bar);
        for( int i = 0; i < map_p_bar_->NumGlobalElements(); ++i ) {
          Thyra::assign(thyra_etaVec.ptr(), 0.0);
          Thyra::set_ele(i, 1.0, thyra_etaVec.ptr());
          dat_B->Multiply(false, *etaVec, *(*DfDp_mv)(i));
        };
      }
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    if(DfDp_op) {
      *dout << "\nDfDp_op =\n\n";
      DfDp_op->Print(Teuchos::OSTab(dout).o());
    }
    if(DfDp_mv) {
      *dout << "\nDfDp_mv =\n\n";
      DfDp_mv->Print(Teuchos::OSTab(dout).o());
    }
#endif
  }
  if(DgDx_trans_out) {
    //
    // DgDx' = H*(x-q)
    //
    Teuchos::TimeMonitor DgDx_TimerMonitor(*DgDx_Timer);
    Epetra_Vector &DgDx_trans = *(*DgDx_trans_out)(0);
    Epetra_Vector xq(x);
    xq.Update(-1.0,*q_,1.0);
    dat_->getH()->Multiply(false,xq,DgDx_trans);
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *dout << "q =\n\n";
    q_->Print(Teuchos::OSTab(dout).o());
    *dout << "x =\n\n";
    x.Print(Teuchos::OSTab(dout).o());
    *dout << "x-q =\n\n";
    xq.Print(Teuchos::OSTab(dout).o());
    *dout << "DgDx_trans = H*(x-q) =\n\n";
    DgDx_trans.Print(Teuchos::OSTab(dout).o());
#endif
  }
  if(DgDp_trans_out) {
    //
    // DgDp' = regBeta*B_bar'*(R*(B_bar*p))
    //
    Teuchos::TimeMonitor DgDp_TimerMonitor(*DgDp_Timer);
    Epetra_Vector &DgDp_trans = *(*DgDp_trans_out)(0);
    if(np_ > 0) {
      DgDp_trans.Multiply('T','N',dat_->getbeta(),*B_bar_,*R_p_bar,0.0);
    }
    else {
      DgDp_trans.Update(dat_->getbeta(),*R_p_bar,0.0);
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *dout << "\nR_p_bar =\n\n";
    R_p_bar->Print(Teuchos::OSTab(dout).o());
    if(B_bar_.get()) {
      *dout << "\nB_bar =\n\n";
      B_bar_->Print(Teuchos::OSTab(dout).o());
    }
    *dout << "\nDgDp_trans =\n\n";
    DgDp_trans.Print(Teuchos::OSTab(dout).o());
#endif
  }
  if(out.get() && trace) *out << "\n*** Leaving AdvDiffReactOptModel::evalModel(...) ...\n"; 
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  *dout << "\nLeaving AdvDiffReactOptModel::evalModel(...) ...\n";
#endif
}

} // namespace GLpApp
