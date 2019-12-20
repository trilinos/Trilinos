// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_THYRARODUCTME_EQUALITYCONSTRAINT_SIMOPT
#define ROL_THYRARODUCTME_EQUALITYCONSTRAINT_SIMOPT

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"

#include "Tpetra_CrsMatrix.hpp"

using namespace ROL;

//! \brief ROL interface wrapper for Sacado SimOpt Constraint
template<class Real>
class ThyraProductME_Constraint_SimOpt : public Constraint_SimOpt<Real> {

public:

  ThyraProductME_Constraint_SimOpt(Thyra::ModelEvaluatorDefaultBase<double>& thyra_model_, int g_index_, const std::vector<int>& p_indices_,Teuchos::RCP<Teuchos::ParameterList> params_ = Teuchos::null, bool recompute = false) :
    thyra_model(thyra_model_), g_index(g_index_), p_indices(p_indices_), params(params_) {
    thyra_solver = Teuchos::null;
    updateValue = updateJacobian1 = true;
    value_ = 0;
    num_responses = -1;
    x_ptr = Teuchos::null;
    grad_ptr = Teuchos::null;
    if(params != Teuchos::null)
      params->set<int>("Optimizer Iteration Number", -1);
  };

  void setExternalSolver(Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> thyra_solver_) {
    thyra_solver = thyra_solver_;
  }

  void setNumResponses(int num_responses_) {
    num_responses = num_responses_;
  }

  void value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
    ThyraVector<Real>  & thyra_f = dynamic_cast<ThyraVector<Real>&>(c);
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    outArgs.set_f(thyra_f.getVector());
    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
    inArgs.set_x(thyra_x.getVector());

    thyra_model.evalModel(inArgs, outArgs);


    /*  {
          const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
          const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
          const ThyraVector<Real>  & thyra_f = dynamic_cast<const ThyraVector<Real>&>(c);


                Thyra::ConstDetachedVectorView<Real> x_view(thyra_x.getVector());
                    Thyra::ConstDetachedVectorView<Real> p_view(thyra_p.getVector());
                    Thyra::ConstDetachedVectorView<Real> f_view(thyra_f.getVector());

                    std::cout << "\nEnd of value... x:" << " ";
                    for (std::size_t i=0; i<x_view.subDim(); ++i)
                      std::cout << x_view(i) << " ";
                    std::cout << "\np:" << " ";
                    for (std::size_t i=0; i<p_view.subDim(); ++i)
                      std::cout << p_view(i) << " ";
                    std::cout << "\nf:" << " ";
                    for (std::size_t i=0; i<f_view.subDim(); ++i)
                      std::cout << f_view(i) << " ";

                    std::cout << "Norm: " << c.norm() <<std::endl;;
          }*/



    updateValue = false;
  }

  void applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
      const Vector<Real> &z, Real &tol) {

    if(updateJacobian1) {
      // Create Jacobian
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

      Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = thyra_model.get_W_factory();
      TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
      Teuchos::RCP< Thyra::LinearOpBase<double> > lop = thyra_model.create_W_op();

      for(std::size_t i=0; i<p_indices.size(); ++i)
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
      inArgs.set_x(thyra_x.getVector());

      outArgs.set_W_op(lop);
      thyra_model.evalModel(inArgs, outArgs);
      jac1 = lop;
      updateJacobian1 = false;
    }

    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_jv = dynamic_cast<ThyraVector<Real>&>(jv);
    jac1->apply(Thyra::NOTRANS, *thyra_v.getVector(), thyra_jv.getVector().ptr(), 1.0, 0.0);

  }

  void applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
      const Vector<Real> &z, Real &tol) {
    //std::cout <<  "Jacobian 2: " << tol <<std::endl;

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_jv = dynamic_cast<ThyraVector<Real>&>(jv);

    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_v = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_v.getVector());
    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    inArgs.set_x(thyra_x.getVector());
    for(std::size_t i=0; i<p_indices.size(); ++i) {
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

      // df/dp

      auto p_space = thyra_model.get_p_space(i);
      auto f_space = thyra_model.get_f_space();
      int num_params = p_space->dim();
      int num_resids = f_space->dim();
      auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Real>>(p_space);
      auto f_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Real>>(f_space);
      bool p_dist = !p_space_plus->isLocallyReplicated();//p_space->DistributedGlobal();
      bool f_dist = !f_space_plus->isLocallyReplicated();//f_space->DistributedGlobal();
      Thyra::ModelEvaluatorBase::DerivativeSupport ds =  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,i);
      // Determine which layout to use for df/dp.  Ideally one would look
      // at num_params, num_resids, what is supported by the underlying
      // model evaluator, and the sensitivity method, and make the best
      // choice to minimze the number of solves.  However this choice depends
      // also on what layout of dg/dx is supported (e.g., if only the operator
      // form is supported for forward sensitivities, then df/dp must be
      // DERIV_MV_BY_COL).  For simplicity, we order the conditional tests
      // to get the right layout in most situations.
      enum DerivativeLayout { OP, COL, ROW } dfdp_layout;
      { // if (sensitivity_method == "Adjoint")
        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP))
          dfdp_layout = OP;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !f_dist)
          dfdp_layout = ROW;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !p_dist)
          dfdp_layout = COL;
        else
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              std::endl << "Piro::NOXSolver::evalModel():  " <<
              "For df/dp(" << i <<") with adjoint sensitivities, " <<
              "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
              "DERIV_MV_BY_COL with p not distributed, or "
              "DERIV_TRANS_MV_BY_ROW with f not distributed." <<
              std::endl);
      }
      if (dfdp_layout == COL) {
        auto dfdp = Thyra::createMembers(f_space, num_params);
        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        outArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == ROW) {
        auto dfdp = Thyra::createMembers(p_space, num_resids);

        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        outArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == OP) {
        auto dfdp_op = thyra_model.create_DfDp_op(i);
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null, std::logic_error,
            std::endl << "Piro::NOXSolver::evalModel():  " <<
            "Needed df/dp operator (" << i << ") is null!" << std::endl);
        outArgs.set_DfDp(i,dfdp_op);
      }
    }

    thyra_model.evalModel(inArgs, outArgs);
    thyra_jv.zero();

    for(std::size_t i=0; i<p_indices.size(); ++i) {
      Thyra::ModelEvaluatorBase::Derivative<Real> dfdp_dv = outArgs.get_DfDp(i);
      auto dfdp_op = dfdp_dv.getLinearOp();
      auto dfdp = dfdp_dv.getMultiVector();

      if (dfdp_op != Teuchos::null) {
        auto temp_jv_ptr = Teuchos::rcp_dynamic_cast<ThyraVector<Real>>(thyra_jv.clone());
        temp_jv_ptr->zero();
        dfdp_op->apply(Thyra::NOTRANS,*thyra_prodvec_v->getVectorBlock(i), temp_jv_ptr->getVector().ptr(),1.0, 0.0);
        thyra_jv.axpy(1.0, *temp_jv_ptr);
      } else {
        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dfdp_orient =
            outArgs.get_DfDp(i).getMultiVectorOrientation();
        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient =
            outArgs.get_DgDx(0).getMultiVectorOrientation();

        Thyra::DetachedVectorView<Real> jv_view(thyra_jv.getVector());
        Thyra::ConstDetachedMultiVectorView<Real> v_view(thyra_v.getVector());
        Thyra::ConstDetachedMultiVectorView<Real> dfdp_view(dfdp);

        if (dgdx_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
          if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.subDim(); ix++)
                jv_view(ix) += v_view(ip,0)*dfdp_view(ip,ix);
          }
          else {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.numSubCols(); ix++)
                jv_view(ix) += v_view(ip,0)*dfdp_view(ix,ip);
          }
        }
        else {
          if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.subDim(); ix++)
                jv_view(ix) += v_view(0,ip)*dfdp_view(ip,ix);
          }
          else {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.numSubCols(); ix++)
                jv_view(ix) += v_view(0,ip)*dfdp_view(ix,ip);
          }
        }
      }
    }
  }


  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) {
    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_ijv = dynamic_cast<ThyraVector<Real>&>(ijv);

    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    Teuchos::RCP<Thyra::MultiVectorBase<Real> > thyra_ijv_ptr = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>(thyra_ijv.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

    inArgs.set_x(thyra_x.getVector());

    // Create Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = thyra_model.get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    Teuchos::RCP< Thyra::LinearOpBase<double> > lop;
    if(updateJacobian1)
      lop = thyra_model.create_W_op();
    else
      lop = jac1;


    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop));
    Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec;

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = thyra_model.create_W_prec();
    }
    const Teuchos::RCP<Thyra::LinearOpWithSolveBase<Real> > jacobian = lows_factory->createOp();

    if(updateJacobian1)
    {
      outArgs.set_W_op(lop);
      thyra_model.evalModel(inArgs, outArgs);
      outArgs.set_W_op(Teuchos::null);
      inArgs.set_x(Teuchos::null);
      jac1 = lop;
      updateJacobian1 = false;
    }

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      inArgs.set_x(thyra_x.getVector());
      outArgs.set_W_prec(prec);
      thyra_model.evalModel(inArgs, outArgs);
    }

    if(Teuchos::nonnull(prec))
      Thyra::initializePreconditionedOp<double>(*lows_factory,
          lop,
          prec,
          jacobian.ptr());
    else
      Thyra::initializeOp<double>(*lows_factory, lop, jacobian.ptr());

    const Thyra::SolveCriteria<Real> solve_criteria;

    Thyra::solve(
        *jacobian,
        Thyra::NOTRANS,
        *thyra_v.getVector(),
        thyra_ijv_ptr.ptr(),
        Teuchos::ptr(&solve_criteria));
  }



  void applyAdjointJacobian_1(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
      const Vector<Real> &z, Real &tol) {

    Teuchos::RCP< Thyra::LinearOpBase<double> > lop;
    if(updateJacobian1){
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);

      //Teuchos::RCP< Thyra::VectorBase<Real> > thyra_f = Thyra::createMember<Real>(thyra_model.get_f_space());
      Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i)
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
      inArgs.set_x(thyra_x.getVector());

      // Create implicitly transpose Jacobian and preconditioner
      Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = thyra_model.get_W_factory();
      TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));

      lop = thyra_model.create_W_op();
      outArgs.set_W_op(lop);
      thyra_model.evalModel(inArgs, outArgs);
      jac1 = lop;
      updateJacobian1 = false;
    }
    else
      lop = jac1;

    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_ajv = dynamic_cast<ThyraVector<Real>&>(ajv);
    lop->apply(Thyra::TRANS, *thyra_v.getVector(), thyra_ajv.getVector().ptr(), 1.0, 0.0);

  }

  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) {

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_iajv = dynamic_cast<ThyraVector<Real>&>(iajv);
    Teuchos::RCP<Thyra::MultiVectorBase<Real> > thyra_iajv_ptr = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>(thyra_iajv.getVector());
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
    inArgs.set_x(thyra_x.getVector());

    // Create implicitly transpose Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = thyra_model.get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    Teuchos::RCP< Thyra::LinearOpBase<double> > lop;

    if(updateJacobian1)
      lop = thyra_model.create_W_op();
    else
      lop = jac1;


    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop));
    Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec;

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = thyra_model.create_W_prec();
    }
    const Teuchos::RCP<Thyra::LinearOpWithSolveBase<Real> > jacobian = lows_factory->createOp();

    if(updateJacobian1)
    {
      outArgs.set_W_op(lop);
      thyra_model.evalModel(inArgs, outArgs);
      outArgs.set_W_op(Teuchos::null);
      jac1 = lop;
      updateJacobian1 = false;
    }

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      outArgs.set_W_prec(prec);
      thyra_model.evalModel(inArgs, outArgs);
    }

    if(Teuchos::nonnull(prec))
      Thyra::initializePreconditionedOp<double>(*lows_factory,
          Thyra::transpose<double>(lop),
          Thyra::unspecifiedPrec<double>(::Thyra::transpose<double>(prec->getUnspecifiedPrecOp())),
          jacobian.ptr());
    else
      Thyra::initializeOp<double>(*lows_factory, Thyra::transpose<double>(lop), jacobian.ptr());

    const Thyra::SolveCriteria<Real> solve_criteria;

    lop->apply(Thyra::NOTRANS, *thyra_v.getVector(), thyra_iajv.getVector().ptr(), 1.0, 0.0);

    thyra_iajv.getVector()->assign(0);
    Thyra::solve(
        *jacobian,
        Thyra::NOTRANS,
        *thyra_v.getVector(),
        thyra_iajv_ptr.ptr(),
        Teuchos::ptr(&solve_criteria));
  };

  void applyAdjointJacobian_2(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
      const Vector<Real> &z, Real &tol) {
    // std::cout << "Adjoint Jacobian 2" <<std::endl;

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
    const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);
    ThyraVector<Real>  & thyra_ajv = dynamic_cast<ThyraVector<Real>&>(ajv);
    //Teuchos::RCP< Thyra::VectorBase<Real> > thyra_f = Thyra::createMember<Real>(thyra_model.get_f_space());
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    Teuchos::RCP<Thyra::ProductVectorBase<Real> > thyra_prodvec_ajv = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Real>>(thyra_ajv.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    inArgs.set_x(thyra_x.getVector());
    for(std::size_t i=0; i<p_indices.size(); ++i) {
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

      // df/dp

      auto p_space = thyra_model.get_p_space(i);
      auto f_space = thyra_model.get_f_space();
      int num_params = p_space->dim();
      int num_resids = f_space->dim();
      auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Real>>(p_space);
      auto f_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Real>>(f_space);
      bool p_dist = !p_space_plus->isLocallyReplicated();//p_space->DistributedGlobal();
      bool f_dist = !f_space_plus->isLocallyReplicated();//f_space->DistributedGlobal();
      Thyra::ModelEvaluatorBase::DerivativeSupport ds =  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,i);
      // Determine which layout to use for df/dp.  Ideally one would look
      // at num_params, num_resids, what is supported by the underlying
      // model evaluator, and the sensitivity method, and make the best
      // choice to minimze the number of solves.  However this choice depends
      // also on what layout of dg/dx is supported (e.g., if only the operator
      // form is supported for forward sensitivities, then df/dp must be
      // DERIV_MV_BY_COL).  For simplicity, we order the conditional tests
      // to get the right layout in most situations.
      enum DerivativeLayout { OP, COL, ROW } dfdp_layout;
      { // if (sensitivity_method == "Adjoint")
        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP))
          dfdp_layout = OP;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !f_dist)
          dfdp_layout = ROW;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !p_dist)
          dfdp_layout = COL;
        else
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              std::endl << "Piro::NOXSolver::evalModel():  " <<
              "For df/dp(" << i <<") with adjoint sensitivities, " <<
              "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
              "DERIV_MV_BY_COL with p not distributed, or "
              "DERIV_TRANS_MV_BY_ROW with f not distributed." <<
              std::endl);
      }
      if (dfdp_layout == COL) {
        auto dfdp = Thyra::createMembers(f_space, num_params);
        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        outArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == ROW) {
        auto dfdp = Thyra::createMembers(p_space, num_resids);

        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        outArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == OP) {
        auto dfdp_op = thyra_model.create_DfDp_op(i);
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null, std::logic_error,
            std::endl << "Piro::NOXSolver::evalModel():  " <<
            "Needed df/dp operator (" << i << ") is null!" << std::endl);
        outArgs.set_DfDp(i,dfdp_op);
      }
    }

    thyra_model.evalModel(inArgs, outArgs);

    for(std::size_t i=0; i<p_indices.size(); ++i) {
      Thyra::ModelEvaluatorBase::Derivative<Real> dfdp_dv = outArgs.get_DfDp(i);
      auto dfdp_op = dfdp_dv.getLinearOp();
      auto dfdp = dfdp_dv.getMultiVector();
      // auto tmp = Thyra::createMembers(dfdp_op->domain(), thyra_v.getVector()->domain()->dim());

      if (dfdp_op != Teuchos::null) {
        dfdp_op->apply(Thyra::TRANS,*thyra_v.getVector(), thyra_prodvec_ajv->getNonconstVectorBlock(i).ptr(),1.0, 0.0);
        // Thyra::update(1.0,  *tmp, thyra_ajv.getMultiVector().ptr());
      } else {
        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dfdp_orient =
            outArgs.get_DfDp(i).getMultiVectorOrientation();
        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient =
            outArgs.get_DgDx(0).getMultiVectorOrientation();

        Thyra::DetachedVectorView<Real> ajv_view(thyra_ajv.getVector());
        Thyra::ConstDetachedMultiVectorView<Real> v_view(thyra_v.getVector());
        Thyra::ConstDetachedMultiVectorView<Real> dfdp_view(dfdp);
        if (dgdx_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
          if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.subDim(); ix++)
                ajv_view(ip) += v_view(ix,0)*dfdp_view(ip,ix);
          }
          else {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.numSubCols(); ix++)
                ajv_view(ip) += v_view(ix,0)*dfdp_view(ix,ip);
          }
        }
        else {
          if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.subDim(); ix++)
                ajv_view(ip) += v_view(0,ix)*dfdp_view(ip,ix);
          }
          else {
            for (std::size_t ip=0; ip<p_indices.size(); ++ip)
              for (int ix=0; ix<dfdp_view.numSubCols(); ix++)
                ajv_view(ip) += v_view(0,ix)*dfdp_view(ix,ip);
          }
        }
      }
    }
  }

  void solve(Vector<Real> &c,
      Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) {

    if(thyra_solver.is_null())
      Constraint_SimOpt<Real>::solve(c,u,z,tol);
    else {
      if(this->zero_) { //SimOpt.Solve.Zero Initial Guess
        u.zero();
        update_1(u);
      }

      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      ThyraVector<Real>  & thyra_x = dynamic_cast<ThyraVector<Real>&>(u);
      ThyraVector<Real>  & thyra_f = dynamic_cast<ThyraVector<Real>&>(c);
      Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > gx = Thyra::createMember<Real>(thyra_solver->get_g_space(num_responses));

      Teuchos::ParameterList origTestList;

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_solver->createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_solver->createOutArgs();

      outArgs.set_f(thyra_f.getVector());
      outArgs.set_g(num_responses, gx);
      for(std::size_t i=0; i<p_indices.size(); ++i)
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

      inArgs.set_x(thyra_x.getVector());

      params->set<bool>("Compute State", true);
      thyra_solver->evalModel(inArgs, outArgs);

      Teuchos::RCP<const Thyra::VectorBase<double> > gx_out = outArgs.get_g(num_responses);
      if (Teuchos::nonnull(gx_out)) {
        Thyra::copy(*gx_out, thyra_x.getVector().ptr());
      }
      this->update_1(u);
    }

    updateValue = false;
  }


  /** \brief Update constraint functions with respect to Sim variable.
                x is the optimization variable,
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
   */
  void update_1( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    if (flag == true) {
      updateValue = true;
      updateJacobian1 = true;
    }
    if(params != Teuchos::null)
      params->set<int>("Optimizer Iteration Number", iter);
  }

  /** \brief Update constraint functions with respect to Opt variable.
                x is the optimization variable,
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
   */
  void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    if (flag == true) {
      updateValue = true;
      updateJacobian1 = true;
    }
    if(Teuchos::nonnull(params)) {
      params->set<int>("Optimizer Iteration Number", iter);
      if(flag == true)
        params->set<bool>("Optimization Variables Changed",true);
    }
  }

public:
  bool updateValue, updateJacobian1;

private:
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> thyra_solver;
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  const int g_index;
  const std::vector<int> p_indices;
  int num_responses;
  Real value_;
  Teuchos::RCP<Vector<Real> > x_ptr, grad_ptr;
  Teuchos::RCP<Teuchos::ParameterList> params;
  Teuchos::RCP< Thyra::LinearOpBase<double> > jac1;

};


#endif
