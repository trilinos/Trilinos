// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_PRODUCTMODELEVAL_HPP
#define PIRO_PRODUCTMODELEVAL_HPP

#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp_decl.hpp"

#ifdef HAVE_PIRO_ROL
#include "ROL_Types.hpp"
#endif

#ifdef HAVE_PIRO_TEKO
#include "Teko_InverseLibrary.hpp"
#include "Teko_PreconditionerFactory.hpp"
#endif

namespace Piro {

/** \brief Product Model Evaluator
 *  Model Evaluator that supports only one parameter which is
 *  a product vector. The evaluator has an internal evaluator that
 *  treats this product parameter as multiple parameters.
 */

template<class Real>
class ProductModelEvaluator : public Thyra::ModelEvaluatorDelegatorBase<Real>
{
public:

    ProductModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model,
                            const std::vector<int>& p_indices);

    ~ProductModelEvaluator();

    /** \brief . */
    int Np() const;
    /** \brief . */
    int Ng() const;

    Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> get_x_space() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> get_f_space() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> get_g_space( int l) const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> get_p_space( int l) const;

    Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;

    Teuchos::RCP<Thyra::LinearOpBase<Real> > create_W_op() const;
    Teuchos::RCP<Thyra::PreconditionerBase<Real> > create_W_prec() const;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real> > get_W_factory() const;

#ifdef HAVE_PIRO_ROL
    void block_diagonal_hessian_22(const Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> H,
                    const ROL::Vector<Real> &u,
                    const ROL::Vector<Real> &z,
                    const int g_idx) const;
#endif

    /** \brief . */
    Teuchos::RCP<Thyra::LinearOpBase<Real> > create_DfDp_op(int l) const;
    /** \brief . */
    Teuchos::RCP<Thyra::LinearOpBase<Real> > create_DgDp_op(int j, int l) const;
    /** \brief . */
    Teuchos::RCP<Thyra::LinearOpBase<Real> > create_DgDp_op(int j, int l, Teuchos::RCP<Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp) const;

    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> getModel() { return thyra_model_; }

    Thyra::ModelEvaluatorBase::InArgs<Real>  createInArgs() const;

    void reportFinalPoint(
        const Thyra::ModelEvaluatorBase::InArgs<Real>& finalPoint,
        const bool wasSolved);
        
    Teuchos::ArrayView<const std::string> get_g_names(int j) const;

    /** \brief . */
    ::Thyra::ModelEvaluatorBase::InArgs<Real> getNominalValues() const;
    /** \brief . */
    ::Thyra::ModelEvaluatorBase::InArgs<Real> getLowerBounds() const;
    /** \brief . */
    ::Thyra::ModelEvaluatorBase::InArgs<Real> getUpperBounds() const;

protected:

    /** \brief . */
    Thyra::ModelEvaluatorBase::OutArgs<Real>
    createOutArgsImpl() const;

    /** \brief . */
    void
    evalModelImpl(
        const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs,
        const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs) const;
    //@}


private:

    void fromInternalInArgs(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs1, Thyra::ModelEvaluatorBase::InArgsSetup<Real>& inArgs2) const;
    void toInternalInArgs(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs1, Thyra::ModelEvaluatorBase::InArgsSetup<Real>& inArgs2) const;
    void fromInternalOutArgs(const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs1, Thyra::ModelEvaluatorBase::OutArgsSetup<Real>& outArgs2) const;
    void toInternalOutArgs(const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs1, Thyra::ModelEvaluatorBase::OutArgsSetup<Real>& outArgs2) const;

    /** \brief . */
    Thyra::ModelEvaluatorBase::InArgs<Real>  createInArgsImpl() const;

    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model_;
    const std::vector<int> p_indices_;
    Teuchos::Array<Thyra::ModelEvaluatorBase::DerivativeSupport> DfDp_op_support_;
    Teuchos::Array<Thyra::ModelEvaluatorBase::DerivativeSupport> DgDp_op_support_;
}; // class ProductModelEvaluator


template <typename Real>
ProductModelEvaluator<Real>::
ProductModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model,
    const std::vector<int>& p_indices) :
    Thyra::ModelEvaluatorDelegatorBase<Real>(thyra_model),
    thyra_model_(thyra_model),
    p_indices_(p_indices)
{
    DfDp_op_support_.clear();
    DgDp_op_support_.clear();
    Thyra::ModelEvaluatorBase::OutArgs<Real> internal_outArgs = thyra_model_->createOutArgs();
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        DfDp_op_support_.push_back(internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, p_indices_[i]));
    }
    for (auto i = 0; i < thyra_model_->Ng(); ++i) {
        for (std::size_t j = 0; j < p_indices_.size(); ++j) {
            DgDp_op_support_.push_back(internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, p_indices_[j]));
        }
    }
}

template <typename Real>
ProductModelEvaluator<Real>::~ProductModelEvaluator()
{
}

template <typename Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>>
ProductModelEvaluator<Real>::get_x_space() const
{
    return thyra_model_->get_x_space();
}

template <typename Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>>
ProductModelEvaluator<Real>::get_f_space() const
{
    return thyra_model_->get_f_space();
}

template <typename Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>>
ProductModelEvaluator<Real>::get_p_space(int l) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                        "Error!  ProductModelEvaluator<Real>::get_p_space() " <<
                        "supports only 1 parameter vector.  Supplied index l = " <<
                        l << std::endl);

    Teuchos::Array<Teuchos::RCP<Thyra::VectorSpaceBase<Real> const>> p_spaces(p_indices_.size());
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        p_spaces[i] = thyra_model_->get_p_space(p_indices_[i]);
    }
    Teuchos::RCP<Thyra::DefaultProductVectorSpace<Real> const> p_space = Thyra::productVectorSpace<Real>(p_spaces);

    return p_space;
}

template <typename Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>>
ProductModelEvaluator<Real>::get_g_space(int l) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(l > thyra_model_->Ng(), std::logic_error,
                        "Error!  ProductModelEvaluator::get_g_space() Supplied index l = " <<
                        l << " is greater or equal to the number of responses of the underlying model " << 
                        thyra_model_->Ng() << std::endl);
    Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> g_space = thyra_model_->get_g_space(l);
    return g_space;
}

template <typename Real>
Teuchos::RCP<const  Teuchos::Array<std::string> >
ProductModelEvaluator<Real>::get_p_names(int l) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                        "Error!  ProductModelEvaluator<Real>::get_p_names() " <<
                        "supports only 1 parameter vector.  Supplied index l = " <<
                        l << std::endl);

    Teuchos::RCP<Teuchos::Array<std::string> > p_names =
        Teuchos::rcp(new Teuchos::Array<std::string>(p_indices_.size()) );
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
    std::stringstream ss;
    ss << "Parameter " << i;
    const std::string name = ss.str();
    (*p_names)[i] = name;
    }
    return thyra_model_->get_p_names(l);
}

template <typename Real>
Teuchos::RCP<Thyra::LinearOpBase<Real>>
ProductModelEvaluator<Real>::create_W_op() const
{
    return thyra_model_->create_W_op();
}

template <typename Real>
Teuchos::RCP<Thyra::PreconditionerBase<Real>>
ProductModelEvaluator<Real>::create_W_prec() const
{
    return thyra_model_->create_W_prec();
}

template <typename Real>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real>>
ProductModelEvaluator<Real>::get_W_factory() const
{
    return thyra_model_->get_W_factory();
}

template <typename Real>
Thyra::ModelEvaluatorBase::InArgs<Real>
ProductModelEvaluator<Real>::createInArgs() const
{
    Thyra::ModelEvaluatorBase::InArgs<Real> internal_inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalInArgs(internal_inArgs, result);

    return result;
}

template <typename Real>
Thyra::ModelEvaluatorBase::OutArgs<Real>
ProductModelEvaluator<Real>::createOutArgsImpl() const
{
    Thyra::ModelEvaluatorBase::OutArgs<Real> internal_outArgs = thyra_model_->createOutArgs();
    Thyra::ModelEvaluatorBase::OutArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalOutArgs(internal_outArgs, result);

    return result;
}

template <typename Real>
void 
ProductModelEvaluator<Real>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs) const
{
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> internal_inArgs;
    Thyra::ModelEvaluatorBase::OutArgsSetup<Real> internal_outArgs;

    internal_outArgs.setModelEvalDescription(outArgs.modelEvalDescription()+"_internal");

    internal_inArgs.set_Np_Ng(thyra_model_->Np(), thyra_model_->Ng());
    internal_outArgs.set_Np_Ng(thyra_model_->Np(), thyra_model_->Ng());

    bool supports_dfdp_op = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,0).supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    
    bool supports_vec_prod_f_xp = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, 0);
    bool supports_vec_prod_f_px = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, 0);
    bool supports_vec_prod_f_pp = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, 0, 0);

    this->toInternalInArgs(inArgs, internal_inArgs);
    this->toInternalOutArgs(outArgs, internal_outArgs);

    internal_outArgs.setArgs(outArgs, true);
    internal_inArgs.setArgs(inArgs, true);

    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(inArgs.get_p(0));
    Teuchos::RCP<const Thyra::ProductMultiVectorBase<Real> > prodvec_direction_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<Real>>(inArgs.get_p_direction(0));

    internal_inArgs.set_x_direction(inArgs.get_x_direction());

    TEUCHOS_TEST_FOR_EXCEPTION(!inArgs.get_p(0).is_null() && prodvec_p.is_null(), std::logic_error,
        std::endl <<
        "Error!  ProductModelEvaluator<Real>::evalModelImpl() " <<
        " prodvec_p is not a ProductVectorBase " << std::endl);

    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        if (!prodvec_p.is_null()) {
            auto tmp = prodvec_p->getVectorBlock(i);

            Teuchos::RCP<const Thyra::ProductVectorBase<Real> > prodvec_p_in
                = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(tmp);

            TEUCHOS_TEST_FOR_EXCEPTION(!prodvec_p_in.is_null(), std::logic_error,
                std::endl <<
                "Error!  ProductModelEvaluator<Real>::evalModelImpl() " <<
                " ProductVectorBase of ProductVectorBase is not supported.  Parameter index i = " <<
                i << std::endl);

            internal_inArgs.set_p(p_indices_[i], prodvec_p->getVectorBlock(i));
        }
        if (!prodvec_direction_p.is_null()) {
            internal_inArgs.set_p_direction(p_indices_[i], prodvec_direction_p->getMultiVectorBlock(i));
        }
    }

    for (auto g_index = 0; g_index < thyra_model_->Ng(); ++g_index) {
        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, 0)) {
            std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > hv_vec(p_indices_.size());

            hv_vec[0] = outArgs.get_hess_vec_prod_g_xp(g_index,0);
            if (!Teuchos::is_null(hv_vec[0])) {
                for(std::size_t j=1; j<p_indices_.size(); ++j) {
                    hv_vec[j] = hv_vec[0]->clone_mv();
                }

                for(std::size_t j=0; j<p_indices_.size(); ++j) {
                    internal_outArgs.set_hess_vec_prod_g_xp(g_index,p_indices_[j], hv_vec[j]);
                }
            }
        }

        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, 0)) {
            Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_hv =
                Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(outArgs.get_hess_vec_prod_g_px(g_index,0));

            if (!Teuchos::is_null(prodvec_hv)) {
                for(std::size_t i=0; i<p_indices_.size(); ++i) {
                    bool supports_deriv_j =   internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, p_indices_[i]);
                    TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv_j, std::logic_error, 
                        "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_g_px is not supported. " <<
                        "Solution index = " << g_index << " Parameter index i = " << i << std::endl);
                    internal_outArgs.set_hess_vec_prod_g_px(g_index,p_indices_[i], prodvec_hv->getNonconstMultiVectorBlock(i));
                }
            }
            else {
                TEUCHOS_TEST_FOR_EXCEPTION( !Teuchos::is_null(outArgs.get_hess_vec_prod_g_px(g_index,0)), std::logic_error, "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_g_px is not a ProductMultiVectorBase. Solution index = " << g_index << std::endl);
            }
        }

        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, 0, 0)) {
            Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_hv =
                Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(outArgs.get_hess_vec_prod_g_pp(g_index,0,0));
            std::vector<std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > > hv_vec(p_indices_.size());

            if (!Teuchos::is_null(prodvec_hv)) {
                for(std::size_t i=0; i<p_indices_.size(); ++i) {
                    hv_vec[i].resize(p_indices_.size());
                    hv_vec[i][0] = prodvec_hv->getNonconstMultiVectorBlock(i);
                    for(std::size_t j=1; j<p_indices_.size(); ++j) {
                        hv_vec[i][j] = hv_vec[i][0]->clone_mv();
                    }
                }

                for(std::size_t i=0; i<p_indices_.size(); ++i) {
                    for(std::size_t j=0; j<p_indices_.size(); ++j) {
                        bool supports_deriv_j =   internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_indices_[i], p_indices_[j]);
                        TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv_j, std::logic_error, 
                            "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_g_pp is not supported. " <<
                            "Solution index = " << g_index << " Parameter index i = " << i << " Parameter index j = " << j << std::endl);

                        internal_outArgs.set_hess_vec_prod_g_pp(g_index,p_indices_[i], p_indices_[j], hv_vec[i][j]);
                    }
                }
            }
            else {
                TEUCHOS_TEST_FOR_EXCEPTION( !Teuchos::is_null(outArgs.get_hess_vec_prod_g_pp(g_index,0,0)), std::logic_error, 
                    "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_g_pp is not a ProductMultiVectorBase. " <<
                    "Solution index = " << g_index << std::endl);
            }
        }
    }

    if (supports_vec_prod_f_xp) {
        Teuchos::RCP< Thyra::MultiVectorBase<Real> > thyra_ahwv = outArgs.get_hess_vec_prod_f_xp(0);
        std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > ahwv_vec(p_indices_.size());

        if (!Teuchos::is_null(thyra_ahwv)) {
            ahwv_vec[0] = thyra_ahwv;
            for(std::size_t j=1; j<p_indices_.size(); ++j) {
                ahwv_vec[j] = thyra_ahwv->clone_mv();
            }

            for(std::size_t j=0; j<p_indices_.size(); ++j) {
                bool supports_deriv_j =   internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, p_indices_[j]);
                TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv_j, std::logic_error, 
                    "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_f_xp is not supported. Parameter index j = " << j << std::endl);
                internal_outArgs.set_hess_vec_prod_f_xp(p_indices_[j], ahwv_vec[j]);
            }
        }
    }

    if (supports_vec_prod_f_px) {
        Teuchos::RCP< Thyra::ProductVectorBase<Real> > prodvec_ahwv =
            Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Real>>(outArgs.get_hess_vec_prod_f_px(0));
        if (!Teuchos::is_null(prodvec_ahwv)) {
            for(std::size_t i=0; i<p_indices_.size(); ++i) {
                bool supports_deriv_i =   internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, p_indices_[i]);
                TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv_i, std::logic_error, 
                    "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_f_px is not supported. Parameter index i = " << i << std::endl);
                internal_outArgs.set_hess_vec_prod_f_px(p_indices_[i], prodvec_ahwv->getNonconstVectorBlock(i));
            }
        }
        else {
            TEUCHOS_TEST_FOR_EXCEPTION( !Teuchos::is_null(outArgs.get_hess_vec_prod_f_px(0)), std::logic_error, 
                "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_f_px is not a ProductMultiVectorBase");
        }
    }

    if (supports_vec_prod_f_pp) {
        Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_ahwv =
            Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(outArgs.get_hess_vec_prod_f_pp(0,0));
        std::vector<std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > > ahwv_vec(p_indices_.size());

        if (!Teuchos::is_null(prodvec_ahwv)) {
            for(std::size_t i=0; i<p_indices_.size(); ++i) {
                ahwv_vec[i].resize(p_indices_.size());
                ahwv_vec[i][0] = prodvec_ahwv->getNonconstMultiVectorBlock(i);
                for(std::size_t j=1; j<p_indices_.size(); ++j) {
                    ahwv_vec[i][j] = ahwv_vec[i][0]->clone_mv();
                }
            }

            for(std::size_t i=0; i<p_indices_.size(); ++i) {
                for(std::size_t j=0; j<p_indices_.size(); ++j) {
                    bool supports_deriv_ij =   internal_outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, p_indices_[i], p_indices_[j]);
                    TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv_ij, std::logic_error, 
                        "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_f_pp is not supported. Parameter index i = " << i << " Parameter index j = " << j << std::endl);

                    internal_outArgs.set_hess_vec_prod_f_pp(p_indices_[i], p_indices_[j], ahwv_vec[i][j]);
                }
            }
        }
        else {
            TEUCHOS_TEST_FOR_EXCEPTION( !Teuchos::is_null(outArgs.get_hess_vec_prod_f_pp(0,0)), std::logic_error, 
                "ProductModelEvaluator<Real>::evalModelImpl(): hess_vec_prod_f_pp is not a ProductMultiVectorBase");
        }
    }

    if (supports_dfdp_op) {
        Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dfdp_op =
            Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(outArgs.get_DfDp(0).getLinearOp());
        if (Teuchos::nonnull(dfdp_op)) {
            for(std::size_t j=0; j<p_indices_.size(); ++j)
                internal_outArgs.set_DfDp( p_indices_[j], Thyra::ModelEvaluatorBase::Derivative<Real>(Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real>>(dfdp_op->getNonconstBlock(0, j))));
        }
    }

    for (auto i = 0; i < thyra_model_->Ng(); ++i) {
        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,0).supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
            Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dgdp_op =
                Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(outArgs.get_DgDp(i,0).getLinearOp());
            if (Teuchos::nonnull(dgdp_op)) {
                for(std::size_t j=0; j<p_indices_.size(); ++j) {
                    auto dgdp_j_mv =
                        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>( Teuchos::rcp_dynamic_cast<Thyra::DefaultScaledAdjointLinearOp<Real>>(dgdp_op->getNonconstBlock(0, j))->getNonconstOp() );
                    if (!Teuchos::is_null(dgdp_j_mv)) {
                        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
                            DgDp_op_support_[i*p_indices_.size()+j];
                        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
                        if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
                            dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
                        }
                        else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
                            dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
                        }
                        else {
                            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "ProductModelEvaluator<Real>::evalModelImpl(): DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
                        }
                        internal_outArgs.set_DgDp(i, p_indices_[j], Thyra::ModelEvaluatorBase::Derivative<Real>(dgdp_j_mv, dgdp_orient));
                    }
                }
            }
        }
    }

    for (auto i = 0; i < thyra_model_->Ng(); ++i) {
        bool suppGradMV = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,0).supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        bool suppJacMV = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,0).supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        if (suppGradMV || suppJacMV) {
            auto dgdp_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(outArgs.get_DgDp(i,0).getMultiVector());
            if (Teuchos::nonnull(dgdp_vec)) {
                auto dgdp_orient = suppGradMV ? Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM : Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
                for(std::size_t j=0; j<p_indices_.size(); ++j) {
                    internal_outArgs.set_DgDp(i, p_indices_[j], Thyra::ModelEvaluatorBase::Derivative<Real>(dgdp_vec->getNonconstMultiVectorBlock(j), dgdp_orient));
               }
            }
        }
    }

    thyra_model_->evalModel(internal_inArgs,internal_outArgs);

    for (auto g_index = 0; g_index < thyra_model_->Ng(); ++g_index) {
        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, 0)) {
            Teuchos::RCP< Thyra::MultiVectorBase<Real> > hv_vec = internal_outArgs.get_hess_vec_prod_g_xp(g_index,p_indices_[0]);
            if (!Teuchos::is_null(hv_vec)) {
                for(std::size_t j=1; j<p_indices_.size(); ++j) {
                    Teuchos::RCP< Thyra::MultiVectorBase<Real> > hv_vec_tmp = internal_outArgs.get_hess_vec_prod_g_xp(g_index,p_indices_[j]);
                    if (!Teuchos::is_null(hv_vec_tmp)) {
                        hv_vec->update(1.0, *hv_vec_tmp);
                    }
                }
            }
        }

        if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, 0, 0)) {
            for(std::size_t i=0; i<p_indices_.size(); ++i) {
                Teuchos::RCP< Thyra::MultiVectorBase<Real> > hv_vec = internal_outArgs.get_hess_vec_prod_g_pp(g_index,p_indices_[i],p_indices_[0]);
                if (!Teuchos::is_null(hv_vec)) {
                    for(std::size_t j=1; j<p_indices_.size(); ++j) {
                        Teuchos::RCP< Thyra::MultiVectorBase<Real> > hv_vec_tmp = internal_outArgs.get_hess_vec_prod_g_pp(g_index,p_indices_[i],p_indices_[j]);
                        if (!Teuchos::is_null(hv_vec_tmp))
                            hv_vec->update(1.0, *hv_vec_tmp);
                    }
                }
            }
        }
    }

    if (supports_vec_prod_f_xp) {
        Teuchos::RCP< Thyra::MultiVectorBase<Real> > ahwv_vec = internal_outArgs.get_hess_vec_prod_f_xp(p_indices_[0]);
        if (!Teuchos::is_null(ahwv_vec)) {
            for(std::size_t j=1; j<p_indices_.size(); ++j) {
                Teuchos::RCP< Thyra::MultiVectorBase<Real> > ahwv_vec_tmp = internal_outArgs.get_hess_vec_prod_f_xp(p_indices_[j]);
                if (!Teuchos::is_null(ahwv_vec_tmp))
                    ahwv_vec->update(1.0, *ahwv_vec_tmp);
            }
        }
    }

    if (supports_vec_prod_f_pp) {
        for(std::size_t i=0; i<p_indices_.size(); ++i) {
            Teuchos::RCP< Thyra::MultiVectorBase<Real> > ahwv_vec = internal_outArgs.get_hess_vec_prod_f_pp(p_indices_[i],p_indices_[0]);
            if (!Teuchos::is_null(ahwv_vec)) {
                for(std::size_t j=1; j<p_indices_.size(); ++j) {
                    Teuchos::RCP< Thyra::MultiVectorBase<Real> > ahwv_vec_tmp = internal_outArgs.get_hess_vec_prod_f_pp(p_indices_[i],p_indices_[j]);
                    if (!Teuchos::is_null(ahwv_vec_tmp))
                        ahwv_vec->update(1.0, *ahwv_vec_tmp);
                }
            }
        }
    }
}

template <typename Real>
Thyra::ModelEvaluatorBase::InArgs<Real>
ProductModelEvaluator<Real>::createInArgsImpl() const
{
    Thyra::ModelEvaluatorBase::InArgs<Real> internal_inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalInArgs(internal_inArgs, result);

    return result;
}

template <typename Real>
Teuchos::ArrayView<const std::string>
ProductModelEvaluator<Real>::get_g_names(int j) const
{
    return thyra_model_->get_g_names(j);
}

template <typename Real>
Thyra::ModelEvaluatorBase::InArgs<Real>
ProductModelEvaluator<Real>::getNominalValues() const
{
    Thyra::ModelEvaluatorBase::InArgs<Real> internal_inArgs = thyra_model_->getNominalValues();
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalInArgs(internal_inArgs, result);
    result.setArgs(internal_inArgs, true);

    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Real>> p_space = Teuchos::rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<Real>>(this->get_p_space(0));

    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> p_vecs(p_indices_.size());
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        p_vecs[i] = internal_inArgs.get_p(p_indices_[i]);
    }
    Teuchos::RCP<Thyra::DefaultProductVector<Real>> p_prod = Thyra::defaultProductVector<double>(p_space, p_vecs());

    result.set_p(0, p_prod);

    return result;
}

template <typename Real>
Thyra::ModelEvaluatorBase::InArgs<Real>
ProductModelEvaluator<Real>::getLowerBounds() const
{
    Thyra::ModelEvaluatorBase::InArgs<Real> internal_inArgs = thyra_model_->getLowerBounds();
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalInArgs(internal_inArgs, result);
    result.setArgs(internal_inArgs, true);

    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Real>> p_space = Teuchos::rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<Real>>(this->get_p_space(0));

    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> p_vecs(p_indices_.size());
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        p_vecs[i] = internal_inArgs.get_p(p_indices_[i]);
    }
    Teuchos::RCP<Thyra::DefaultProductVector<Real>> p_prod = Thyra::defaultProductVector<double>(p_space, p_vecs());

    result.set_p(0, p_prod);

    return result;
}

template <typename Real>
Thyra::ModelEvaluatorBase::InArgs<Real>
ProductModelEvaluator<Real>::getUpperBounds() const
{
    Thyra::ModelEvaluatorBase::InArgs<Real> internal_inArgs = thyra_model_->getUpperBounds();
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> result; 
    result.setModelEvalDescription(this->description());
    result.set_Np_Ng(1, thyra_model_->Ng());

    this->fromInternalInArgs(internal_inArgs, result);
    result.setArgs(internal_inArgs, true);

    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Real>> p_space = Teuchos::rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<Real>>(this->get_p_space(0));

    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> p_vecs(p_indices_.size());
    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        p_vecs[i] = internal_inArgs.get_p(p_indices_[i]);
    }
    Teuchos::RCP<Thyra::DefaultProductVector<Real>> p_prod = Thyra::defaultProductVector<double>(p_space, p_vecs());

    result.set_p(0, p_prod);

    return result;
}

template <typename Real>
int
ProductModelEvaluator<Real>::Np() const
{
    return 1;
}

template <typename Real>
int
ProductModelEvaluator<Real>::Ng() const
{
    return thyra_model_->Ng();
}

template <typename Real>
void
ProductModelEvaluator<Real>::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<Real>& finalPoint,
    const bool wasSolved)
{
    return thyra_model_->reportFinalPoint(finalPoint, wasSolved);
}

template <typename Real>
Teuchos::RCP<Thyra::LinearOpBase<Real> > 
ProductModelEvaluator<Real>::create_DfDp_op(int l) const {
    
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> J = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<Real>());
    J->beginBlockFill(1, p_indices_.size());
    for(std::size_t i=0; i<p_indices_.size(); ++i) {
        Teuchos::RCP<Thyra::LinearOpBase<Real> > dfdp_op = thyra_model_->create_DfDp_op(p_indices_[i]);
        J->setNonconstBlock(0, i, dfdp_op);
    }
    J->endBlockFill();
    return J;
}

template <typename Real>
Teuchos::RCP<Thyra::LinearOpBase<Real> > 
ProductModelEvaluator<Real>::create_DgDp_op(int j, int l) const {    
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> J = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<Real>());
    J->beginBlockFill(1, p_indices_.size());
    for(std::size_t i=0; i<p_indices_.size(); ++i) {
        Teuchos::RCP<Thyra::LinearOpBase<Real> > dgdp_op = thyra_model_->create_DgDp_op(j, p_indices_[i]);
        J->setNonconstBlock(0, i, dgdp_op);
    }
    J->endBlockFill();
    return J;
}

template <typename Real>
Teuchos::RCP<Thyra::LinearOpBase<Real> > 
ProductModelEvaluator<Real>::create_DgDp_op(int j, int l, Teuchos::RCP<Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp) const {    
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dgdp_op = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<Real>());
    dgdp_op->beginBlockFill(1, p_indices_.size());
    for(std::size_t i=0; i<p_indices_.size(); ++i) {
        Teuchos::RCP<Thyra::DefaultScaledAdjointLinearOp<Real>> dgdp_i_ALOP;

        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support = DgDp_op_support_[j*p_indices_.size()+i];


        if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
            dgdp_i_ALOP = Teuchos::rcp<Thyra::DefaultScaledAdjointLinearOp<Real>>( 
                new Thyra::DefaultScaledAdjointLinearOp<Real>(1, Thyra::EOpTransp::TRANS, 
                Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real> > (prodvec_dgdp->getNonconstMultiVectorBlock(i))));
        }
        else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
            dgdp_i_ALOP = Teuchos::rcp<Thyra::DefaultScaledAdjointLinearOp<Real>>( 
                new Thyra::DefaultScaledAdjointLinearOp<Real>(1, Thyra::EOpTransp::NOTRANS, 
                Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real> > (prodvec_dgdp->getNonconstMultiVectorBlock(i))));
        }
        dgdp_op->setNonconstBlock(0, i, dgdp_i_ALOP);
    }
    dgdp_op->endBlockFill();
    return dgdp_op;
}

template <typename Real>
void
ProductModelEvaluator<Real>::fromInternalInArgs(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs1, Thyra::ModelEvaluatorBase::InArgsSetup<Real>& inArgs2) const
{
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_dot, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_dot));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_poly, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_poly));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_poly, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_poly));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_mp, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_mp));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_mp, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_mp)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_t, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_t)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_alpha, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_alpha)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_beta, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_beta)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_W_x_dot_dot_coeff, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_W_x_dot_dot_coeff)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_step_size, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_step_size)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_stage_number, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_stage_number)); 
}

template <typename Real>
void
ProductModelEvaluator<Real>::toInternalInArgs(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs1, Thyra::ModelEvaluatorBase::InArgsSetup<Real>& inArgs2) const
{
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_dot, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_dot));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_poly, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_poly));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_poly, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_poly));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_mp, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_dot_mp));
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_x_mp, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_x_mp)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_t, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_t)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_alpha, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_alpha)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_beta, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_beta)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_W_x_dot_dot_coeff, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_W_x_dot_dot_coeff)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_step_size, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_step_size)); 
    inArgs2.setSupports(Thyra::ModelEvaluator<Real>::IN_ARG_stage_number, inArgs1.supports(Thyra::ModelEvaluator<Real>::IN_ARG_stage_number)); 
}

template <typename Real>
void
ProductModelEvaluator<Real>::fromInternalOutArgs(const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs1, Thyra::ModelEvaluatorBase::OutArgsSetup<Real>& outArgs2) const
{
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_mp, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_mp));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_mp, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_mp));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_op, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_op));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_prec, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_prec));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_poly, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_poly));

    for (auto g_index = 0; g_index < outArgs1.Ng(); ++g_index) {
        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx, g_index, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx, g_index));
        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx_dot, g_index, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx_dot, g_index));
        
        bool suppLinOp(true), suppGradMV(true), suppJacMV(true);        
        for (std::size_t i = 0; i < p_indices_.size(); ++i) {
            const auto& dgdpi_support = outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index, p_indices_[i]);
            suppLinOp = dgdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP) ?  suppLinOp : false;
            suppGradMV = dgdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) ?  suppGradMV : false;
            suppJacMV = dgdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) ?  suppJacMV : false;
        }

        Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support;
        if(suppLinOp) dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        if(suppGradMV) dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        if(suppJacMV) dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDp, g_index, 0, dgdp_support);
    }

    { //DfDp
        bool suppLinOp(true), suppGradMV(true), suppJacMV(true);        
        for (std::size_t i = 0; i < p_indices_.size(); ++i) {
            const auto& dfdpi_support = outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, p_indices_[i]);
            suppLinOp = dfdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP) ?  suppLinOp : false;
            suppGradMV = dfdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) ?  suppGradMV : false;
            suppJacMV = dfdpi_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) ?  suppJacMV : false;
        }

        Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support;
        if(suppLinOp) dfdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        if(suppGradMV) dfdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        if(suppJacMV) dfdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DfDp, 0, dfdp_support);
    }

    for (auto g_index = 0; g_index < outArgs1.Ng(); ++g_index) {

        bool all_hess_g_pp = false;
        for (std::size_t i = 0; i < p_indices_.size(); ++i) {
            if (outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_g_pp, g_index, p_indices_[i], p_indices_[i])) {
                if (i == 0) all_hess_g_pp = true;
                if (!all_hess_g_pp)
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                            std::endl <<
                            "ProductModelEvaluator<Real>::fromInternalOutArgs(): hess_g_pp for Solution index = " << 
                            g_index << " is supported for Parameter index i = " << i << " but was not supported for the " <<
                            "previous parameter indices; a non consistent support is not supported. " << std::endl);
            }
            else {
                if (all_hess_g_pp)
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                            std::endl <<
                            "ProductModelEvaluator<Real>::fromInternalOutArgs(): hess_g_pp for Solution index = " << 
                            g_index << " is not supported for Parameter index i = " << i << " but was supported for the " <<
                            "previous parameter indices; a non consistent support is not supported. " << std::endl);
            }
        }
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_g_pp, g_index, 0, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_g_pp, g_index, p_indices_[0], p_indices_[0]));

        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index));
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, p_indices_[0]));
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, p_indices_[0]));
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, 0, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_indices_[0], p_indices_[0]));
    }
    outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xx, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xx));
    outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, p_indices_[0]));
    outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, p_indices_[0]));
    outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, 0, 0, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, p_indices_[0], p_indices_[0]));
}

template <typename Real>
void
ProductModelEvaluator<Real>::toInternalOutArgs(const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs1, Thyra::ModelEvaluatorBase::OutArgsSetup<Real>& outArgs2) const
{
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_mp, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_mp));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_mp, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_mp));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_op, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_op));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_prec, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_W_prec));
    outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_poly, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_f_poly));

    for (auto g_index = 0; g_index < outArgs1.Ng(); ++g_index) {
        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx, g_index, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx, g_index));
        outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx_dot, g_index, outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDx_dot, g_index));
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index));
        for (std::size_t i = 0; i < p_indices_.size(); ++i) {
            outArgs2.setSupports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDp, g_index, p_indices_[i], outArgs1.supports(Thyra::ModelEvaluator<Real>::OUT_ARG_DgDp, g_index, 0));
            outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, p_indices_[i], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, 0));
            outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, p_indices_[i], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, 0));
            for (std::size_t j = 0; j < p_indices_.size(); ++j) {
                outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_indices_[i], p_indices_[j], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, 0, 0));
            }
        }
    }

    outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xx, outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xx));

    for (std::size_t i = 0; i < p_indices_.size(); ++i) {
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, p_indices_[i], DfDp_op_support_[i]);
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, p_indices_[i], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, 0));
        outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, p_indices_[i], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, 0));
        for (std::size_t j = 0; j < p_indices_.size(); ++j) {
            outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, p_indices_[i], p_indices_[j], outArgs1.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, 0, 0));
        }
    }

    for (auto i = 0; i < thyra_model_->Ng(); ++i) {
        for (std::size_t j = 0; j < p_indices_.size(); ++j) {
            outArgs2.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, p_indices_[j], DgDp_op_support_[i*p_indices_.size()+j]);
        }
    }
}

#ifdef HAVE_PIRO_ROL
template <typename Real>
void
ProductModelEvaluator<Real>::block_diagonal_hessian_22(const Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> H,
                    const ROL::Vector<Real> &u,
                    const ROL::Vector<Real> &z,
                    const int g_idx) const
{
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = true;
    for(std::size_t i=0; i<p_indices_.size(); ++i)
      supports_deriv = supports_deriv &&  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_g_pp, g_idx, p_indices_[i], p_indices_[i]);
    
    TEUCHOS_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, 
        "ProductModelEvaluator<Real>::block_diagonal_hessian_22: H_pp is not supported");

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
    unew->set(u);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);

    Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

    H->beginBlockFill(p_indices_.size(), p_indices_.size());

    for(std::size_t i=0; i<p_indices_.size(); ++i) {
      inArgs.set_p(p_indices_[i], thyra_prodvec_p->getVectorBlock(i));
    }
    inArgs.set_x(thyra_x.getVector());

    Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model_->get_g_multiplier_space(g_idx));
    Thyra::put_scalar(1.0, multiplier_g.ptr());
    inArgs.set_g_multiplier(g_idx, multiplier_g);

    for(std::size_t i=0; i<p_indices_.size(); ++i) {
      Teuchos::RCP<Thyra::LinearOpBase<Real>> hess_g_pp = thyra_model_->create_hess_g_pp(g_idx, p_indices_[i], p_indices_[i]);
      outArgs.set_hess_g_pp(g_idx, p_indices_[i], p_indices_[i], hess_g_pp);
      H->setBlock(i, i, hess_g_pp);
    }
    H->endBlockFill();

    thyra_model_->evalModel(inArgs, outArgs);
}
#endif

template <typename Real>
Teuchos::RCP<Piro::ProductModelEvaluator<Real>> getNonconstProductModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<Real>> model) {
    Teuchos::RCP<Piro::ProductModelEvaluator<Real>> model_PME = 
        Teuchos::rcp_dynamic_cast<Piro::ProductModelEvaluator<Real>>(model);
    if (!model_PME.is_null()) {
        return model_PME;
    }
    Teuchos::RCP<Thyra::ModelEvaluator<Real>> model_tmp = model;
    while (true) {
        Teuchos::RCP<Thyra::ModelEvaluatorDelegatorBase<Real>> model_MEDB =
            Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDelegatorBase<Real>>(model_tmp);
        if (!model_MEDB.is_null()) {
            model_tmp = model_MEDB->getNonconstUnderlyingModel();
            //std::cout << model_MEDB->description() << std::endl;
            model_PME = Teuchos::rcp_dynamic_cast<Piro::ProductModelEvaluator<Real>>(model_tmp);
            if (!model_PME.is_null()) {
                return model_PME;
            }
        }
        else
            return Teuchos::null;
    }
}

template <typename Real>
Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> getProductModelEvaluator(const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> model) {
    Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> model_PME = 
        Teuchos::rcp_dynamic_cast<const Piro::ProductModelEvaluator<Real>>(model);
    if (!model_PME.is_null()) {
        return model_PME;
    }
    Teuchos::RCP<const Thyra::ModelEvaluator<Real>> model_tmp = model;
    while (true) {
        Teuchos::RCP<const Thyra::ModelEvaluatorDelegatorBase<Real>> model_MEDB =
            Teuchos::rcp_dynamic_cast<const Thyra::ModelEvaluatorDelegatorBase<Real>>(model_tmp);
        if (!model_MEDB.is_null()) {
            model_tmp = model_MEDB->getUnderlyingModel();
            //std::cout << model_MEDB->description() << std::endl;
            model_PME = Teuchos::rcp_dynamic_cast<const Piro::ProductModelEvaluator<Real>>(model_tmp);
            if (!model_PME.is_null()) {
                return model_PME;
            }
        }
        else
            return Teuchos::null;
    }
}

template <typename Real>
Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> getProductModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Real>> model) {
    return getProductModelEvaluator(Teuchos::rcp_dynamic_cast<const Thyra::ModelEvaluator<Real>>(model));
}

template <typename Real>
Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> getProductModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<Real>> model) {
    return getProductModelEvaluator(Teuchos::rcp_dynamic_cast<const Thyra::ModelEvaluator<Real>>(model));
}

} // namespace Piro

#endif