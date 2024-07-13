// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "LOCA_Tpetra_ConstraintModelEvaluator.hpp" // class declaration
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "NOX_Thyra_MultiVector.H"
#include "NOX_Thyra_Vector.H"
#include "NOX_TpetraTypedefs.hpp"

namespace LOCA {
  namespace MultiContinuation {

    ConstraintModelEvaluator::ConstraintModelEvaluator(const Teuchos::RCP<::Thyra::ModelEvaluator<double>>& model,
                                                       const LOCA::ParameterVector& pVec,
                                                       const std::vector<std::string>& constraintResponseNames,
                                                       const NOX::Abstract::Vector& cloneVec) :
      model_(model),
      pVec_(pVec),
      gNames_(constraintResponseNames),
      constraints_(constraintResponseNames.size(),1),
      isValidConstraints_(false),
      isValidDx_(false),
      printDebug_(false)
    {
      x_ = cloneVec.clone(NOX::ShapeCopy);
      dgdx_ = cloneVec.createMultiVector(constraintResponseNames.size(),NOX::ShapeCopy);

      TEUCHOS_ASSERT(constraintResponseNames.size() <= size_t(pVec_.length()));

      // Get the parameter indices
      meParameterIndices_.clear();
      const auto pNames = pVec_.getNamesVector();
      for (size_t i=0; i < pNames.size(); ++i) {
        const auto& p_name = pNames[i];
        bool found = false;
        for (int j=0; j < model_->Np(); ++j) {
          const auto& names_vec = model_->get_p_names(j);
          const auto search = std::find(names_vec->begin(),names_vec->end(),p_name);
          if (search != names_vec->end()) {
            meParameterIndices_.push_back(j);
            found = true;
            break;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(!found,std::runtime_error,"ERROR LOCA::TpetraConstraintModelEvaluator::CTOR - could not find parameter named \"" << p_name << "\" in the model evaluator!");
      }
      TEUCHOS_ASSERT(pNames.size() == meParameterIndices_.size());

      // Get the response indices
      meResponseIndices_.clear();
      for (size_t i=0; i < gNames_.size(); ++i) {
        const auto& g_name = gNames_[i];
        bool found = false;
        for (int j=0; j < model_->Ng(); ++j) {
          const auto& names_vec = model_->get_g_names(j);
          const auto search = std::find(names_vec.begin(),names_vec.end(),g_name);
          if (search != names_vec.end()) {
            meResponseIndices_.push_back(j);
            found = true;
            break;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(!found,std::runtime_error,"ERROR LOCA::TpetraConstraintModelEvaluator::CTOR - could not find response named \"" << g_name << "\" in the model evaluator!");
      }
      TEUCHOS_ASSERT(gNames_.size() == meResponseIndices_.size());

      // Allocate all objects to pass into model evaluator.  The
      // parameters and responses are treated as separate objects.
      me_p_.resize(pVec_.length());
      for (size_t l=0; l < me_p_.size(); ++l)
        me_p_[l] = ::Thyra::createMember(*model->get_p_space(meParameterIndices_[l]),"p_l");

      const auto numResponses = gNames_.size();
      me_g_.resize(numResponses);
      me_dgdx_.resize(numResponses);
      me_dgdp_.resize(numResponses);
      for (auto& i : me_dgdp_)
        i.resize(me_p_.size());
      for (size_t j=0; j < numResponses; ++j) {
        me_g_[j] = ::Thyra::createMember(*model->get_g_space(meResponseIndices_[j]),"g_j");
        me_dgdx_[j] = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<double>>(model->create_DgDx_op(meResponseIndices_[j]),false);
        if (me_dgdx_[j].is_null()) {
          // It might be wrapped in an Adjoint linear op for the transpose.
          auto ptr = Teuchos::rcp_dynamic_cast<::Thyra::DefaultScaledAdjointLinearOp<double>>(model->create_DgDx_op(meResponseIndices_[j]),true);
          me_dgdx_[j] = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<double>>(ptr->getNonconstOp(),true);
        }
        for (size_t l=0; l < me_p_.size(); ++l)
          me_dgdp_[j][l] = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<double>>(model->create_DgDp_op(meResponseIndices_[j],
                                                                                                             meParameterIndices_[l]),
                                                                                       true);
      }

    }

    ConstraintModelEvaluator::ConstraintModelEvaluator(const LOCA::MultiContinuation::ConstraintModelEvaluator& src,
                                                       NOX::CopyType type) :
      ConstraintModelEvaluator(src.model_,src.pVec_,src.gNames_,*src.x_)
    {
      if (type == NOX::DeepCopy) {
        *x_ = *src.x_;
        pVec_ = src.pVec_;
        for (int i=0; i < pVec_.length(); ++i)
          ::Thyra::V_S(me_p_[i].ptr(),pVec_[i]);
        constraints_ = src.constraints_;
        *dgdx_ = *src.dgdx_;
        isValidConstraints_ = src.isValidConstraints_;
        isValidDx_ = src.isValidDx_;
        printDebug_ = src.printDebug_;
      }
    }

    ConstraintModelEvaluator::~ConstraintModelEvaluator(){}

    void ConstraintModelEvaluator::copy(const ConstraintInterface& source)
    {
      const auto& src = dynamic_cast<const LOCA::MultiContinuation::ConstraintModelEvaluator&>(source);
      if (this != &src) {
        *x_ = *src.x_;
        pVec_ = src.pVec_;
        for (int i=0; i < pVec_.length(); ++i)
          ::Thyra::V_S(me_p_[i].ptr(),pVec_[i]);
        constraints_ = src.constraints_;
        *dgdx_ = *src.dgdx_;
        isValidConstraints_ = src.isValidConstraints_;
        isValidDx_ = src.isValidDx_;
        meParameterIndices_ = src.meParameterIndices_;
        meResponseIndices_ = src.meResponseIndices_;
        printDebug_ = src.printDebug_;
      }
    }

    Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
    ConstraintModelEvaluator::clone(NOX::CopyType type) const
    {
      return Teuchos::rcp(new LOCA::MultiContinuation::ConstraintModelEvaluator(*this,type));
    }

    int ConstraintModelEvaluator::numConstraints() const
    {
      return int(constraints_.numRows());
    }

    void ConstraintModelEvaluator::setX(const NOX::Abstract::Vector& x)
    {
      *x_ = x; // copy values into local storage
      isValidConstraints_ = false;
      isValidDx_ = false;
    }

    void ConstraintModelEvaluator::setParam(int paramID, double val)
    {
      (pVec_)[paramID] = val;
      ::Thyra::V_S(me_p_[paramID].ptr(),val);
      isValidConstraints_ = false;
      isValidDx_ = false;
    }

    void ConstraintModelEvaluator::setParams(const std::vector<int>& paramIDs,
                                             const NOX::Abstract::MultiVector::DenseMatrix& vals)
    {
      for (unsigned int i=0; i<paramIDs.size(); i++) {
        (pVec_)[paramIDs[i]] = vals(i,0);
        ::Thyra::V_S(me_p_[paramIDs[i]].ptr(),vals(i,0));
      }
      isValidConstraints_ = false;
      isValidDx_ = false;
    }

    NOX::Abstract::Group::ReturnType ConstraintModelEvaluator::computeConstraints()
    {
      auto inArgs = model_->createInArgs();
      auto x_thyra = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(x_)->getThyraRCPVector();
      inArgs.set_x(x_thyra);
      for (int i=0; i < pVec_.length(); ++i)
        inArgs.set_p(meParameterIndices_[i],me_p_[i]);

      auto outArgs = model_->createOutArgs();
      for (size_t i=0; i < me_g_.size(); ++i)
        outArgs.set_g(meResponseIndices_[i],::Thyra::ModelEvaluatorBase::Evaluation<::Thyra::VectorBase<double>>(me_g_[i]));

      model_->evalModel(inArgs,outArgs);

      // Now copy responses into the aggregated constraint object for
      // loca.
      using extractor = ::Thyra::TpetraOperatorVectorExtraction<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
      for (size_t i=0; i < me_g_.size(); ++i) {
        auto tmp = extractor::getTpetraMultiVector(me_g_[i]);
        auto val = tmp->getLocalViewHost(Tpetra::Access::ReadOnly);
        constraints_(i,0) = val(0,0);
        if (printDebug_)
          std::cout << "LOCA::ConstraintME: constraints_(" << i << ")=" << val(0,0) << std::endl;
      }

      isValidConstraints_ = true;
      return NOX::Abstract::Group::Ok;
    }

    NOX::Abstract::Group::ReturnType ConstraintModelEvaluator::computeDX()
    {
      auto inArgs = model_->createInArgs();
      auto x_thyra = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(x_)->getThyraRCPVector();
      inArgs.set_x(x_thyra);
      for (int i=0; i < pVec_.length(); ++i)
        inArgs.set_p(meParameterIndices_[i],me_p_[i]);

      auto outArgs = model_->createOutArgs();
      for (size_t i=0; i < me_dgdx_.size(); ++i)
        outArgs.set_DgDx(meResponseIndices_[i],::Thyra::ModelEvaluatorBase::Derivative<NOX::Scalar>(me_dgdx_[i],::Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM));

      model_->evalModel(inArgs,outArgs);

      // Aggregate derivative columns into a single
      // multivector. Copies entries out of temporary me_dgdx_ into
      // the aggregated object loca expects.
      for (size_t i=0; i < me_dgdx_.size(); ++i) {
        auto tmp = NOX::Thyra::MultiVector(me_dgdx_[i]);
        (*dgdx_)[i] = tmp[0];

        if (printDebug_) {
          std::cout << "LOCA::ConstraintME: computeDX:" << std::endl;
          me_dgdx_[i]->describe(std::cout,Teuchos::VERB_EXTREME);
        }
      }

      isValidDx_ = true;
      return NOX::Abstract::Group::Ok;
    }

    NOX::Abstract::Group::ReturnType
    ConstraintModelEvaluator::computeDP(const std::vector<int>& paramIDs,
                                        NOX::Abstract::MultiVector::DenseMatrix& dgdp,
                                        bool isValidG)
    {
      if (!isValidG) {
        if (!this->isConstraints())
          this->computeConstraints();

        for (size_t i=0; i < me_g_.size(); ++i)
          dgdp(i,0) = constraints_(i,0);

        if (printDebug_) {
          for (size_t i=0; i < me_g_.size(); ++i)
            std::cout << "LOCA::ConstraintME: computeDP: dgdp_(" << i << ",0)=" << dgdp(i,0) << std::endl;
        }
      }

      auto inArgs = model_->createInArgs();
      auto x_thyra = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(x_)->getThyraRCPVector();
      inArgs.set_x(x_thyra);

      // Set all parameters from the parameter vector
      for (int i=0; i < pVec_.length(); ++i)
        inArgs.set_p(meParameterIndices_[i],me_p_[i]);

      // Only request derivatives for the incoming paramIDs. NOTE:
      // this could be, and most likely is, a subset of the total
      // number of parameters in the pVec_.
      TEUCHOS_ASSERT(size_t(dgdp.numRows()) == me_dgdp_.size());
      TEUCHOS_ASSERT(size_t(dgdp.numCols()) == (1+paramIDs.size()));
      auto outArgs = model_->createOutArgs();
      for (size_t j=0; j < me_dgdp_.size(); ++j) {
        for (auto& l : paramIDs) {
          outArgs.set_DgDp(meResponseIndices_[j],meParameterIndices_[l],me_dgdp_[j][l]);

          if (printDebug_) {
            std::cout << "LOCA::ConstraintME: computeDP: dgdp_(" << j << "," << l << "):" << std::endl;
            me_dgdp_[j][l]->describe(std::cout,Teuchos::VERB_EXTREME);
          }
        }
      }

      model_->evalModel(inArgs,outArgs);

      // Extract values of individual elements and copy into
      // aggregated object for loca.
      using extractor = ::Thyra::TpetraOperatorVectorExtraction<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
      for (size_t j=0; j < me_dgdp_.size(); ++j) {
        for (size_t l=0; l < paramIDs.size(); ++l) {
          auto tmp = extractor::getTpetraMultiVector(me_dgdp_[j][paramIDs[l]]);
          auto val = tmp->getLocalViewHost(Tpetra::Access::ReadOnly);
          // first col contains g, so we shift the columns by one
          dgdp(j,l+1) = val(0,0);
        }
      }

      return NOX::Abstract::Group::Ok;
    }

    bool ConstraintModelEvaluator::isConstraints() const
    {return isValidConstraints_;}

    bool ConstraintModelEvaluator::isDX() const
    {return isValidDx_;}

    const NOX::Abstract::MultiVector::DenseMatrix&
    ConstraintModelEvaluator::getConstraints() const
    {return constraints_;}

    bool ConstraintModelEvaluator::isDXZero() const
    {return false;}

    NOX::Abstract::MultiVector *
    ConstraintModelEvaluator::getDX () const
    {return dgdx_.get();}

    const LOCA::ParameterVector
    ConstraintModelEvaluator::getParams() const
    {return pVec_;}

  } // namespace MultiContinuation
} //  namespace LOCA
