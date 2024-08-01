// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef LOCA_TPETRA_CONSTRAINT_MODEL_EVALUATOR_HPP
#define LOCA_TPETRA_CONSTRAINT_MODEL_EVALUATOR_HPP

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_RCP.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include <string>
#include <vector>

namespace Thyra {
  template<typename T> class ModelEvaluator;
}

namespace LOCA {

  namespace MultiContinuation {

    /** \brief Generic object that provides constraints through model evaluator responses. */
    class ConstraintModelEvaluator : virtual public ConstraintInterfaceMVDX {
    public:

      /** Constructor

          \param model Model evaluator that provides constraints as responses.
          \param pVec The independent parameters for constraints.
          \param constrantResponseNames The names of the responses used as constraint equations.
          \param cloneVec NOX vector used to clone data structures with same space/map.
       */
      ConstraintModelEvaluator(const Teuchos::RCP<::Thyra::ModelEvaluator<double>>& model,
                               const LOCA::ParameterVector& pVec,
                               const std::vector<std::string>& constraintResponseNames,
                               const NOX::Abstract::Vector& cloneVec);

      ConstraintModelEvaluator(const LOCA::MultiContinuation::ConstraintModelEvaluator& cme,
                               NOX::CopyType type = NOX::DeepCopy);

      ~ConstraintModelEvaluator();

      void copy(const ConstraintInterface& source);

      Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      int numConstraints() const;

      void setX(const NOX::Abstract::Vector& x);

      /// Set parameter value given a parameter indices corresponding to the LOCA::ParameterVector.
      void setParam(int paramID, double val);

      /// Set parameter value given a parameter indices corresponding to the LOCA::ParameterVector.
      void setParams(const std::vector<int>& paramIDs,
                     const NOX::Abstract::MultiVector::DenseMatrix& vals);

      NOX::Abstract::Group::ReturnType computeConstraints();

      NOX::Abstract::Group::ReturnType computeDX();

      NOX::Abstract::Group::ReturnType
      computeDP(const std::vector<int>& paramIDs,
                NOX::Abstract::MultiVector::DenseMatrix& dgdp,
                bool isValidG);

      bool isConstraints() const;

      bool isDX() const;

      const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const;

      bool isDXZero() const;

      NOX::Abstract::MultiVector * getDX () const;

      // Convenience method for unit testing
      const LOCA::ParameterVector getParams() const;

    private:

      // The underlying model evaluator used to compute responses
      Teuchos::RCP<::Thyra::ModelEvaluator<double>> model_;

      // Parameter vector. Not sure why this is here unless it comes from the application code?
      LOCA::ParameterVector pVec_;

      // Name of responses in model evaluator.
      std::vector<std::string> gNames_;

      // The indices used to set the parameters in the ModelEvalautor inArgs.
      std::vector<int> meParameterIndices_;

      // The indices used to set the responses in the ModelEvalutor outArgs.
      std::vector<int> meResponseIndices_;

      // LOCA data
      NOX::Abstract::MultiVector::DenseMatrix constraints_;
      Teuchos::RCP<NOX::Abstract::Vector> x_;
      Teuchos::RCP<NOX::Abstract::MultiVector> dgdx_;
      bool isValidConstraints_;
      bool isValidDx_;

      // Model Evaluator parameter/response objects
      std::vector<Teuchos::RCP<::Thyra::VectorBase<double>>> me_p_;
      std::vector<Teuchos::RCP<::Thyra::VectorBase<double>>> me_g_;
      std::vector<Teuchos::RCP<::Thyra::MultiVectorBase<double>>> me_dgdx_;
      std::vector<std::vector<Teuchos::RCP<::Thyra::MultiVectorBase<double>>>> me_dgdp_;

      bool printDebug_;
    };

  } // namespace MultiContinuation
} //  namespace LOCA

#endif
