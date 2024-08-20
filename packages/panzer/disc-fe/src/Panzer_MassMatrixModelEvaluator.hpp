// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MASS_MATRIX_MODEL_EVALUATOR_DECL_HPP
#define PANZER_MASS_MATRIX_MODEL_EVALUATOR_DECL_HPP

#include "PanzerDiscFE_config.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_ResponseMESupportBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_AbstractFactory.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include <Panzer_NodeType.hpp>

namespace panzer {

// Class for explicit model evaluators that construct a mass matrix
// A mass matrix inverse may or may not be required within the explicit evaluator at different parts of an IMEX algorithm
// This class allows this option to be switched on and off

template<typename Scalar>
class MassMatrixModelEvaluator
  : public Thyra::ModelEvaluatorBase
{
public:

void setApplyMassInverse(const bool applyMassInverse) const
{
  applyMassInverse_ = applyMassInverse;
}

virtual void applyMassMatrix(const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > input, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > output) const = 0;
virtual void applyInverseMassMatrix(const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > input, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > output) const = 0;

protected:

//! Apply mass matrix inverse within the evaluator
mutable bool applyMassInverse_;

private:

};


}


#endif 
