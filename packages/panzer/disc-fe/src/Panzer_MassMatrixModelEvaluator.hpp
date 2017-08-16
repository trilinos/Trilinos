// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
