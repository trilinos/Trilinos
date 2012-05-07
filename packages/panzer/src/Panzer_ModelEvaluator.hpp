// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_MODEL_EVALUATOR_DECL_HPP
#define PANZER_MODEL_EVALUATOR_DECL_HPP

#include "Panzer_config.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_ParameterLibrary.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_AbstractFactory.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include <Kokkos_DefaultNode.hpp>

namespace panzer {

template<typename, typename>  class FieldManagerBuilder;
template<typename> class ResponseLibrary;
template<typename> class LinearObjFactory;
class GlobalData;

template<typename Scalar, typename LO, typename GO, typename NODE>
class ModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

//   typedef typename panzer::Traits<T>::lid_type LO;
//   typedef typename panzer::Traits<T>::gid_type GO;
//   typedef typename panzer::Traits<T>::node_type NODE;

public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  ModelEvaluator(const Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> >& fmb,
                 const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
		 const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
		 const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
                 const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
		 const Teuchos::RCP<panzer::GlobalData>& global_data,
		 bool build_transient_support);

  /** \brief . */
  ModelEvaluator();

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;

  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  //@}

private: // data members

  double t_init_;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  // parameters
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > p_space_;
  std::vector<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > p_init_;

  // responses
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > g_space_;

  Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> > fmb_;
  mutable Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > responseLibrary_; // These objects are basically the same
  mutable panzer::AssemblyEngine_TemplateManager<panzer::Traits,LO,GO> ae_tm_;     // they control and provide access to evaluate
  std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names_;

  mutable Teuchos::Array<panzer::ParamVec> parameter_vector_;
  Teuchos::RCP<panzer::GlobalData> global_data_;
  bool build_transient_support_;

  // basic specific linear object objects
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof_;
  mutable Teuchos::RCP<panzer::LinearObjContainer> ghostedContainer_;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory_;
};


}

#include "Panzer_ModelEvaluator_impl.hpp"

#endif 
