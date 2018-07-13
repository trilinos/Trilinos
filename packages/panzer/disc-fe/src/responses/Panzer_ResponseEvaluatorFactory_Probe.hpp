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

#ifndef __Panzer_ResponseEvaluatorFactory_Probe_hpp__
#define __Panzer_ResponseEvaluatorFactory_Probe_hpp__

#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ResponseMESupportBuilderBase.hpp"

#include <mpi.h>

namespace panzer {

/** This class defines a response based on a point value.
  */
template <typename EvalT,typename LO,typename GO>
class ResponseEvaluatorFactory_Probe : public ResponseEvaluatorFactory<EvalT> {
public:

  ResponseEvaluatorFactory_Probe(
    MPI_Comm comm,
    const Teuchos::Array<double>& point,
    int fieldComponent = 0,
    int cubatureDegree=1,
    const std::string & fieldName="",
    const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory=Teuchos::null,
    const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & globalIndexer=Teuchos::null,
    bool applyDirichletToDerivative=false)
    : comm_(comm), point_(point), fieldComponent_(fieldComponent), cubatureDegree_(cubatureDegree)
    , fieldName_(fieldName), linearObjFactory_(linearObjFactory), globalIndexer_(globalIndexer)
    , applyDirichletToDerivative_(applyDirichletToDerivative)
    {
      TEUCHOS_ASSERT((linearObjFactory==Teuchos::null && globalIndexer==Teuchos::null) ||
                     (linearObjFactory!=Teuchos::null && globalIndexer!=Teuchos::null));
    }

  virtual ~ResponseEvaluatorFactory_Probe() {}

  /** Build the response object used by this factory. This object
   * assumes the role of the scatter target and will be accessible
   * by all the evaluators in the field managers.
   *
   * \param[in] responseName Name of response to be built. This
   *                         name will be used for looking up
   *                         the response in the <code>GlobalEvaluationDataContainer</code>
   *                         object.
   */
  virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName) const;

  virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName,
                                                         const std::vector<WorksetDescriptor> & wkstDesc) const
    { return buildResponseObject(responseName); }

  /** Build and register evaluators for a response on a particular physics
   * block.
   *
   * \param[in] responseName The name of the response to be constructed
   *                         by these evaluators.
   * \param[in,out] fm Field manager to be fuild with the evaluators.
   * \param[in] physicsBlock What physics block is being used for constructing
   *                         the evaluators
   * \param[in] user_data The user data parameter list, this stores things
   *                      that the user may find useful.
   */
  virtual void buildAndRegisterEvaluators(const std::string & responseName,
                                          PHX::FieldManager<panzer::Traits> & fm,
                                          const panzer::PhysicsBlock & physicsBlock,
                                          const Teuchos::ParameterList & user_data) const;

  /** Is this evaluation type supported by the factory. This is used to determine cases
   * where a response may support a particular evaluation type, however at runtime the user
   * decides not to enable the (say) Jacobian evaluation of this response.
   *
   * Note that use of this mechanism is complementary to having the builder return
   * <code>Teuchos::null</code> for a particular evaluation type.
   */
  virtual bool typeSupported() const;

protected:
  //! Accessor method for Cubature degree (can be used by sub classes)
  int getCubatureDegree() const { return cubatureDegree_; }

private:
  MPI_Comm comm_;
  Teuchos::Array<double> point_;
  int fieldComponent_;
  int cubatureDegree_;
  std::string fieldName_;
  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linearObjFactory_;
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  bool applyDirichletToDerivative_; // do we need this???
};

template <typename LO,typename GO>
struct ProbeResponse_Builder : public ResponseMESupportBuilderBase {
  MPI_Comm comm;
  Teuchos::Array<double> point;
  int fieldComponent;
  int cubatureDegree;
  std::string fieldName;
  bool applyDirichletToDerivative; // if this is set to true, then the dirichlet values will be zerod out in
                                   // the DgDx vector

  ProbeResponse_Builder() : applyDirichletToDerivative(false) {}

  virtual ~ProbeResponse_Builder() {}

  void setDerivativeInformation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & in_linearObjFactory,
                                const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & in_globalIndexer)
  {
    linearObjFactory = in_linearObjFactory;
    globalIndexer = in_globalIndexer;

    TEUCHOS_ASSERT((linearObjFactory==Teuchos::null && globalIndexer==Teuchos::null) ||
                   (linearObjFactory!=Teuchos::null && globalIndexer!=Teuchos::null));
  }

  virtual void setDerivativeInformation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & in_linearObjFactory)
  {
    using Teuchos::rcp_dynamic_cast;

    setDerivativeInformation(in_linearObjFactory,
                             rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<LO,GO> >(in_linearObjFactory->getDomainGlobalIndexer(),true));
  }

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_Probe<T,LO,GO>(comm,point,fieldComponent,cubatureDegree,fieldName,
                                                                    linearObjFactory,globalIndexer,
                                                                    applyDirichletToDerivative)); }

  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildValueFactory() const
  { return build<panzer::Traits::Residual>(); }

  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildDerivativeFactory() const
  { return build<panzer::Traits::Jacobian>(); }

  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildTangentFactory() const
  { return build<panzer::Traits::Tangent>(); }

private:
  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer;
};


}

#include "Panzer_ResponseEvaluatorFactory_Probe_impl.hpp"

#endif
