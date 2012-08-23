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

#ifndef __Panzer_ResponseAggregator_IPCoordinates_decl_hpp__
#define __Panzer_ResponseAggregator_IPCoordinates_decl_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_Traits.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_IPCoordinates.hpp"

namespace panzer {

template <typename EvalT,typename TraitsT> class ResponseAggregator_IPCoordinates_Data;
template <typename EvalT,typename TraitsT> class ResponseAggregator_IPCoordinates;

template <typename TraitsT> 
class ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> 
   : public ResponseDataDefault<TraitsT> {
public:
   ResponseAggregator_IPCoordinates_Data() {}

  /** \defgroup Required from ResponseData
   * @{
   */
  
  virtual void allocateAndInitializeData(const std::vector<std::string> & fields)
  { 
    block_id = Teuchos::rcp(new std::string);
    coords = Teuchos::rcp(new std::vector<panzer::Traits::Residual::ScalarT>);
    reinitializeData();
  }
  
  virtual void reinitializeData() {}
  
  virtual void fillResponse(const std::string & field,Response<TraitsT> & response) const
  { 
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Element Block ID",*block_id);
    pl->set("IP Coordinates",coords);
    response.setParameterList(pl);
  }
  
   //! @}
  Teuchos::RCP<const std::string> getBlockID() const
  { return block_id; }

  Teuchos::RCP<std::string> getNonconstBlockID()
  { return block_id; }

  Teuchos::RCP<const std::vector<panzer::Traits::Residual::ScalarT> > getCoords() const
  { return coords; }

  Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > getNonconstCoords()
  { return coords; }

private:
  Teuchos::RCP<std::string> block_id;
  Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > coords;
};

template <typename TraitsT>
class ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>
   : public ResponseAggregator<panzer::Traits::Residual,TraitsT> {
public:
   // useful for cloning and the factory mechanism
   ResponseAggregator_IPCoordinates();

   ResponseAggregator_IPCoordinates(const Teuchos::ParameterList & p);

   /** \defgroup Methods required by ResponseAggregatorBase
     * @{
     */

   //! Clone this aggregator
   Teuchos::RCP<ResponseAggregatorBase<TraitsT> > clone(const Teuchos::ParameterList & p) const;

   //! Build response data for a specified set of fields (ResponseAggregator decides data layout)
   Teuchos::RCP<ResponseData<TraitsT> > buildResponseData(const std::vector<std::string> & fields) const;

   //! Register and build evaluator required by this aggregator and this set of data.
   virtual void registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<ResponseData<TraitsT> > & data, 
                                             const PhysicsBlock & pb,
                                             const Teuchos::ParameterList & p) const;

   //! perform global reduction on this set of response data
   void globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const;

   //! Aggregate a set of responses locally
   virtual void aggregateResponses(Response<TraitsT> & dest,const std::list<Teuchos::RCP<const Response<TraitsT> > > & sources) const;

   //! @}

private:

  //! set to true first time coordinates are evaluated.
  mutable bool first;
};

// Specialized for panzer::Traits
class ResponseAggregator_IPCoordinates_Builder {
public:
   void setGlobalIndexer(const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi)
   { globalIndexer_ = ugi; }

   void setLinearObjFactory(const Teuchos::RCP<LinearObjFactory<panzer::Traits> > & lof)
   { linObjFactory_ = lof; }

   Teuchos::RCP<UniqueGlobalIndexerBase> getGlobalIndexer() const
   { return globalIndexer_; }

   Teuchos::RCP<LinearObjFactory<panzer::Traits> > getLinearObjFactory() const
   { return linObjFactory_; }

   template <typename EvalT>
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > build() const
   { return Teuchos::null; }

private:
   Teuchos::RCP<UniqueGlobalIndexerBase> globalIndexer_;
   Teuchos::RCP<LinearObjFactory<panzer::Traits> > linObjFactory_;
};

// declaration so methods are not inlined
template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseAggregator_IPCoordinates_Builder::build<panzer::Traits::Residual>() const;

}

#endif
