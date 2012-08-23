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

#ifndef __Panzer_ResponseAggregatorBase_hpp__
#define __Panzer_ResponseAggregatorBase_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_Response.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Comm.hpp"

#include "Phalanx_FieldManager.hpp"

namespace panzer {

/** Abstract data object for response storage. Also manages
  * setting up and transfering to the responses. This works
  * in concert with the response aggregator.
  */
template <typename TraitsT>
class ResponseData {
public:
   virtual ~ResponseData() {}

   //! Allocate and initialize required storage based on the fields passed in
   virtual void allocateAndInitializeData(const std::vector<std::string> & fields) = 0;

   /** Reinitialize the data based on the fields original requested by
     * <code>allocateAndInitializeData</code>.
     */
   virtual void reinitializeData() = 0;

   /** Fill a response object with data from a particular field.
     */
   virtual void fillResponse(const std::string & field,Response<TraitsT> & response) const = 0;

   //! Is this field contain in the response data
   virtual bool contains(const std::string & field) const = 0;

   //! Get a vector of fields in this data object
   virtual const std::vector<std::string> & getFields() const = 0;

   //! Set a vector of fields in this data object
   virtual void setFields(const std::vector<std::string> & fields) = 0;
};

/** A default implementation that handles basic field string management
  */
template <typename TraitsT>
class ResponseDataDefault : public ResponseData<TraitsT> {
public:
   virtual ~ResponseDataDefault() {}

   //! Is this field contain in the response data
   virtual bool contains(const std::string & field) const
   { return this->fieldIndex(field)<fields_.size(); }
   
   //! Get a vector of fields in this data object
   virtual const std::vector<std::string> & getFields() const
   { return fields_; }

   //! Set a vector of fields in this data object
   virtual void setFields(const std::vector<std::string> & fields) 
   { fields_ = fields; }

   //! Get field index
   std::size_t fieldIndex(const std::string & field) const
   {
      std::vector<std::string>::const_iterator itr 
            = std::find(fields_.begin(),fields_.end(),field);
      return itr-fields_.begin();
   }

private:
   std::vector<std::string> fields_; // Storing fields in this data object
};

template <typename TraitsT>
class ResponseAggregatorBase {
public:
   virtual ~ResponseAggregatorBase() {}

   //! Clone this aggregator
   virtual Teuchos::RCP<ResponseAggregatorBase<TraitsT> > clone(const Teuchos::ParameterList & p) const = 0;

   //! Build response data for a specified set of fields (ResponseAggregator decides data layout)
   virtual Teuchos::RCP<ResponseData<TraitsT> > buildResponseData(const std::vector<std::string> & fields) const = 0;

   //! Register and build evaluator required by this aggregator and this set of data.
   virtual void registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<ResponseData<TraitsT> > & data,
                                             const PhysicsBlock & pb,
                                             const Teuchos::ParameterList & p) const = 0;

   //! perform global reduction on this set of response data
   virtual void globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const = 0;

   /** \brief Aggregate a set of responses locally over a set of element blocks
    *
    *  This aggregation takes place after the globalReduction().  It
    *  will take the responses from each element block and aggregate
    *  them into a single destination response that is created by
    *  reserving a labeled response on the library.  Therefore this
    *  method will only be called if the user calls
    *  getBlockAggregatedVolumeResponseByLabel() on the response library.
    */
   virtual void aggregateResponses(Response<TraitsT> & dest,const std::list<Teuchos::RCP<const Response<TraitsT> > > & sources) const = 0;

   virtual void setGlobalIndexer(const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi) = 0;

   virtual void setLinearObjFactory(const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof) = 0;
};

template <typename EvalT,typename TraitsT>
class ResponseAggregator : public ResponseAggregatorBase<TraitsT> {
public:
   virtual ~ResponseAggregator() {}

   virtual void setGlobalIndexer(const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi)
   { globalIndexer_ = ugi; }

   virtual void setLinearObjFactory(const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof)
   { linObjFactory_ = lof; }

   Teuchos::RCP<UniqueGlobalIndexerBase> getGlobalIndexer() const
   { return globalIndexer_; }

   Teuchos::RCP<LinearObjFactory<TraitsT> > getLinearObjFactory() const
   { return linObjFactory_; }

private:
   Teuchos::RCP<UniqueGlobalIndexerBase> globalIndexer_;
   Teuchos::RCP<LinearObjFactory<TraitsT> > linObjFactory_;
};

}

#endif
