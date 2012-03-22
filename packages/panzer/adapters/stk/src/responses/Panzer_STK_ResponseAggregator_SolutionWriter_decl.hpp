#ifndef __Panzer_ResponseAggregator_SolutionWriter_decl_hpp__
#define __Panzer_ResponseAggregator_SolutionWriter_decl_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_Traits.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseData_Action.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace panzer_stk {

template <typename EvalT,typename TraitsT>
class ResponseAggregator_SolutionWriter
   : public panzer::ResponseAggregator<EvalT,TraitsT> {
public:
   // useful for cloning and the factory mechanism
   ResponseAggregator_SolutionWriter(const Teuchos::RCP<STK_Interface> & mesh);

   /** \defgroup Methods required by ResponseAggregatorBase
     * @{
     */

   //! Clone this aggregator
   Teuchos::RCP<panzer::ResponseAggregatorBase<TraitsT> > clone(const Teuchos::ParameterList & p) const
   { return Teuchos::rcp(new ResponseAggregator_SolutionWriter<EvalT,TraitsT>(mesh_)); }

   Teuchos::RCP<panzer::ResponseData<TraitsT> > buildResponseData(const std::vector<std::string> & fields) const
   { return Teuchos::rcp(new panzer::ResponseData_Action<TraitsT>()); }

   //! Register and build evaluators required by this aggregator and this set of data.
   virtual void registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<panzer::ResponseData<TraitsT> > & data, 
                                             const panzer::PhysicsBlock & pb,
                                             const Teuchos::ParameterList & p) const;

   //! perform global reduction on this set of response data
   void globalReduction(const Teuchos::Comm<int> & comm,panzer::ResponseData<TraitsT>  & rd) const {}

   //! Aggregate a set of responses locally
   virtual void aggregateResponses(panzer::Response<TraitsT> & dest,const std::list<Teuchos::RCP<const panzer::Response<TraitsT> > > & sources) const {}

   //! @}
  
   /** Take a vector of (std::string (field name), RCP<PureBasis>) pairs and bucket them
     * by basis name. What is returned is a map pairing the basis to a vector of field names.
     */
   static void bucketByBasisType(const std::vector<panzer::StrPureBasisPair> & providedDofs,
                                 std::map<std::string,std::vector<std::string> > & basisBucket);


private:
   void computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases,
                                 int baseDimension,
                                 Intrepid::FieldContainer<double> & centroid) const;

   Teuchos::RCP<STK_Interface> mesh_;
};

// Specialized for panzer::Traits
class ResponseAggregator_SolutionWriter_Builder {
public:
   ResponseAggregator_SolutionWriter_Builder(const Teuchos::RCP<STK_Interface> & mesh)
      : mesh_(mesh) {}

   void setGlobalIndexer(const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > & ugi)
   { globalIndexer_ = ugi; }

   void setLinearObjFactory(const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof)
   { linObjFactory_ = lof; }

   Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > getGlobalIndexer() const
   { return globalIndexer_; }

   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > getLinearObjFactory() const
   { return linObjFactory_; }

   template <typename EvalT>
   Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > build() const
   { return Teuchos::null; }

private:
   Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > globalIndexer_;
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;

   Teuchos::RCP<STK_Interface> mesh_;
};

// declaration so methods are not inlined
template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > ResponseAggregator_SolutionWriter_Builder::build<panzer::Traits::Residual>() const;

#ifdef HAVE_STOKHOS
template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > ResponseAggregator_SolutionWriter_Builder::build<panzer::Traits::SGResidual>() const;
#endif

}

#endif
