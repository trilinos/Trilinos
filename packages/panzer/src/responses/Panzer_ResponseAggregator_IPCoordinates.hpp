#ifndef __Panzer_ResponseAggregator_IPCoordinates_hpp__
#define __Panzer_ResponseAggregator_IPCoordinates_hpp__

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
    coords = Teuchos::rcp(new std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate>);
    reinitializeData();
  }
  
  virtual void reinitializeData() {}
  
  virtual void fillResponse(const std::string & field,Response<TraitsT> & response) const
  { 
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("IP Coordinates",coords);
    response.setParameterList(pl);
  }
  
   //! @}
 
  Teuchos::RCP<const std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > getCoords() const
  {
    return coords;
  }

  Teuchos::RCP<std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > getNonconstCoords()
  {
    return coords;
  }
    

private:
  Teuchos::RCP<std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > coords;
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
   void setGlobalIndexer(const Teuchos::RCP<UniqueGlobalIndexer<int,int> > & ugi)
   { globalIndexer_ = ugi; }

   void setLinearObjFactory(const Teuchos::RCP<LinearObjFactory<panzer::Traits> > & lof)
   { linObjFactory_ = lof; }

   Teuchos::RCP<UniqueGlobalIndexer<int,int> > getGlobalIndexer() const
   { return globalIndexer_; }

   Teuchos::RCP<LinearObjFactory<panzer::Traits> > getLinearObjFactory() const
   { return linObjFactory_; }

   template <typename EvalT>
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > build() const
   { return Teuchos::null; }

private:
   Teuchos::RCP<UniqueGlobalIndexer<int,int> > globalIndexer_;
   Teuchos::RCP<LinearObjFactory<panzer::Traits> > linObjFactory_;
};

// declaration so methods are not inlined
template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseAggregator_IPCoordinates_Builder::build<panzer::Traits::Residual>() const;

}

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION
#include "Panzer_ResponseAggregator_IPCoordinatesT.hpp"
#endif

#endif
