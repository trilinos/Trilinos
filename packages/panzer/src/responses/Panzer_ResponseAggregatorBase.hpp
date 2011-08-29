#ifndef __Panzer_ResponseAggregatorBase_hpp__
#define __Panzer_ResponseAggregatorBase_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_Response.hpp"

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
   virtual void registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,Teuchos::RCP<ResponseData<TraitsT> > & data,
                                             const Teuchos::ParameterList & p) const = 0;

   //! perform global reduction on this set of response data
   virtual void globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const = 0;
};

}

#endif
