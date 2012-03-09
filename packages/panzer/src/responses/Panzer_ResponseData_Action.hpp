#ifndef __Panzer_ResponseData_Action_hpp__
#define __Panzer_ResponseData_Action_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_ResponseAggregatorBase.hpp"

namespace panzer {

/** A data object that is really only useful as a place holder
  * for a response that is primarily action (not data). This would
  * be the case for instance when writing data to a mesh, or file.
  */
template <typename TraitsT>
class ResponseData_Action : public ResponseDataDefault<TraitsT> {
public:
   ResponseData_Action() {}
   ResponseData_Action(const ResponseData_Action &) {}

   virtual ~ResponseData_Action() {}

   //! Allocate and initialize required storage based on the fields passed in
   virtual void allocateAndInitializeData(const std::vector<std::string> & fields) { }

   /** Reinitialize the data based on the fields original requested by
     * <code>allocateAndInitializeData</code>.
     */
   virtual void reinitializeData() {} 

   /** Fill a response object with data from a particular field.
     */
   virtual void fillResponse(const std::string & field,Response<TraitsT> & response) const {}
};

}

#endif
