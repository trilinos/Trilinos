#ifndef __Panzer_GlobalEvaluationDataContainer_hpp__
#define __Panzer_GlobalEvaluationDataContainer_hpp__

#include "Teuchos_RCP.hpp"

#include <boost/unordered_map.hpp>

namespace panzer {

class GlobalEvaluationData {
public:
   virtual ~GlobalEvaluationData() = 0;
};

class GlobalEvaluationDataContainer {
public:
   /** Add a data object to be used in evaluation loop.
     */
   void addDataObject(const std::string & key,
                      const Teuchos::RCP<GlobalEvaluationData> & ged);

   /** Does this containe have a match to a certain key.
     */
   bool containsDataObject(const std::string & key) const;

   /** Get the data object associated with the key.
     */
   Teuchos::RCP<GlobalEvaluationData> getDataObject(const std::string & key) const;

private:
   boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> > lookupTable_;
};

}

#endif
