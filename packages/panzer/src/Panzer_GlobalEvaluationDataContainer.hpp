#ifndef __Panzer_GlobalEvaluationDataContainer_hpp__
#define __Panzer_GlobalEvaluationDataContainer_hpp__

#include "Teuchos_RCP.hpp"

#include <boost/unordered_map.hpp>

#include "Panzer_GlobalEvaluationData.hpp"

namespace panzer {

class GlobalEvaluationDataContainer {
public:
   typedef boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator const_iterator;
   typedef boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> >::iterator iterator;

   /** Add a data object to be used in evaluation loop.
     */
   void addDataObject(const std::string & key,
                      const Teuchos::RCP<GlobalEvaluationData> & ged);

   /** Does this container have a match to a certain key.
     */
   bool containsDataObject(const std::string & key) const;

   /** Get the data object associated with the key.
     */
   Teuchos::RCP<GlobalEvaluationData> getDataObject(const std::string & key) const;

   //! Call ghost to global on all the containers
   void ghostToGlobal();

   //! Call global to ghost on all the containers
   void globalToGhost();

   const_iterator begin() const { return lookupTable_.begin(); }
   const_iterator end() const { return lookupTable_.end(); }

   iterator begin() { return lookupTable_.begin(); }
   iterator end() { return lookupTable_.end(); }

private:
   boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> > lookupTable_;
};

}

#endif
