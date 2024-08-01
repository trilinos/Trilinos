// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_GlobalEvaluationDataContainer.hpp"

namespace panzer {

GlobalEvaluationData::~GlobalEvaluationData() {}

/** Add a data object to be used in evaluation loop.
  */
void GlobalEvaluationDataContainer::addDataObject(const std::string & key,
                                                  const Teuchos::RCP<GlobalEvaluationData> & ged)
{
   lookupTable_[key] = ged;
}

/** Does this containe have a match to a certain key.
  */
bool GlobalEvaluationDataContainer::containsDataObject(const std::string & key) const
{
   return lookupTable_.find(key)!=lookupTable_.end();
}

/** Get the data object associated with the key.
  */
Teuchos::RCP<GlobalEvaluationData> GlobalEvaluationDataContainer::getDataObject(const std::string & key) const
{
   std::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator itr;
   for(itr=begin();itr!=end();++itr) {
     if(itr->first==key) 
       return itr->second;
   }

   {
     std::stringstream ss;
     ss << "Valid keys = ";
     for(const_iterator litr=begin();litr!=end();++litr)
       ss << "\"" << litr->first << "\" ";

     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                        "In GlobalEvaluationDataContainer::getDataObject(key) failed to find the data object specified by \""+key+"\"\n   " + ss.str());
   }

   return Teuchos::null;
/*
   std::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator itr = lookupTable_.find(key); 
   
   if(itr==lookupTable_.end()) {
     std::stringstream ss;
     ss << "Valid keys = ";
     for(const_iterator litr=begin();litr!=end();++litr)
       ss << "\"" << litr->first << "\" ";

     TEUCHOS_TEST_FOR_EXCEPTION(itr==lookupTable_.end(),std::logic_error,
                        "In GlobalEvaluationDataContainer::getDataObject(key) failed to find the data object specified by \""+key+"\"\n   " + ss.str());
   }

   return itr->second; 
*/
}

//! Call ghost to global on all the containers
void GlobalEvaluationDataContainer::ghostToGlobal(int p)
{
   for(iterator itr=begin();itr!=end();++itr)
     itr->second->ghostToGlobal(p);
}

//! Call global to ghost on all the containers
void GlobalEvaluationDataContainer::globalToGhost(int p)
{
   for(iterator itr=begin();itr!=end();++itr)
     itr->second->globalToGhost(p);
}

//! Call global to ghost on all the containers
void GlobalEvaluationDataContainer::initialize()
{
   for(iterator itr=begin();itr!=end();++itr)
     itr->second->initializeData();
}

}
