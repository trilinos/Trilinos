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
   boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator itr = lookupTable_.find(key); 
   if(itr==lookupTable_.end()) {
     std::stringstream ss;
     ss << "Valid keys = ";
     for(const_iterator litr=begin();litr!=end();++litr)
       ss << "\"" << litr->first << "\" ";

     TEUCHOS_TEST_FOR_EXCEPTION(itr==lookupTable_.end(),std::logic_error,
                        "In GlobalEvaluationDataContainer::getDataObject(key) failed to find the data object specified by \""+key+"\"\n   " + ss.str());
   }

   return itr->second; 
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
