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
   void ghostToGlobal(int p);

   //! Call global to ghost on all the containers
   void globalToGhost(int p);

   //! Call initialize on all containers
   void initialize();

   const_iterator begin() const { return lookupTable_.begin(); }
   const_iterator end() const { return lookupTable_.end(); }

   iterator begin() { return lookupTable_.begin(); }
   iterator end() { return lookupTable_.end(); }

private:
   boost::unordered_map<std::string,Teuchos::RCP<GlobalEvaluationData> > lookupTable_;
};

}

#endif
