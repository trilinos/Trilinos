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

#include "Panzer_WorksetContainer.hpp"

#include "Panzer_WorksetFactoryBase.hpp"

namespace panzer
{


WorksetContainer::
WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory):
  factory_(factory)
{

}

WorksetContainer::
WorksetContainer(const WorksetContainer & wc):
  factory_(wc.getFactory())
{
  // We do not copy the worksets themselves
}

void
WorksetContainer::
setFactory(const Teuchos::RCP<const WorksetFactoryBase> & factory)
{
  // Clear all worksets
  clearWorksets();

  // Now set the factory
  factory_ = factory;
}

Teuchos::RCP<const WorksetFactoryBase>
WorksetContainer::
getFactory() const
{
  return factory_;
}

Teuchos::RCP<std::vector<Workset> >
WorksetContainer::
getWorksets(const WorksetDescriptor & wd)
{
  auto itr = worksets_.find(wd);
  if(itr == worksets_.end()){

    // Make sure the factory exists
    TEUCHOS_TEST_FOR_EXCEPT_MSG(factory_.is_null(), "WorksetContainer::getWorksets : Factory has not been set");

    // Generate the worksets using the factory
    auto worksets = factory_->getWorksets(wd);

    // Set a somewhat unique identifier per workset
    std::size_t hash = std::hash<WorksetDescriptor>()(wd);
    for(auto & workset : *worksets)
      workset.setIdentifier(hash++);

    // Even if worksets is empty, we will still store it
    worksets_[wd] = worksets;

    return worksets;
  }
  return itr->second;
}

void
WorksetContainer::
clearWorksets()
{
  worksets_.clear();
}

}
