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

#include "Panzer_config.hpp"

#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

using Teuchos::RCP;

namespace panzer {
namespace orientation_helpers {

void computePatternEdgeIndices(const FieldPattern & pattern,std::vector<std::pair<int,int> > & edgeIndices)
{
   unsigned dim = 1;
   shards::CellTopology cellTopo = pattern.getCellTopology();
   for(unsigned e=0;e<cellTopo.getEdgeCount();e++) {
      // get local vertex ids for a this edge
      unsigned local_v0 = cellTopo.getNodeMap(dim,e,0);
      unsigned local_v1 = cellTopo.getNodeMap(dim,e,1);

      // get sub cell indices for geometric pattern
      const std::vector<int> & v0_indices = pattern.getSubcellIndices(0,local_v0);
      const std::vector<int> & v1_indices = pattern.getSubcellIndices(0,local_v1);

      TEUCHOS_ASSERT(v0_indices.size()>0); // there must be a node
      TEUCHOS_ASSERT(v1_indices.size()>0); // there must be a node

      // take the first index on each vertex and make a edge lookup
      edgeIndices.push_back(std::make_pair(v0_indices[0],v1_indices[0]));
   }
}

} // end orientation_helpers
} // end panzer
