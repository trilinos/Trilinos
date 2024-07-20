// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_NodalFieldPattern.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

NodalFieldPattern::NodalFieldPattern(const shards::CellTopology & ct)
{
   setCellTopology(ct);
}

void NodalFieldPattern::setCellTopology(const shards::CellTopology & ct)
{
   cellTopo_ = ct;

   // allocate the space and setup the indices
   nodeIndices_.clear();
   nodeIndices_.resize(cellTopo_.getNodeCount());
   for(std::size_t n=0;n<nodeIndices_.size();n++)
      nodeIndices_[n].push_back(n);
}

int NodalFieldPattern::getSubcellCount(int dim) const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getSubcellCount(dim);
}

const std::vector<int> & NodalFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(dim==0)
      return nodeIndices_[cellIndex];
   
   // only nodes
   return empty_;
}

void NodalFieldPattern::getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int>& /* indices */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                      "NodalFieldPattern::getSubcellClosureIndices should not be called"); 
}

int NodalFieldPattern::getDimension() const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getDimension();
}

}
