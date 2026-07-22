// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_EdgeFieldPattern.hpp"

#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

EdgeFieldPattern::EdgeFieldPattern(const shards::CellTopology & ct)
{
   setCellTopology(ct);
}

void EdgeFieldPattern::setCellTopology(const shards::CellTopology & ct)
{
   cellTopo_ = ct;

   // allocate the space and setup the indices
   edgeIndices_.clear();
   // edgeIndices_.resize(cellTopo_.getEdgeCount());
   edgeIndices_.resize(cellTopo_.getSubcellCount(1));
   for(std::size_t n=0;n<edgeIndices_.size();n++)
      edgeIndices_[n].push_back(n);
}

int EdgeFieldPattern::getSubcellCount(int dim) const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getSubcellCount(dim);
}

const std::vector<int> & EdgeFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(dim==1)
      return edgeIndices_[cellIndex];
   
   // only edges
   return empty_;
}

void EdgeFieldPattern::getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int>& /* indices */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                      "EdgeFieldPattern::getSubcellClosureIndices should not be called"); 
}

int EdgeFieldPattern::getDimension() const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getDimension();
}

}
