// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_ElemFieldPattern.hpp"

#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

ElemFieldPattern::ElemFieldPattern(const shards::CellTopology & ct)
{
   setCellTopology(ct);
}

void ElemFieldPattern::setCellTopology(const shards::CellTopology & ct)
{
   cellTopo_ = ct;

   // allocate the space and setup the indices
   ElemIndices_.clear();
   ElemIndices_.resize(1);
   for(std::size_t n=0;n<ElemIndices_.size();n++)
      ElemIndices_[n].push_back(n);
}

int ElemFieldPattern::getSubcellCount(int dim) const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getSubcellCount(dim);
}

const std::vector<int> & ElemFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(dim==getDimension())
      return ElemIndices_[cellIndex];
   
   // only Elems
   return empty_;
}

void ElemFieldPattern::getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int>& /* indices */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                      "ElemFieldPattern::getSubcellClosureIndices should not be called");
}

int ElemFieldPattern::getDimension() const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getDimension();
}

}
