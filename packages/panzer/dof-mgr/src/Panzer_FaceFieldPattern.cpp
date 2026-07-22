// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_FaceFieldPattern.hpp"

#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

FaceFieldPattern::FaceFieldPattern(const shards::CellTopology & ct)
{
   setCellTopology(ct);
}

void FaceFieldPattern::setCellTopology(const shards::CellTopology & ct)
{
   cellTopo_ = ct;

   // allocate the space and setup the indices
   FaceIndices_.clear();
   if ( cellTopo_.getDimension() == 3 ) 
     FaceIndices_.resize(cellTopo_.getFaceCount());
   else
     FaceIndices_.resize(1);
   for(std::size_t n=0;n<FaceIndices_.size();n++)
      FaceIndices_[n].push_back(n);
}

int FaceFieldPattern::getSubcellCount(int dim) const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getSubcellCount(dim);
}

const std::vector<int> & FaceFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(dim==2)
      return FaceIndices_[cellIndex];
   
   // only Faces
   return empty_;
}

void FaceFieldPattern::getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int>& /* indices */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                      "FaceFieldPattern::getSubcellClosureIndices should not be called");
}

int FaceFieldPattern::getDimension() const
{
   const shards::CellTopology ct = getCellTopology();
   return ct.getDimension();
}

}
