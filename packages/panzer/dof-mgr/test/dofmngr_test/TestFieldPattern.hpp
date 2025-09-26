// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TestFieldPattern_hpp__
#define __TestFieldPattern_hpp__

#include "Panzer_FieldPattern.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

class TestFieldPattern : public FieldPattern {
public:
   TestFieldPattern() {}

  Teuchos::RCP<panzer::FieldPattern> clone() const override
  {return Teuchos::rcp(new TestFieldPattern(*this));}
  
   /* This function has no functionality in this case.
    * If called it will throw an assertion failure
    */
   virtual void getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int> & /* indices */) const override
   { TEUCHOS_ASSERT(false); }

   virtual int getSubcellCount(int dim) const override
   {  return subcellIndices[dim].size(); }

   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const override
   {  return subcellIndices[dim][cellIndex]; }

   virtual int getDimension() const override
   { return subcellIndices.size()-1; }

   std::vector<std::vector<int> > & operator[](int v)
   { return subcellIndices[v]; } 

   virtual shards::CellTopology getCellTopology() const override
   { return cellTopo; }

public:
   std::vector<std::vector<std::vector<int> > > subcellIndices;
   shards::CellTopology cellTopo;
};

}

#endif
