// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
#include "PanzerDofMgr_config.hpp"

#include "UnitTest_ConnManager.hpp"

#include <vector>

namespace panzer {
namespace unit_test {

ConnManager::ConnManager(int rank,int procCount)
{
   TEUCHOS_ASSERT(procCount==2);
   TEUCHOS_ASSERT(rank==0 || rank==1);

   procRank_ = rank; 

   if(procRank_==0) {
      elements_["block_0"].push_back(0);
      elements_["block_0"].push_back(1);
      elements_["block_1"].push_back(4);
      elements_["block_1"].push_back(5);
      elements_["block_2"].push_back(2);
      elements_["block_2"].push_back(3);
   }

   if(procRank_==1) {
      elements_["block_0"].push_back(0);
      elements_["block_0"].push_back(1);
      elements_["block_1"].push_back(2);
      elements_["block_1"].push_back(3);
      elements_["block_2"].resize(0);
   }
}

  void GLOBAL_CONN(std::vector<std::vector<panzer::GlobalOrdinal> > & conn,
                        int le,int a,int b,int c,int d)
{ conn[le].push_back(a); conn[le].push_back(b); conn[le].push_back(c); conn[le].push_back(d); }

void ConnManager::buildConnectivity(const FieldPattern & fp) 
{
   if(callback_!=Teuchos::null) 
      callback_->buildConnectivity(fp);

   if(procRank_==0) {
      connectivity_.resize(6);

      GLOBAL_CONN(connectivity_, 0,   0,  1,  6,  5);
      GLOBAL_CONN(connectivity_, 1,   5,  6, 11, 10);
      GLOBAL_CONN(connectivity_, 2,  10, 11, 16, 15);
      GLOBAL_CONN(connectivity_, 3,  11, 12, 17, 16);
      GLOBAL_CONN(connectivity_, 4,   7,  8, 13, 12);
      GLOBAL_CONN(connectivity_, 5,   8,  9, 14, 13);

      return;
   }

   if(procRank_==1) {
      connectivity_.resize(4);

      GLOBAL_CONN(connectivity_, 0,   1,  2,  7,  6);
      GLOBAL_CONN(connectivity_, 1,   6,  7, 12, 11);
      GLOBAL_CONN(connectivity_, 2,   2,  3,  8,  7);
      GLOBAL_CONN(connectivity_, 3,   3,  4,  9,  8);

      return;
   }

   TEUCHOS_ASSERT(false);
}

const ConnManager::GlobalOrdinal * ConnManager::getConnectivity(ConnManager::LocalOrdinal localElmtId) const
{
   return &connectivity_[localElmtId][0];
}

ConnManager::LocalOrdinal ConnManager::getConnectivitySize(ConnManager::LocalOrdinal /* localElmtId */) const
{ return 4; }

std::string ConnManager::getBlockId(ConnManager::LocalOrdinal localElmtId) const
{
   if(procRank_==0) {
      switch(localElmtId) {
      case 0: return "block_0";
      case 1: return "block_0";
      case 2: return "block_2";
      case 3: return "block_2";
      case 4: return "block_1";
      case 5: return "block_1";
      default: break;
      }
   }
   else if(procRank_==1) {
      switch(localElmtId) {
      case 0: return "block_0";
      case 1: return "block_0";
      case 2: return "block_1";
      case 3: return "block_1";
      default: break;
      }
   }
   TEUCHOS_ASSERT(false);
}

std::size_t ConnManager::numElementBlocks() const { return 3; }

void ConnManager::getElementBlockIds(std::vector<std::string> & elementBlockIds) const
{
   elementBlockIds.push_back("block_0");
   elementBlockIds.push_back("block_1");
   elementBlockIds.push_back("block_2");
}

void ConnManager::getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const
{
  const CellTopologyData & myCellData = *shards::getCellTopologyData<shards::Quadrilateral<4> >();
  struct shards::CellTopology my_topo(&myCellData);
  elementBlockTopologies = std::vector<shards::CellTopology>(3,my_topo);
}

const std::vector<ConnManager::LocalOrdinal> & ConnManager::getElementBlock(const std::string & blockID) const
{
   std::map<std::string,std::vector<int> >::const_iterator itr = elements_.find(blockID);
   if(itr!=elements_.end())
      return itr->second;

   TEUCHOS_ASSERT(false);
}

} // end unit test
} // end panzer
