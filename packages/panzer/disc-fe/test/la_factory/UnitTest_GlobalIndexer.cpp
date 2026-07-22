// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "UnitTest_GlobalIndexer.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {
namespace unit_test {

// GlobalIndexer
/////////////////////////////////////////////////////////////////////////////////////////

GlobalIndexer::GlobalIndexer(int rank,int procCount)
   : procRank_(rank)
{
   TEUCHOS_TEST_FOR_EXCEPTION(procCount!=2,std::runtime_error,"unit_test::GlobalIndexer runs on only two processors!");
}

int GlobalIndexer::getFieldNum(const std::string & str) const
{
   if(str=="U") 
      return 0;
   else if(str=="T") 
      return 1;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << str << "\" in unit_test::GlobalIndexer, try \'U\' or \'T\'");
}

const std::string & GlobalIndexer::getFieldString(int field) const
{
   static std::string u = "U";
   static std::string t = "T";

   if(field==0)
      return u;
   else if(field==1)
      return t;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << field << "\" in unit_test::GlobalIndexer, try \'0\' or \'1\'");
}

void GlobalIndexer::getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
{
   elementBlockIds.clear();
   elementBlockIds.push_back("block_0");
}

bool GlobalIndexer::fieldInBlock(const std::string & field, const std::string & block) const
{
   if(block!="block_0") 
      return false;
 
   if(field=="U" || field=="T")
      return true;

   return false;
}

const std::vector<panzer::LocalOrdinal> & GlobalIndexer::getElementBlock(const std::string & blockId) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(blockId!="block_0",std::runtime_error,
                      "Can't find block ID \"" << blockId << "\" in unit_test::GlobalIndexer");

   if(elements_==Teuchos::null) { 
      elements_ = Teuchos::rcp(new std::vector<panzer::LocalOrdinal>);
      switch(procRank_) {
      case 0:
         elements_->push_back(0);
         break;
      case 1:
         elements_->push_back(1);
         break;
      default:
         TEUCHOS_ASSERT(false);
      }
   }

   return *elements_;
}

void GlobalIndexer::getElementGIDs(panzer::LocalOrdinal /* localElmtId */, std::vector<panzer::GlobalOrdinal>& gids, const std::string& /* blockId */) const
{
   gids.resize(8);

   switch(procRank_) {
   case 0:
      gids[0] = 0; gids[1] = 1; 
      gids[2] = 2; gids[3] = 3;
      gids[4] = 4; gids[5] = 5; 
      gids[6] = 6; gids[7] = 7;
      break;
   case 1:
      gids[0] = 2; gids[1] = 3; 
      gids[2] = 8; gids[3] = 9;
      gids[4] = 10; gids[5] = 11; 
      gids[6] = 4; gids[7] = 5;
      break;
   default:
      break;
   }
}

const std::vector<int> & GlobalIndexer::getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(not ((fieldNum==0 || fieldNum==1) && blockId=="block_0"), std::runtime_error,
                   "unit_test::GlobalIndexer - Invalid field or block id specified");

   if(field0Offset_==Teuchos::null || field1Offset_==Teuchos::null) {
      field0Offset_ = Teuchos::rcp(new std::vector<int>(4)); 
      (*field0Offset_)[0] = 0;
      (*field0Offset_)[1] = 2;
      (*field0Offset_)[2] = 4;
      (*field0Offset_)[3] = 6;

      field1Offset_ = Teuchos::rcp(new std::vector<int>(4)); 
      (*field1Offset_)[0] = 1;
      (*field1Offset_)[1] = 3;
      (*field1Offset_)[2] = 5;
      (*field1Offset_)[3] = 7;
   }

   if(fieldNum==0) return *field0Offset_;
   else return *field1Offset_;
}

const std::pair<std::vector<int>,std::vector<int> > & 
GlobalIndexer::getGIDFieldOffsets_closure(const std::string& /* blockId */, int /* fieldNum */,
                                                int /* subcellDim */, int /* subcellId */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "unit_test::GlobalIndexer::getGIDFieldOffsets_closure is not implemented yet.");
}

///////////////////////////////////////////////////////////////////////////////
//
//  getOwnedIndices()
//
///////////////////////////////////////////////////////////////////////////////
void
GlobalIndexer::
getOwnedIndices(
  std::vector<panzer::GlobalOrdinal>& indices) const
{
  indices.resize(6);
  switch (procRank_)
  {
    case 0:
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      indices[3] = 3;
      indices[4] = 6;
      indices[5] = 7;
      break;
    case 1:
      indices[0] = 8;
      indices[1] = 9;
      indices[2] = 10;
      indices[3] = 11;
      indices[4] = 4;
      indices[5] = 5;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getOwnedIndices()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedIndices()
//
///////////////////////////////////////////////////////////////////////////////
void
GlobalIndexer::
getGhostedIndices(
  std::vector<panzer::GlobalOrdinal>& indices) const
{
  indices.resize(2);
  switch (procRank_)
  {
    case 0:
      indices[0] = 4;
      indices[1] = 5;
      break;
    case 1:
      indices[0] = 2;
      indices[1] = 3;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getGhostedIndices()

///////////////////////////////////////////////////////////////////////////////
//
//  getOwnedAndGhostedIndices()
//
///////////////////////////////////////////////////////////////////////////////
void
GlobalIndexer::
getOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
{
  indices.resize(8);
  switch (procRank_)
  {
    case 0:
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      indices[3] = 3;
      indices[4] = 4;
      indices[5] = 5;
      indices[6] = 6;
      indices[7] = 7;
      break;
    case 1:
      indices[0] = 2;
      indices[1] = 3;
      indices[2] = 8;
      indices[3] = 9;
      indices[4] = 10;
      indices[5] = 11;
      indices[6] = 4;
      indices[7] = 5;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getOwnedAndGhostedIndices()

///////////////////////////////////////////////////////////////////////////////
void GlobalIndexer::getElementGIDsAsInt(panzer::LocalOrdinal /* localElmtId */, std::vector<int>& gids, const std::string& /* blockId */) const
{
   gids.resize(8);

   switch(procRank_) {
   case 0:
      gids[0] = 0; gids[1] = 1; 
      gids[2] = 2; gids[3] = 3;
      gids[4] = 4; gids[5] = 5; 
      gids[6] = 6; gids[7] = 7;
      break;
   case 1:
      gids[0] = 2; gids[1] = 3; 
      gids[2] = 8; gids[3] = 9;
      gids[4] = 10; gids[5] = 11; 
      gids[6] = 4; gids[7] = 5;
      break;
   default:
      break;
   }
}
///////////////////////////////////////////////////////////////////////////////
void GlobalIndexer::
getOwnedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(6);
  switch (procRank_)
  {
    case 0:
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      indices[3] = 3;
      indices[4] = 6;
      indices[5] = 7;
      break;
    case 1:
      indices[0] = 8;
      indices[1] = 9;
      indices[2] = 10;
      indices[3] = 11;
      indices[4] = 4;
      indices[5] = 5;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getOwnedIndices()

///////////////////////////////////////////////////////////////////////////////
void GlobalIndexer::
getGhostedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(2);
  switch (procRank_)
  {
    case 0:
      indices[0] = 4;
      indices[1] = 5;
      break;
    case 1:
      indices[0] = 2;
      indices[1] = 3;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getGhostedIndices()

///////////////////////////////////////////////////////////////////////////////
void GlobalIndexer::
getOwnedAndGhostedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(8);
  switch (procRank_)
  {
    case 0:
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      indices[3] = 3;
      indices[4] = 4;
      indices[5] = 5;
      indices[6] = 6;
      indices[7] = 7;
      break;
    case 1:
      indices[0] = 2;
      indices[1] = 3;
      indices[2] = 8;
      indices[3] = 9;
      indices[4] = 10;
      indices[5] = 11;
      indices[6] = 4;
      indices[7] = 5;
      break;
    default:
      TEUCHOS_ASSERT(false);
  } // end switch (procRank_)
} // end of getOwnedAndGhostedIndices()


///////////////////////////////////////////////////////////////////////////////
//
//  getNumOwned()
//
///////////////////////////////////////////////////////////////////////////////
int
GlobalIndexer::
getNumOwned() const
{
  return 6;
} // end of getNumOwned()

///////////////////////////////////////////////////////////////////////////////
//
//  getNumGhosted()
//
///////////////////////////////////////////////////////////////////////////////
int
GlobalIndexer::
getNumGhosted() const
{
  return 2;
} // end of getNumGhosted()

///////////////////////////////////////////////////////////////////////////////
//
//  getNumOwnedAndGhosted()
//
///////////////////////////////////////////////////////////////////////////////
int
GlobalIndexer::
getNumOwnedAndGhosted() const
{
  return 8;
} // end of getNumOwnedAndGhosted()

void GlobalIndexer::ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const
{
   std::vector<panzer::GlobalOrdinal> owned;
   getOwnedIndices(owned);

   isOwned.resize(indices.size(),false);
   for(std::size_t i=0;i<indices.size();i++) 
      isOwned[i] = (std::find(owned.begin(),owned.end(),indices[i])!=owned.end());
}

/** Get field numbers associated with a particular element block.
  */
const std::vector<int> & GlobalIndexer::getBlockFieldNumbers(const std::string & /* blockId */) const
{
   static std::vector<int> fieldNums;
   if(fieldNums.size()==0) {
      fieldNums.resize(8);

      fieldNums[0] = getFieldNum("U"); fieldNums[1] = getFieldNum("T");
      fieldNums[2] = getFieldNum("U"); fieldNums[3] = getFieldNum("T");
      fieldNums[4] = getFieldNum("U"); fieldNums[5] = getFieldNum("T");
      fieldNums[6] = getFieldNum("U"); fieldNums[7] = getFieldNum("T");
   }

   return fieldNums;
}

void GlobalIndexer::getCoordinates(panzer::LocalOrdinal /* localElementId */, Kokkos::DynRankView<double,PHX::Device>& vertices)
{
  vertices = Kokkos::DynRankView<double,PHX::Device>("vertices",1,4,2);
   switch(procRank_) {
   case 0:
      vertices(0,0,0) = 0.0; vertices(0,0,1) = 0.0;
      vertices(0,1,0) = 1.0; vertices(0,1,1) = 0.0;
      vertices(0,2,0) = 1.0; vertices(0,2,1) = 1.0;
      vertices(0,3,0) = 0.0; vertices(0,3,1) = 1.0;
      break;
   case 1:
      vertices(0,0,0) = 1.0; vertices(0,0,1) = 0.0;
      vertices(0,1,0) = 2.0; vertices(0,1,1) = 0.0;
      vertices(0,2,0) = 2.0; vertices(0,2,1) = 1.0;
      vertices(0,3,0) = 1.0; vertices(0,3,1) = 1.0;
      break;
   default:
      TEUCHOS_ASSERT(false); // fail!
   };
}

int GlobalIndexer::
getElementBlockGIDCount(const std::string &) const
{
  return 8;
}

int GlobalIndexer::
getElementBlockGIDCount(const std::size_t &) const
{
  return 8;
}

// GlobalIndexer_Element
/////////////////////////////////////////////////////////////////////////////////////////

GlobalIndexer_Element::GlobalIndexer_Element(int rank,int procCount)
   : procRank_(rank)
{
   TEUCHOS_TEST_FOR_EXCEPTION(procCount!=2,std::runtime_error,"unit_test::GlobalIndexer_Element runs on only two processors!");
}

int GlobalIndexer_Element::getFieldNum(const std::string & str) const
{
   if(str=="U") 
      return 0;
   else if(str=="T") 
      return 1;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << str << "\" in unit_test::GlobalIndexer_Element, try \'U\' or \'T\'");
}

const std::string & GlobalIndexer_Element::getFieldString(int field) const
{
   static std::string u = "U";
   static std::string t = "T";

   if(field==0)
      return u;
   else if(field==1)
      return t;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << field << "\" in unit_test::GlobalIndexer, try \'0\' or \'1\'");
}

void GlobalIndexer_Element::getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
{
   elementBlockIds.clear();
   elementBlockIds.push_back("block_0");
}

bool GlobalIndexer_Element::fieldInBlock(const std::string & field, const std::string & block) const
{
   if(block!="block_0") 
      return false;
 
   if(field=="U" || field=="T")
      return true;

   return false;
}

const std::vector<panzer::LocalOrdinal> & GlobalIndexer_Element::getElementBlock(const std::string & blockId) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(blockId!="block_0",std::runtime_error,
                      "Can't find block ID \"" << blockId << "\" in unit_test::GlobalIndexer_Element");

   if(elements_==Teuchos::null) { 
      elements_ = Teuchos::rcp(new std::vector<panzer::LocalOrdinal>);
      switch(procRank_) {
      case 0:
         elements_->push_back(0);
         break;
      case 1:
         elements_->push_back(1);
         break;
      default:
         TEUCHOS_ASSERT(false);
      }
   }

   return *elements_;
}

void GlobalIndexer_Element::getElementGIDs(panzer::LocalOrdinal /* localElmtId */, std::vector<panzer::GlobalOrdinal>& gids, const std::string& /* blockId */) const
{
   gids.resize(2);

   switch(procRank_) {
   case 0:
      gids[0] = 0; 
      gids[1] = 1;
      break;
   case 1:
      gids[0] = 2; 
      gids[1] = 3;
      break;
   default:
      break;
   }
}

const std::vector<int> & GlobalIndexer_Element::getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(not ((fieldNum==0 || fieldNum==1) && blockId=="block_0"), std::runtime_error,
                   "unit_test::GlobalIndexer_Element - Invalid field or block id specified");

   if(field0Offset_==Teuchos::null || field1Offset_==Teuchos::null) {
      field0Offset_ = Teuchos::rcp(new std::vector<int>(1)); 
      (*field0Offset_)[0] = 0;

      field1Offset_ = Teuchos::rcp(new std::vector<int>(1)); 
      (*field1Offset_)[0] = 1;
   }

   if(fieldNum==0) return *field0Offset_;
   else return *field1Offset_;
}

const std::pair<std::vector<int>,std::vector<int> > & 
GlobalIndexer_Element::getGIDFieldOffsets_closure(const std::string& /* blockId */, int /* fieldNum */,
                                                int /* subcellDim */, int /* subcellId */) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "unit_test::GlobalIndexer_Element::getGIDFieldOffsets_closure is not implemented yet.");
}

void GlobalIndexer_Element::getOwnedIndices(std::vector<panzer::GlobalOrdinal> & indices) const
{
   indices.resize(2);
   switch(procRank_) {
   case 0:
      indices[0] = 0;
      indices[1] = 1;
      break;
   case 1:
      indices[0] = 2;
      indices[1] = 3;
      break;
   default:
      TEUCHOS_ASSERT(false);
   }
}

void GlobalIndexer_Element::getOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal> & indices) const
{
   indices.resize(2);
   switch(procRank_) {
   case 0:
      indices[0] = 0;
      indices[1] = 1;
      break;
   case 1:
      indices[0] = 2;
      indices[1] = 3;
      break;
   default:
      TEUCHOS_ASSERT(false);
   }
}

void GlobalIndexer_Element::ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const
{
   std::vector<panzer::GlobalOrdinal> owned;
   getOwnedIndices(owned);

   isOwned.resize(indices.size(),false);
   for(std::size_t i=0;i<indices.size();i++) 
      isOwned[i] = (std::find(owned.begin(),owned.end(),indices[i])!=owned.end());
}

/** Get field numbers associated with a particular element block.
  */
const std::vector<int> & GlobalIndexer_Element::getBlockFieldNumbers(const std::string & /* blockId */) const
{
   static std::vector<int> fieldNums;
   if(fieldNums.size()==0) {
      fieldNums.resize(2);

      fieldNums[0] = getFieldNum("U"); fieldNums[1] = getFieldNum("T");
   }

   return fieldNums;
}

void GlobalIndexer_Element::getCoordinates(panzer::LocalOrdinal /* localElementId */, Kokkos::DynRankView<double,PHX::Device>& vertices)
{
   vertices = Kokkos::DynRankView<double,PHX::Device>("vertices",1,1,2);
   switch(procRank_) {
   case 0:
      vertices(0,0,0) = 0.5; vertices(0,0,1) = 0.5;
      break;
   case 1:
      vertices(0,0,0) = 1.5; vertices(0,0,1) = 0.5;
      break;
   default:
      TEUCHOS_ASSERT(false); // fail!
   };
}

int GlobalIndexer_Element::
getElementBlockGIDCount(const std::string &) const
{
  return 2;
}

int GlobalIndexer_Element::
getElementBlockGIDCount(const std::size_t &) const
{
  return 2;
}

// BlockGlobalIndexer
/////////////////////////////////////////////////////////////////////////////////////////

BlockGlobalIndexer::BlockGlobalIndexer(int /* blocks */, int rank, int procCount)
   : procRank_(rank)
{
   TEUCHOS_TEST_FOR_EXCEPTION(procCount!=2,std::runtime_error,"unit_test::BlockGlobalIndexer runs on only two processors!");
}

void BlockGlobalIndexer::getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
{
   elementBlockIds.clear();
   elementBlockIds.push_back("block_0");
}

const std::vector<panzer::LocalOrdinal> & BlockGlobalIndexer::getElementBlock(const std::string & blockId) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(blockId!="block_0",std::runtime_error,
                      "Can't find block ID \"" << blockId << "\" in unit_test::BlockGlobalIndexer");

   if(elements_==Teuchos::null) { 
      elements_ = Teuchos::rcp(new std::vector<panzer::LocalOrdinal>);
      switch(procRank_) {
      case 0:
         elements_->push_back(0);
         break;
      case 1:
         elements_->push_back(1);
         break;
      default:
         TEUCHOS_ASSERT(false);
      }
   }

   return *elements_;
}

} // end unit test
} // end panzer
