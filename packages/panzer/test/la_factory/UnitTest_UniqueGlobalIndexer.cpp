#include "UnitTest_UniqueGlobalIndexer.hpp"

namespace panzer {
namespace unit_test {

UniqueGlobalIndexer::UniqueGlobalIndexer(int rank,int procCount)
   : procRank_(rank)
{
   TEST_FOR_EXCEPTION(procCount!=2,std::runtime_error,"unit_test::UniqueGlobalIndexer runs on only two processors!");
}

int UniqueGlobalIndexer::getFieldNum(const std::string & str) const
{
   if(str=="U") 
      return 0;
   else if(str=="T") 
      return 1;
   else  
      TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << str << "\" in unit_test::UniqueGlobalIndexer, try \'U\' or \'T\'");
}

void UniqueGlobalIndexer::getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
{
   elementBlockIds.clear();
   elementBlockIds.push_back("block_0");
}

const std::vector<short> & UniqueGlobalIndexer::getElementBlock(const std::string & blockId) const
{
   TEST_FOR_EXCEPTION(blockId!="block_0",std::runtime_error,
                      "Can't find block ID \"" << blockId << "\" in unit_test::UniqueGlobalIndexer");

   if(elements_==Teuchos::null) { 
      elements_ = Teuchos::rcp(new std::vector<short>);
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
}

void UniqueGlobalIndexer::getElementGIDs(short localElmtId,std::vector<int> & gids) const
{
   gids.resize(8);

   switch(localElmtId) {
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

const std::vector<int> & UniqueGlobalIndexer::getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
{
   TEST_FOR_EXCEPTION(not ((fieldNum==0 || fieldNum==1) && blockId=="block_0"), std::runtime_error,
                   "unit_test::UniqueGlobalIndexer - Invalid field or block id specified");

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
UniqueGlobalIndexer::getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                                                int subcellDim,int subcellId) const
{
   TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "unit_test::UniqueGlobalIndexer::getGIDFieldOffsets_closure is not implemented yet.");
}

void UniqueGlobalIndexer::getOwnedIndices(std::vector<int> & indices) const
{
   indices.resize(6);
   switch(procRank_) {
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
   }
}

void UniqueGlobalIndexer::getOwnedAndSharedIndices(std::vector<int> & indices) const
{
   indices.resize(8);
   switch(procRank_) {
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
   }
}

void UniqueGlobalIndexer::ownedIndices(const std::vector<int> & indices,std::vector<bool> & isOwned) const
{
   std::vector<int> owned;
   getOwnedIndices(owned);

   isOwned.resize(indices.size(),false);
   for(std::size_t i=0;i<indices.size();i++) 
      isOwned[i] = (std::find(owned.begin(),owned.end(),indices[i])!=owned.end());
}

} // end unit test
} // end panzer
