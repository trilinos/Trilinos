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

#include "UnitTest_UniqueGlobalIndexer.hpp"

namespace panzer {
namespace unit_test {

UniqueGlobalIndexer::UniqueGlobalIndexer(int rank,int procCount)
   : procRank_(rank)
{
   TEUCHOS_TEST_FOR_EXCEPTION(procCount!=2,std::runtime_error,"unit_test::UniqueGlobalIndexer runs on only two processors!");
}

int UniqueGlobalIndexer::getFieldNum(const std::string & str) const
{
   if(str=="TEMPERATURE") 
      return 0;
   else if(str=="ION_TEMPERATURE") 
      return 1;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << str << "\" in unit_test::UniqueGlobalIndexer, try \'U\' or \'T\'");
}

const std::string & UniqueGlobalIndexer::getFieldString(int field) const
{
   static std::string u = "U";
   static std::string t = "T";

   if(field==0)
      return u;
   else if(field==1)
      return t;
   else  
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Can't find field \"" << field << "\" in unit_test::UniqueGlobalIndexer, try \'0\' or \'1\'");
}

void UniqueGlobalIndexer::getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
{
   elementBlockIds.clear();
   elementBlockIds.push_back("block_0");
}

bool UniqueGlobalIndexer::fieldInBlock(const std::string & field, const std::string & block) const
{
   if(block!="block_0") 
      return false;
 
   if(field=="TEMPERATURE" || field=="ION_TEMPERATURE")
      return true;

   return false;
}

const std::vector<int> & UniqueGlobalIndexer::getElementBlock(const std::string & blockId) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(blockId!="block_0",std::runtime_error,
                      "Can't find block ID \"" << blockId << "\" in unit_test::UniqueGlobalIndexer");

   if(elements_==Teuchos::null) { 
      elements_ = Teuchos::rcp(new std::vector<int>);
      switch(procRank_) {
      case 0:
         elements_->push_back(0);
         break;
      case 1:
         elements_->push_back(0);
         break;
      default:
         TEUCHOS_ASSERT(false);
      }
   }

   return *elements_;
}

void UniqueGlobalIndexer::getElementGIDs(int localElmtId,std::vector<int> & gids,const std::string & blockId) const
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

const std::vector<int> & UniqueGlobalIndexer::getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(not ((fieldNum==0 || fieldNum==1) && blockId=="block_0"), std::runtime_error,
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
   typedef std::pair<std::vector<int>,std::vector<int> > VectorPair;

   // TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
   //                    "unit_test::UniqueGlobalIndexer::getGIDFieldOffsets_closure is not implemented yet.");
   TEUCHOS_TEST_FOR_EXCEPTION(not ((fieldNum==0 || fieldNum==1) && blockId=="block_0"), std::runtime_error,
                   "unit_test::UniqueGlobalIndexer - Invalid field or block id specified");

   if(field0_closures_==Teuchos::null || field1_closures_==Teuchos::null) {
      field0_closures_ = Teuchos::rcp(new std::map<std::pair<int,int>,VectorPair>); 

      {
         VectorPair vp;

         vp = std::make_pair(std::vector<int>(2),std::vector<int>(2));
         vp.first[0] = 0; vp.first[1] = 2;
         vp.second[0] = 0; vp.second[1] = 1;
         (*field0_closures_)[std::make_pair(1,0)] = vp;

         vp = std::make_pair(std::vector<int>(2),std::vector<int>(2));
         vp.first[0] = 2; vp.first[1] = 4;
         vp.second[0] = 1; vp.second[1] = 2;
         (*field0_closures_)[std::make_pair(1,1)] = vp;
      }

      field1_closures_ = Teuchos::rcp(new std::map<std::pair<int,int>,VectorPair>); 
      {
         VectorPair vp;

         vp = std::make_pair(std::vector<int>(2),std::vector<int>(2));
         vp.first[0] = 1; vp.first[1] = 3;
         vp.second[0] = 0; vp.second[1] = 1;
         (*field1_closures_)[std::make_pair(1,0)] = vp;

         vp = std::make_pair(std::vector<int>(2),std::vector<int>(2));
         vp.first[0] = 3; vp.first[1] = 5;
         vp.second[0] = 1; vp.second[1] = 2;
         (*field1_closures_)[std::make_pair(1,1)] = vp;
      }
   }

   std::stringstream ss;
   ss << "(subcellDim,subcellId) = (" << subcellDim << "," << subcellId << ") on fieldNum = " << fieldNum;

   if(fieldNum==0) {
      std::map<std::pair<int,int>, VectorPair>::const_iterator itr 
            = field0_closures_->find(std::make_pair(subcellDim,subcellId));
      TEUCHOS_TEST_FOR_EXCEPTION(itr==field0_closures_->end(),std::runtime_error,
                         "unit_test::UniqueGlobalIndexer::getGIDFieldOffsets_closure has not added "+ss.str()+" yet");

      return itr->second;
   }
   else {
      std::map<std::pair<int,int>, VectorPair>::const_iterator itr 
            = field1_closures_->find(std::make_pair(subcellDim,subcellId));
      TEUCHOS_TEST_FOR_EXCEPTION(itr==field1_closures_->end(),std::runtime_error,
                         "unit_test::UniqueGlobalIndexer::getGIDFieldOffsets_closure has not added "+ss.str()+" yet");

      return itr->second;
   } 
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

/** Get field numbers associated with a particular element block.
  */
const std::vector<int> & UniqueGlobalIndexer::getBlockFieldNumbers(const std::string & blockId) const
{
   static std::vector<int> fieldNums;
   if(fieldNums.size()==0) {
      fieldNums.resize(8);

      fieldNums[0] = getFieldNum("TEMPERATURE"); fieldNums[1] = getFieldNum("ION_TEMPERATURE");
      fieldNums[2] = getFieldNum("TEMPERATURE"); fieldNums[3] = getFieldNum("ION_TEMPERATURE");
      fieldNums[4] = getFieldNum("TEMPERATURE"); fieldNums[5] = getFieldNum("ION_TEMPERATURE");
      fieldNums[6] = getFieldNum("TEMPERATURE"); fieldNums[7] = getFieldNum("ION_TEMPERATURE");
   }

   return fieldNums;
}

void UniqueGlobalIndexer::getCoordinates(int localElementId,Intrepid::FieldContainer<double> & vertices)
{
   vertices.resize(1,4,2);
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

} // end unit test
} // end panzer
