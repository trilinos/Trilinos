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

#include "Panzer_FieldAggPattern.hpp"

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {
 
FieldAggPattern::FieldAggPattern() 
{
}

FieldAggPattern::FieldAggPattern(std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns,
                                 const Teuchos::RCP<const FieldPattern> & geomAggPattern)
{
   buildPattern(patterns,geomAggPattern);
}

Teuchos::RCP<const FieldPattern> FieldAggPattern::getGeometricAggFieldPattern() const
{
   return geomAggPattern_;
}

void FieldAggPattern::buildPattern(const std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns,
                                   const Teuchos::RCP<const FieldPattern> & geomAggPattern)
{
   TEUCHOS_ASSERT(patterns.size()>0);

   // build geometric information
   if(geomAggPattern==Teuchos::null) {
      std::vector<FPPtr> patternVec;

      // convert map into vector
      std::vector<std::pair<int,FPPtr> >::const_iterator itr;
      for(itr=patterns.begin();itr!=patterns.end();++itr)
         patternVec.push_back(itr->second);

      // build geometric aggregate field pattern
      geomAggPattern_ = rcp(new GeometricAggFieldPattern(patternVec));
   }
   else
      geomAggPattern_ = geomAggPattern;

   patterns_ = patterns;

   buildFieldIdToPatternIdx(); // builds look up from fieldId to index into patterns_ vector
   buildFieldIdsVector();      // meat of matter: build field pattern information
   buildFieldPatternData();    // this hanldes the "getSubcell*" information
}

int FieldAggPattern::getDimension() const
{
   FPPtr geomPattern = getGeometricAggFieldPattern();
   TEUCHOS_TEST_FOR_EXCEPTION(geomPattern==Teuchos::null,std::logic_error,
                      "Geometric field pattern not yet set, call buildPatterns first");

   return geomPattern->getDimension();
}

shards::CellTopology FieldAggPattern::getCellTopology() const
{
   FPPtr geomPattern = getGeometricAggFieldPattern();
   TEUCHOS_TEST_FOR_EXCEPTION(geomPattern==Teuchos::null,std::logic_error,
                      "Geometric field pattern not yet set, call buildPatterns first");

   return geomPattern->getCellTopology();
}

int FieldAggPattern::getSubcellCount(int dimension) const
{
   return patternData_[dimension].size();
}

const std::vector<int> & FieldAggPattern::getSubcellIndices(int dimension, int subcell) const
{
   return patternData_[dimension][subcell];
}

void FieldAggPattern::getSubcellClosureIndices(int, int, std::vector<int> &) const 
{ 
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                      "FieldAggPattern::getSubcellClosureIndices should not be called"); 
}

void FieldAggPattern::print(std::ostream & os) const
{
   FieldPattern::print(os);

   os << "FieldPattern: FieldAggPattern" << std::endl;
   os << "FieldPattern:    |numFields| = " << numFields_.size() << std::endl;
   os << "FieldPattern:    numFields = [ ";
   int total=0;
   for(std::size_t i=0;i<numFields_.size();i++)  {
      os << numFields_[i] << " ";
      total += numFields_[i]; 
   }
   os << "]" << std::endl;
   os << "FieldPattern:    |fieldIds| = " << fieldIds_.size() << " (" << total << ")" << std::endl;
   os << "FieldPattern:    fieldIds = [ ";
   for(std::size_t i=0;i<fieldIds_.size();i++) 
      os << fieldIds_[i] << " ";
   os << "]" << std::endl;
   os << "FieldPattern:    local offsets\n";

   std::map<int,int>::const_iterator itr;
   for(itr=fieldIdToPatternIdx_.begin();itr!=fieldIdToPatternIdx_.end();++itr) {
      int fieldId = itr->first;
      const std::vector<int> & offsets = localOffsets(fieldId);
      os << "FieldPattern:       field " << itr->first << " = [ ";
      for(std::size_t i=0;i<offsets.size();i++)
         os << offsets[i] << " ";
      os << "]" << std::endl;
   }
}

Teuchos::RCP<const FieldPattern> FieldAggPattern::getFieldPattern(int fieldId) const
{
   std::map<int,int>::const_iterator idxIter = fieldIdToPatternIdx_.find(fieldId);
   TEUCHOS_TEST_FOR_EXCEPTION(idxIter==fieldIdToPatternIdx_.end(),std::logic_error,
                     "FieldID = " << fieldId << " not defined in this pattern");

   return patterns_[idxIter->second].second;
}

void FieldAggPattern::buildFieldIdToPatternIdx()
{
   // convert map into vector
   int index = 0;
   std::vector<std::pair<int,FPPtr> >::const_iterator itr;
   for(itr=patterns_.begin();itr!=patterns_.end();++itr,++index)
      fieldIdToPatternIdx_[itr->first] = index;
}

void FieldAggPattern::buildFieldIdsVector()
{
   FPPtr geomAggPattern = getGeometricAggFieldPattern();

   // numFields_: stores number of fields per geometric ID
   // fieldIds_:  stores field IDs for each entry field indicated by numFields_
   numFields_.resize(geomAggPattern->numberIds());
   fieldIds_.clear();

   // build numFields_ and fieldIds_ vectors
   for(int dim=0;dim<getDimension()+1;dim++) {
      int numSubcell = geomAggPattern->getSubcellCount(dim);
      for(int sc=0;sc<numSubcell;sc++) {
         // merge the field patterns for multiple fields
         // on a specific dimension and subcell
         mergeFieldPatterns(dim,sc);
      }
   }
}

void FieldAggPattern::mergeFieldPatterns(int dim,int subcell)
{
   // make sure that the geometric IDs count is equal to the maxDOFs
   const std::vector<int> & geomIndices = getGeometricAggFieldPattern()->getSubcellIndices(dim,subcell);

   std::vector<std::pair<int,FPPtr> >::const_iterator itr;
   for(std::size_t i=0;i<geomIndices.size();i++) {
      numFields_[geomIndices[i]] = 0; // no fields initially      

      for(itr=patterns_.begin();itr!=patterns_.end();++itr) {
         std::size_t fieldSize = itr->second->getSubcellIndices(dim,subcell).size();

         // add field ID if their are enough in current pattern 
         if(i<fieldSize) { 
            fieldIds_.push_back(itr->first);
            numFields_[geomIndices[i]]++;
         }
      }

      TEUCHOS_ASSERT(numFields_[geomIndices[i]]>=0);
   }

}

void FieldAggPattern::buildFieldPatternData()
{
   int dimension = getDimension();
   FPPtr geomIdsPattern = getGeometricAggFieldPattern();

   // build patternData_ vector for implementation of the FieldPattern
   // functions (indicies will index into fieldIds_)
   patternData_.resize(dimension+1);
   int nextIndex = 0;
   for(int d=0;d<dimension+1;d++) {
      int numSubcell = geomIdsPattern->getSubcellCount(d);
      patternData_[d].resize(numSubcell);

      // loop through sub cells
      for(int sc=0;sc<numSubcell;sc++) {
         const std::vector<int> & geomIds = geomIdsPattern->getSubcellIndices(d,sc);
         for(std::size_t i=0;i<geomIds.size();i++) {
            int numFields = numFields_[geomIds[i]];
            for(int k=0;k<numFields;k++,nextIndex++) 
               patternData_[d][sc].push_back(nextIndex);
         }
      }
   }
}

const std::vector<int> & FieldAggPattern::localOffsets(int fieldId) const
{
   // lazy evaluation
   std::map<int,std::vector<int> >::const_iterator itr = fieldOffsets_.find(fieldId);
   if(itr!=fieldOffsets_.end())
      return itr->second;

   std::vector<int> & offsets = fieldOffsets_[fieldId];
   localOffsets_build(fieldId,offsets); 
   return offsets;
}

bool FieldAggPattern::LessThan::operator()(const Teuchos::Tuple<int,3> & a,const Teuchos::Tuple<int,3> & b) const 
{
   if(a[0] < b[0]) return true;
   if(a[0] > b[0]) return false;

   // a[0]==b[0]  
   if(a[1] < b[1]) return true;
   if(a[1] > b[1]) return false;

   // a[1]==b[1] && a[0]==b[0] 
   if(a[2] < b[2]) return true;
   if(a[2] > b[2]) return false;

   // a[2]==b[2] && a[1]==b[1] && a[0]==b[0]
   return false; // these are equal to, but not less than!
}

//const std::vector<int> & 
const std::pair<std::vector<int>,std::vector<int> > &
FieldAggPattern::localOffsets_closure(int fieldId,int subcellDim,int subcellId) const
{
   // lazy evaluation
   typedef std::map<Teuchos::Tuple<int,3>, std::pair<std::vector<int>,std::vector<int> >,LessThan> OffsetMap;

   Teuchos::Tuple<int,3> subcellTuple = Teuchos::tuple<int>(fieldId,subcellDim,subcellId);

   OffsetMap::const_iterator itr
         = fieldSubcellOffsets_closure_.find(subcellTuple);
   if(itr!=fieldSubcellOffsets_closure_.end()) {
      return itr->second;
   }

   TEUCHOS_TEST_FOR_EXCEPTION(subcellDim>=getDimension(),std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellDim<getDimension() failed");
   TEUCHOS_TEST_FOR_EXCEPTION(subcellId<0,std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellId>=0 failed");
   TEUCHOS_TEST_FOR_EXCEPTION(subcellId>=getSubcellCount(subcellDim),std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellId<getSubcellCount(subcellDim) failed");

   // build vector for sub cell closure indices
   ///////////////////////////////////////////////
  
   // grab field offsets
   const std::vector<int> & fieldOffsets = localOffsets(fieldId);

   // get offsets into field offsets for the closure indices (from field pattern)
   std::vector<int> closureOffsets; 
   FPPtr fieldPattern = getFieldPattern(fieldId);
   fieldPattern->getSubcellClosureIndices(subcellDim,subcellId,closureOffsets);

   // build closure indices into the correct location in lazy evaluation map.
   std::pair<std::vector<int>,std::vector<int> > & indicesPair
         = fieldSubcellOffsets_closure_[subcellTuple];

   std::vector<int> & closureIndices = indicesPair.first;
   for(std::size_t i=0;i<closureOffsets.size();i++)
      closureIndices.push_back(fieldOffsets[closureOffsets[i]]);

   std::vector<int> & basisIndices = indicesPair.second;
   basisIndices.assign(closureOffsets.begin(),closureOffsets.end());

   TEUCHOS_ASSERT(fieldSubcellOffsets_closure_[subcellTuple].first.size()==closureIndices.size());
   return fieldSubcellOffsets_closure_[subcellTuple];
}

void FieldAggPattern::localOffsets_build(int fieldId,std::vector<int> & offsets) const
{
   // This function makes the assumption that if there are multiple indices
   // shared by a subcell then the GeometricAggFieldPattern reflects this.
   // This is a fine assumption on edges and faces because the symmetries require
   // extra information about ordering. However, on nodes and "volumes" the
   // assumption appears to be stupid. For consistency we will make it.
   //
   // This code needs to be tested with a basis that has multple IDs at
   // a node or "volume" sub cell.

   FPPtr fieldPattern = getFieldPattern(fieldId);

   offsets.clear();
   offsets.resize(fieldPattern->numberIds(),-111111); // fill with some negative number
                                                     // for testing

   // this will offsets for all IDs associated with the field
   // but using a geometric ordering
   std::vector<int> fieldIdsGeomOrder;
   for(std::size_t i=0;i<fieldIds_.size();++i) {
      if(fieldIds_[i]==fieldId) 
         fieldIdsGeomOrder.push_back(i);
   }
   TEUCHOS_ASSERT((int) fieldIdsGeomOrder.size()==fieldPattern->numberIds());

   // built geometric ordering for this pattern: will index into fieldIdsGeomOrder
   GeometricAggFieldPattern geomPattern(fieldPattern);
   TEUCHOS_ASSERT((int) fieldIdsGeomOrder.size()==geomPattern.numberIds());
    
   for(int dim=0;dim<geomPattern.getDimension()+1;dim++) {
      for(int sc=0;sc<geomPattern.getSubcellCount(dim);sc++) {
         const std::vector<int> & gIndices = geomPattern.getSubcellIndices(dim,sc);
         const std::vector<int> & fIndices = fieldPattern->getSubcellIndices(dim,sc);

         TEUCHOS_ASSERT(gIndices.size()==fIndices.size());
         for(std::size_t i=0;i<gIndices.size();i++) { 
            offsets[fIndices[i]] = fieldIdsGeomOrder[gIndices[i]];
         }
      }
   }

   // check for failure/bug
   for(std::size_t i=0;i<offsets.size();i++) {
      TEUCHOS_ASSERT(offsets[i]>=0);
   }
}

}
