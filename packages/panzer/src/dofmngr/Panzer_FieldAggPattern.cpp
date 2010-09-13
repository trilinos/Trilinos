#include "Panzer_FieldAggPattern.hpp"

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {
 
FieldAggPattern::FieldAggPattern() 
{
}

FieldAggPattern::FieldAggPattern(std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns)
{
   buildPattern(patterns);
}

Teuchos::RCP<const FieldPattern> FieldAggPattern::getGeometricAggFieldPattern() const
{
   return geomAggPattern_;
}

void FieldAggPattern::buildPattern(const std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns)
{
   TEUCHOS_ASSERT(patterns.size()>0);

   // build geometric information
   {
      std::vector<FPPtr> patternVec;

      // convert map into vector
      std::vector<std::pair<int,FPPtr> >::const_iterator itr;
      for(itr=patterns.begin();itr!=patterns.end();++itr)
         patternVec.push_back(itr->second);

      // build geometric aggregate field pattern
      geomAggPattern_ = rcp(new GeometricAggFieldPattern(patternVec));
   }

   patterns_ = patterns;

   buildFieldIdToPatternIdx(); // builds look up from fieldId to index into patterns_ vector
   buildFieldIdsVector();      // meat of matter: build field pattern information
   buildFieldPatternData();    // this hanldes the "getSubcell*" information
}

int FieldAggPattern::getDimension() const
{
   FPPtr geomPattern = getGeometricAggFieldPattern();
   TEST_FOR_EXCEPTION(geomPattern==Teuchos::null,std::logic_error,
                      "Geometric field pattern not yet set, call buildPatterns first");

   return geomPattern->getDimension();
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
   TEST_FOR_EXCEPTION(true,std::logic_error,
                      "FieldAggPattern::getSubcellClousreIndices should not be called"); 
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
   TEST_FOR_EXCEPTION(idxIter==fieldIdToPatternIdx_.end(),std::logic_error,
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

      TEUCHOS_ASSERT(numFields_[geomIndices[i]]>0);
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
