// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_FieldAggPattern.hpp"

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {

FieldAggPattern::FieldAggPattern()
{
}

FieldAggPattern::FieldAggPattern(std::vector<std::tuple<int,panzer::FieldType,Teuchos::RCP<const FieldPattern> > > & patterns,
                                 const Teuchos::RCP<const FieldPattern> & geomAggPattern)
{
   buildPattern(patterns,geomAggPattern);
}

Teuchos::RCP<const FieldPattern> FieldAggPattern::getGeometricAggFieldPattern() const
{
   return geomAggPattern_;
}

void FieldAggPattern::buildPattern(const std::vector<std::tuple<int,panzer::FieldType,Teuchos::RCP<const FieldPattern> > > & patterns,
                                   const Teuchos::RCP<const FieldPattern> & geomAggPattern)
{
   TEUCHOS_ASSERT(patterns.size()>0);

   // build geometric information
   if(geomAggPattern==Teuchos::null) {
      std::vector<std::pair<FieldType,FPPtr>> patternVec;

      // convert map into vector
      auto itr = patterns.cbegin();
      for(;itr!=patterns.cend();++itr)
        patternVec.push_back(std::make_pair(std::get<1>(*itr),std::get<2>(*itr)));

      // build geometric aggregate field pattern
      geomAggPattern_ = rcp(new GeometricAggFieldPattern(patternVec));
   }
   else
      geomAggPattern_ = geomAggPattern;

   patterns_ = patterns;

   buildFieldIdToPatternIdx(); // builds look up from fieldId to index into patterns_ vector
   buildFieldIdsVector();      // meat of matter: build field pattern information
   buildFieldPatternData();    // this handles the "getSubcell*" information
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

   return std::get<2>(patterns_[idxIter->second]);
}

FieldType FieldAggPattern::getFieldType(int fieldId) const
{
   std::map<int,int>::const_iterator idxIter = fieldIdToPatternIdx_.find(fieldId);
   TEUCHOS_TEST_FOR_EXCEPTION(idxIter==fieldIdToPatternIdx_.end(),std::logic_error,
                     "FieldID = " << fieldId << " not defined in this pattern");

   return std::get<1>(patterns_[idxIter->second]);
}

void FieldAggPattern::buildFieldIdToPatternIdx()
{
   // convert map into vector
   int index = 0;
   auto itr = patterns_.cbegin();
   for(;itr!=patterns_.cend();++itr,++index)
     fieldIdToPatternIdx_[std::get<0>(*itr)] = index;
}

void FieldAggPattern::buildFieldIdsVector()
{
   FPPtr geomAggPattern = getGeometricAggFieldPattern();

   // numFields_: stores number of fields per geometric ID
   // fieldIds_:  stores field IDs for each entry field indicated by numFields_
   numFields_.resize(geomAggPattern->numberIds(),0);
   fieldIds_.clear();

   // build numFields_ and fieldIds_ vectors

   // Merge the field patterns for multiple fields for each dimension
   // and subcell. We do all the CG first, then all DG. This allows us
   // to use one offset for mapping DOFs to subcells later on.
   this->mergeFieldPatterns(FieldType::CG);
   this->mergeFieldPatterns(FieldType::DG);
}

void FieldAggPattern::mergeFieldPatterns(const FieldType& fieldType)
{
  auto geomAggPattern = getGeometricAggFieldPattern();

  // For DG, all DOFs are added to the internal cell
  const int cellDim = this->getDimension();
  const int numDimensions = getDimension()+1;

  for(int dim=0;dim<numDimensions;dim++) {
    int numSubcell = geomAggPattern->getSubcellCount(dim);
    for(int subcell=0;subcell<numSubcell;subcell++) {

      // Get the geometric index to add the DOF indices to. CG adds to
      // the subcell we are iterating over, (dim,subcel), while DG
      // adds to the internal cell (cellDim,0)
      const std::vector<int> * geomIndices = nullptr;
      if (fieldType == FieldType::CG)
        geomIndices = &(getGeometricAggFieldPattern()->getSubcellIndices(dim,subcell));
      else
        geomIndices = &(getGeometricAggFieldPattern()->getSubcellIndices(cellDim,0));

      if (geomIndices->size() > 0) {
        const int geomIndex = (*geomIndices)[0];

        auto itr = patterns_.cbegin();
        for(;itr!=patterns_.cend();++itr) {
          // Only merge in if pattern matches the FieldType.
          if (std::get<1>(*itr) == fieldType) {
            const std::size_t fieldSize = std::get<2>(*itr)->getSubcellIndices(dim,subcell).size();

            // add field ID if there are enough in current pattern
            for (std::size_t i=0;i<fieldSize;++i)
              fieldIds_.push_back(std::get<0>(*itr));
            numFields_[geomIndex]+=fieldSize;
          }
        }
        TEUCHOS_ASSERT(numFields_[geomIndex]>=0);
      }

    }
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
         // a single geometric entity is assigned to field pattern
         // geomIds size should be either 0 or 1.
         const std::vector<int> & geomIds = geomIdsPattern->getSubcellIndices(d,sc);
         TEUCHOS_ASSERT(geomIds.size()<=1);
         if (geomIds.size() > 0) {
           const int geomId = geomIds[0];
           const int numFields = numFields_[geomId];
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

const PHX::View<const int*> FieldAggPattern::localOffsetsKokkos(int fieldId) const
{
   // lazy evaluation
   std::map<int,PHX::View<int*> >::const_iterator itr = fieldOffsetsKokkos_.find(fieldId);
   if(itr!=fieldOffsetsKokkos_.end())
      return itr->second;

   const auto hostOffsetsStdVector = this->localOffsets(fieldId);
   PHX::View<int*> offsets("panzer::FieldAggPattern::localOffsetsKokkos",hostOffsetsStdVector.size());
   auto hostOffsets = Kokkos::create_mirror_view(offsets);
   for (size_t i=0; i < hostOffsetsStdVector.size(); ++i)
     hostOffsets(i) = hostOffsetsStdVector[i];
   Kokkos::deep_copy(offsets,hostOffsets);
   fieldOffsetsKokkos_[fieldId] = offsets;
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

   TEUCHOS_TEST_FOR_EXCEPTION(subcellDim >= getDimension(),std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellDim<getDimension() failed");
   TEUCHOS_TEST_FOR_EXCEPTION(subcellId < 0,std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellId>=0 failed");
   TEUCHOS_TEST_FOR_EXCEPTION(subcellId>=getSubcellCount(subcellDim),std::logic_error,
                         "FieldAggPattern::localOffsets_closure precondition subcellId<getSubcellCount(subcellDim) failed");

   // build vector for sub cell closure indices
   ///////////////////////////////////////////////
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
   FieldType fieldType = getFieldType(fieldId);

   offsets.clear();
   offsets.resize(fieldPattern->numberIds(),-111111); // fill with negative number for testing

   // this will offsets for all IDs associated with the field
   // but using a geometric ordering
   std::vector<int> fieldIdsGeomOrder;
   for(std::size_t i=0;i<fieldIds_.size();++i) {
      if(fieldIds_[i]==fieldId)
         fieldIdsGeomOrder.push_back(i);
   }

   // Check that number of DOFs line up
   TEUCHOS_ASSERT((int) fieldIdsGeomOrder.size()==fieldPattern->numberIds());

   // Build geometric ordering for this pattern: will index into fieldIdsGeomOrder
   GeometricAggFieldPattern geomPattern(fieldType,fieldPattern);
   int cnt = 0;
   for(int dim=0;dim<geomPattern.getDimension()+1;dim++) {
       for(int sc=0;sc<geomPattern.getSubcellCount(dim);sc++) {
          const std::vector<int> & fIndices = fieldPattern->getSubcellIndices(dim,sc);

          for(std::size_t i=0;i<fIndices.size();i++)
            offsets[fIndices[i]] = fieldIdsGeomOrder[cnt++];
       }
   }

   // check for failure/bug
   for(std::size_t i=0;i<offsets.size();i++) {
     TEUCHOS_ASSERT(offsets[i]>=0);
   }
}

}
