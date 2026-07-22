// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;

namespace panzer {

GeometricAggFieldPattern::GeometricAggFieldPattern()
   : patternBuilt_(false), dimension_(0)
{}

  GeometricAggFieldPattern::GeometricAggFieldPattern(std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> & patterns) 
   : patternBuilt_(false), dimension_(0)
{ 
   buildPattern(patterns); 
}

GeometricAggFieldPattern::GeometricAggFieldPattern(const FieldType& fieldType,
                                                   const Teuchos::RCP<const FieldPattern> & pattern) 
   : patternBuilt_(false), dimension_(0)
{ 
  buildPattern(fieldType,pattern);
}

void GeometricAggFieldPattern::buildPattern(const std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> & patterns)
{
   std::size_t numPat = patterns.size();

   // must be at least one field to do something
   if(numPat<1) {
      bool no_patterns_to_construct = true;
      TEUCHOS_TEST_FOR_EXCEPTION(no_patterns_to_construct,std::logic_error,
                         "GeometricAggFieldPattern::buildPattern requires at least one field pattern");
      return;
   }

   bool sameGeometry=true;
   for(std::size_t i=1;i<patterns.size();i++)
     sameGeometry &= patterns[0].second->sameGeometry(*(patterns[i].second));       
   TEUCHOS_TEST_FOR_EXCEPTION(not sameGeometry,std::logic_error,
             "GeometricAggFieldPattern::buildPattern(): Patterns must "
             "have the same geometry!");

   // copy cell topology
   cellTopo_ = patterns[0].second->getCellTopology();

   // grab the dimension
   dimension_ = patterns[0].second->getDimension();
   patternData_.resize(dimension_+1);

   // build space for subcells
   std::vector<int> subcellCount(dimension_+1);
   for(std::size_t d=0;d<dimension_+1;d++) {
      subcellCount[d] = patterns[0].second->getSubcellCount(d);
      patternData_[d].resize(subcellCount[d]);
   }

   // Build geometric pattern: increment it logically over all the
   // subcells. Start with CG fields first. Then all DG fields. This
   // is done so that we can use individual field pattern offsets when
   // mapping DOFs to subcells.
   int counter = 0;
   for(std::size_t d=0;d<dimension_+1;d++) {
      for(int s=0;s<subcellCount[d];s++) {
         std::vector<int> & current = patternData_[d][s];
         for(std::size_t p=0;p<patterns.size();p++) {
           if ( (patterns[p].first) == FieldType::CG) {
              RCP<const FieldPattern> field = (patterns[p]).second;
              // if dofs exist, we have a geometric entity 
              const std::size_t num = ( (field->getSubcellIndices(d,s).size() > 0) ? 1 : 0 );
              if(current.size()<num) { 
                for(int i=num-current.size();i>0;i--,counter++) 
                  current.push_back(counter);
              }
            }
         } 
      }
   }

   // Add DG fields. These fields are considered internal "cell"
   // fields for DOF numbering.
   const int cellDim = dimension_;
   for(std::size_t d=0;d<dimension_+1;d++) {
      for(int s=0;s<subcellCount[d];s++) {
         std::vector<int> & current = patternData_[cellDim][0];
         for(std::size_t p=0;p<patterns.size();p++) {
           if ( (patterns[p].first) == FieldType::DG) {
              RCP<const FieldPattern> field = (patterns[p]).second;
              // if dofs exist, we have a geometric entity
              const std::size_t num = ( (field->getSubcellIndices(d,s).size() > 0) ? 1 : 0 );
              if(current.size()<num) { 
                for(int i=num-current.size();i>0;i--,counter++) 
                  current.push_back(counter);
              }
            }
         } 
      }
   }
   
   // record that the pattern has been built
   patternBuilt_ = true;
}

void GeometricAggFieldPattern::buildPattern(const FieldType& fieldType,
                                            const Teuchos::RCP<const FieldPattern> & pattern)
{
  std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> patterns;
  patterns.push_back(std::make_pair(fieldType,pattern));
  buildPattern(patterns);
}

int GeometricAggFieldPattern::getSubcellCount(int dim) const
{
   if(patternBuilt_) return patternData_[dim].size();

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getSubcellCount() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

const std::vector<int> & GeometricAggFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(patternBuilt_) return patternData_[dim][cellIndex];

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getSubcellIndices() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

int GeometricAggFieldPattern::getDimension() const
{
   if(patternBuilt_) return dimension_;

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getDimension() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

shards::CellTopology GeometricAggFieldPattern::getCellTopology() const
{
   if(patternBuilt_) return cellTopo_;

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getCellTopology() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

}
