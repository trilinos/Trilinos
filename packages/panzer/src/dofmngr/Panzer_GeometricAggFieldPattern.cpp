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

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;

namespace panzer {

GeometricAggFieldPattern::GeometricAggFieldPattern()
   : patternBuilt_(false), dimension_(0)
{}

GeometricAggFieldPattern::GeometricAggFieldPattern(std::vector<Teuchos::RCP<const FieldPattern> > & patterns) 
   : patternBuilt_(false), dimension_(0)
{ 
   buildPattern(patterns); 
}

GeometricAggFieldPattern::GeometricAggFieldPattern(const Teuchos::RCP<const FieldPattern> & pattern) 
   : patternBuilt_(false), dimension_(0)
{ 
   buildPattern(pattern); 
}

void GeometricAggFieldPattern::buildPattern(const std::vector<Teuchos::RCP<const FieldPattern> > & patterns)
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
      sameGeometry &= patterns[0]->sameGeometry(*patterns[i]);       
   TEUCHOS_TEST_FOR_EXCEPTION(not sameGeometry,std::logic_error,
             "GeometricAggFieldPattern::buildPattern(): Patterns must "
             "have the same geometry!");

   // copy cell topology
   cellTopo_ = patterns[0]->getCellTopology();

   // grab the dimension
   dimension_ = patterns[0]->getDimension();
   patternData_.resize(dimension_+1);

   // build space for subcells
   std::vector<int> subcellCount(dimension_+1);
   for(std::size_t d=0;d<dimension_+1;d++) {
      subcellCount[d] = patterns[0]->getSubcellCount(d);
      patternData_[d].resize(subcellCount[d]);
   }

   // build geometric pattern: increment it logically
   // over all the subcells.
   int counter = 0;
   for(std::size_t d=0;d<dimension_+1;d++) {
      for(int s=0;s<subcellCount[d];s++) {
         std::vector<int> & current = patternData_[d][s];
         for(std::size_t p=0;p<patterns.size();p++) {
            RCP<const FieldPattern> field = patterns[p];
            std::size_t num = field->getSubcellIndices(d,s).size();
 
            if(current.size()<num) { 
               for(int i=num-current.size();i>0;i--,counter++) 
                  current.push_back(counter);
            }
         } 
      }
   }

   // record that the pattern has been built
   patternBuilt_ = true;
}

void GeometricAggFieldPattern::buildPattern(const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::vector<Teuchos::RCP<const FieldPattern> > patterns;
   patterns.push_back(pattern);
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
